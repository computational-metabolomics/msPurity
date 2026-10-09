context("Parquet: library datasets from msp2db and MSP")

.lib_file <- function(f) system.file("extdata", "tests", "library",
                                     f, package = "msPurity")

.convert_lib <- function(x, ...) {
    skip_if_no_arrow()
    p <- tempfile("lib-")
    convertLibraryToParquet(x, p, ...)
    p
}

.lib_spectra <- function(p) {
    m <- .manifest(p)
    out <- lapply(m$runs, function(r) {
        sp <- as.data.frame(arrow::read_parquet(file.path(p, r$path,
                                                          "part-0.parquet")))
        sp$run_id <- r$run_id
        sp
    })
    as.data.frame(dplyr::bind_rows(out))
}

test_that("an msp2db library becomes a library dataset (R-005, R-006)", {
    p <- .convert_lib(.lib_file("mini_library.sqlite"))
    m <- .manifest(p)
    expect_identical(m$role, "library")
    expect_match(m$library$digest, "^sha256:[0-9a-f]{64}$")
    expect_identical(m$library$digest, paste0("sha256:", unname(
        tools::sha256sum(.lib_file("mini_library.sqlite")))))
    expect_false(is.null(m$library$version))
    ## One run per library source.
    expect_identical(length(m$runs), 2L)
    expect_identical(sort(unname(unlist(m$provenance[[1]]$parameters$run_names))),
                     c("gnps", "massbank"))
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character())
    expect_identical(.mzs(".mzs_schema_check")(m, "manifest.schema.json"),
                     character())
})

test_that("within a run, spectra are ordered by polarity and precursor m/z", {
    sp <- .lib_spectra(.convert_lib(.lib_file("mini_library.sqlite")))
    for (r in unique(sp$run_id)) {
        s <- sp[sp$run_id == r, ]
        expect_identical(order(s$scan_polarity, s$selected_ion_mz),
                         seq_len(nrow(s)))
    }
})

test_that("peaks keep their source order and values exactly", {
    sp <- .lib_spectra(.convert_lib(.lib_file("mini_library.sqlite")))
    con <- DBI::dbConnect(RSQLite::SQLite(), .lib_file("mini_library.sqlite"))
    on.exit(DBI::dbDisconnect(con))
    meta <- DBI::dbGetQuery(con, "SELECT id, accession FROM library_spectra_meta")
    pk <- DBI::dbGetQuery(con, "SELECT * FROM library_spectra ORDER BY id")
    for (i in seq_len(nrow(sp))) {
        id <- meta$id[meta$accession == sp$id[i]]
        expect_identical(sp$mz[[i]], pk$mz[pk$library_spectra_meta_id == id])
        expect_identical(sp$intensity[[i]], pk$i[pk$library_spectra_meta_id == id])
    }
    ## MSP: the values are the file's own literals, in file order.
    lines <- readLines(.lib_file("mini_mona.msp"))
    first <- which(startsWith(lines, "Num Peaks"))[1]
    n <- as.integer(sub("^.*: ", "", lines[first]))
    tok <- do.call(rbind, strsplit(lines[first + seq_len(n)], " "))
    msp <- .lib_spectra(.convert_lib(.lib_file("mini_mona.msp"),
                                     format = "msp", dialect = "mona"))
    acc <- sub("^DB#: ", "", lines[which(startsWith(lines, "DB#"))[1]])
    expect_identical(msp$mz[[match(acc, msp$id)]], as.numeric(tok[, 1]))
    expect_identical(msp$intensity[[match(acc, msp$id)]], as.numeric(tok[, 2]))
})

test_that("msp2db, MassBank MSP and MoNA MSP give the same spectra (B-057)", {
    a <- .lib_spectra(.convert_lib(.lib_file("mini_library.sqlite")))
    b <- .lib_spectra(.convert_lib(.lib_file("mini_massbank.msp"),
                                   format = "msp", dialect = "massbank"))
    c <- .lib_spectra(.convert_lib(.lib_file("mini_mona.msp"),
                                   format = "msp", dialect = "mona"))
    cols <- c("ms_level", "scan_polarity", "selected_ion_mz",
              "collision_energy", "mz", "intensity",
              "x_mspurity_precursor_type", "x_mspurity_inchikey",
              "x_mspurity_collision_energy_text", "x_mspurity_formula")
    for (d in list(b, c)) {
        k <- match(a$id, d$id)
        expect_false(anyNA(k))
        for (col in cols)
            expect_identical(d[[col]][k], a[[col]], info = col)
    }
})

test_that("an MSP file round-trips stably (parse idempotence)", {
    for (dialect in c("massbank", "mona")) {
        p1 <- .convert_lib(.lib_file(paste0("mini_", dialect, ".msp")),
                           format = "msp", dialect = dialect)
        f <- tempfile(fileext = ".msp")
        .mzs(".mzs_write_msp")(p1, f, dialect)
        p2 <- .convert_lib(f, format = "msp", dialect = dialect)
        f2 <- tempfile(fileext = ".msp")
        .mzs(".mzs_write_msp")(p2, f2, dialect)
        p3 <- .convert_lib(f2, format = "msp", dialect = dialect)
        a <- .lib_spectra(p2)
        b <- .lib_spectra(p3)
        a$data_origin <- b$data_origin <- a$run_id <- b$run_id <- NULL
        expect_identical(b, a, info = dialect)
        expect_identical(readLines(f2), readLines(f), info = dialect)
    }
})

test_that("MSP keys outside the declared dialect are refused, not skipped", {
    skip_if_no_arrow()
    expect_parquet_error(convertLibraryToParquet(
        .lib_file("mini_mona.msp"), tempfile(), format = "msp"),
        "unsupported")
    ## MoNA keys are not MassBank keys.
    e <- tryCatch(convertLibraryToParquet(
        .lib_file("mini_mona.msp"), tempfile(), format = "msp",
        dialect = "massbank"), error = function(e) e)
    expect_s3_class(e, "msPurity_parquet_unsupported")
    expect_identical(e$line, 1L)
    f <- tempfile(fileext = ".msp")
    writeLines(c(readLines(.lib_file("mini_mona.msp"))[1:3],
                 "Mystery_key: 1", "Num Peaks: 1", "100 1", ""), f)
    e <- tryCatch(convertLibraryToParquet(f, tempfile(), format = "msp",
                                          dialect = "mona"),
                  error = function(e) e)
    expect_s3_class(e, "msPurity_parquet_unsupported")
    expect_identical(e$constructs, "MYSTERY_KEY")
    expect_identical(e$line, 4L)
    ## A peak count that disagrees with the peaks listed.
    writeLines(c("Name: x", "PrecursorMZ: 100", "Num Peaks: 2", "100 1", ""), f)
    expect_parquet_error(convertLibraryToParquet(f, tempfile(), format = "msp",
                                                 dialect = "mona"), "semantic")
})

test_that("collision energy is numeric only where unambiguous", {
    ce <- .mzs(".mzs_collision_energy")
    expect_identical(ce(c("35", "35 eV", "5.0eV", "35%", "HCD 35 (NCE)",
                          "ramp 20-40", "Ramp 5-60 V", NA)),
                     c(35, 35, 5, NA, NA, NA, NA, NA))
})

test_that("retention time is converted only when its unit is declared", {
    msp <- system.file("extdata", "tests", "msp", "av_all.msp",
                       package = "msPurity")
    a <- .lib_spectra(.convert_lib(msp, format = "msp", dialect = "mspurity"))
    expect_false("time" %in% names(a))
    expect_false(anyNA(a$x_mspurity_retention_time))
    b <- .lib_spectra(.convert_lib(msp, format = "msp", dialect = "mspurity",
                                   rtUnit = "s"))
    expect_identical(b$time, b$x_mspurity_retention_time / 60)
})

test_that("every createMSP() file converts in the mspurity dialect", {
    for (f in list.files(system.file("extdata", "tests", "msp",
                                     package = "msPurity"), full.names = TRUE)) {
        p <- .convert_lib(f, format = "msp", dialect = "mspurity")
        v <- validateParquet(p)
        expect_identical(v$message[v$level == "MUST"], character(),
                         info = basename(f))
    }
})

test_that("matching against a converted library scores as SQLite does", {
    skip_if_no_arrow()
    lib <- .lib_file("mini_library.sqlite")
    q <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite",
                     package = "msPurity")
    qd <- tempfile("q-")
    convertSqliteToParquet(q, qd)
    ld <- .convert_lib(lib)
    for (types in list("av_all", c("av_all", "inter"))) {
        a <- suppressMessages(spectralMatching(q, lib, q_spectraTypes = types,
                                               cores = 1))$matchedResults
        b <- suppressMessages(spectralMatching(qd, ld, q_spectraTypes = types,
                                               cores = 1,
                                               format = "parquet"))$matchedResults
        si <- .read_table(qd, "source_identifier")
        sp <- si[si$target_table == "spectra", ]
        b$qpid <- as.numeric(sp$source_value[match(b$qpid, sp$target_key)])
        ka <- paste(a$qpid, a$library_accession)
        kb <- paste(b$qpid, b$library_accession)
        expect_true(nrow(a) > 0)
        expect_setequal(kb, ka)
        k <- match(ka, kb)
        for (s in c("dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent"))
            expect_identical(as.numeric(b[[s]][k]), as.numeric(a[[s]]),
                             info = s)
    }
    ## And the evidence references the library by its release.
    suppressMessages(spectralMatching(qd, ld, q_spectraTypes = "av_all",
                                      cores = 1, format = "parquet",
                                      updateDb = TRUE))
    m <- .manifest(qd)
    lib_src <- Filter(function(s) identical(s$role, "library"), m$sources)
    expect_identical(lib_src[[1]]$digest, .manifest(ld)$library$digest)
    ev <- .read_table(qd, "evidence")
    expect_true(all(ev$reference_source == lib_src[[1]]$key))
    expect_false(anyNA(ev$reference_spectrum_id_))
    v <- validateParquet(qd)
    expect_identical(v$message[v$level == "MUST"], character())
})
