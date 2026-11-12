context("Parquet: spectra read and written through Spectra")

.lib_file <- function(f) system.file("extdata", "tests", "library",
                                     f, package = "msPurity")

.mona_spectra <- function() {
    testthat::skip_if_not_installed("MsBackendMsp")
    Spectra::Spectra(.lib_file("mini_mona.msp"),
                     source = MsBackendMsp::MsBackendMsp(),
                     mapping = c(name = "Name", accession = "DB#",
                                 precursorMz = "PrecursorMZ",
                                 adduct = "Precursor_type",
                                 inchikey = "InChIKey", formula = "Formula",
                                 polarity = "Ion_mode"))
}

# Every run's stored rows, in manifest run order.
.stored <- function(p) {
    m <- MsBackendParquet::readManifest(p)
    as.data.frame(dplyr::bind_rows(lapply(m$runs, function(r) {
        sp <- .mzs(".mzs_read_run")(p, r)
        sp$run_id <- rep(r$run_id, nrow(sp))
        sp
    })))
}

test_that("spectra read through MsBackendParquet are the stored ones", {
    skip_if_no_arrow()
    lib <- tempfile("lib-")
    suppressMessages(convertLibraryToParquet(.lib_file("mini_library.sqlite"),
                                             lib))
    flags <- c("x_mspurity_snr_pass_flag", "x_mspurity_minnum_pass_flag",
               "x_mspurity_minfrac_pass_flag", "x_mspurity_ra_pass_flag")
    for (p in c(lib, .results_dataset())) {
        st <- .stored(p)
        expect_gt(nrow(st), 0)
        cols <- c("id", "data_origin", "ms_level", "scan_polarity",
                  "selected_ion_mz")
        meta <- .mzs(".mzs_read_spectra_meta")(
            p, MsBackendParquet::readManifest(p), cols)
        expect_identical(meta$run_id, st$run_id)
        expect_equal(meta$spectrum_id_, st$spectrum_id_)
        ## A core variable with no values is not stored, and reads as NA.
        for (c in cols) {
            if (is.null(st[[c]]))
                expect_true(all(is.na(meta[[c]])), info = c)
            else
                expect_equal(meta[[c]], st[[c]], info = c)
        }

        have <- intersect(flags, names(st))
        pk <- .mzs(".mzs_read_peaks")(p, rev(st$spectrum_id_), flags)
        expect_identical(names(pk), as.character(st$spectrum_id_))
        want <- lapply(seq_len(nrow(st)), function(k) {
            d <- data.frame(mz = st$mz[[k]], i = st$intensity[[k]])
            for (f in have) {
                v <- st[[f]][[k]]
                d[[f]] <- if (is.null(v)) rep(NA, nrow(d)) else v
            }
            d
        })
        expect_identical(unname(pk), want)
    }
})

test_that("a Spectra library converts like the MSP file it was read from", {
    skip_if_no_arrow()
    a <- tempfile("lib-")
    suppressMessages(convertLibraryToParquet(.lib_file("mini_mona.msp"), a,
                                             format = "msp", dialect = "mona"))
    b <- tempfile("lib-")
    suppressMessages(convertLibraryToParquet(.mona_spectra(), b))
    read <- function(p) .stored(p)[, c(
        "id", "ms_level", "scan_polarity", "selected_ion_mz", "mz",
        "intensity", "x_mspurity_precursor_type", "x_mspurity_inchikey",
        "x_mspurity_formula")]
    sa <- read(a)
    sb <- read(b)
    k <- match(sa$id, sb$id)
    expect_false(anyNA(k))
    for (col in setdiff(names(sa), "id"))
        expect_identical(sb[[col]][k], sa[[col]], info = col)
    v <- validateParquet(b)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("matching against a Spectra library scores as against the MSP file", {
    skip_if_no_arrow()
    q <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite",
                     package = "msPurity")
    qd <- tempfile("q-")
    .without_identity_warning(
        suppressMessages(convertSqliteToParquet(q, qd)))
    a <- tempfile("lib-")
    suppressMessages(convertLibraryToParquet(.lib_file("mini_mona.msp"), a,
                                             format = "msp", dialect = "mona"))
    b <- tempfile("lib-")
    suppressMessages(convertLibraryToParquet(.mona_spectra(), b))
    ra <- suppressMessages(spectralMatching(qd, a, q_spectraTypes = "av_all",
                                            cores = 1, format = "parquet"))
    rb <- suppressMessages(spectralMatching(qd, b, q_spectraTypes = "av_all",
                                            cores = 1, format = "parquet"))
    ma <- ra$matchedResults
    mb <- rb$matchedResults
    expect_true(nrow(ma) > 0)
    ka <- paste(ma$qpid, ma$library_accession)
    kb <- paste(mb$qpid, mb$library_accession)
    expect_setequal(kb, ka)
    k <- match(ka, kb)
    for (s in c("dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent"))
        expect_identical(mb[[s]][k], ma[[s]], info = s)
})

test_that("a Spectra library needs format = \"spectra\"", {
    skip_if_no_arrow()
    sp <- .mona_spectra()
    expect_error(convertLibraryToParquet(sp, tempfile(), format = "msp"),
                 "Spectra object")
    expect_error(convertLibraryToParquet(.lib_file("mini_mona.msp"),
                                         tempfile(), format = "spectra"),
                 "needs a Spectra object")
    expect_error(convertLibraryToParquet(sp, tempfile(),
                                         mapping = c(not_a_field = "name")),
                 "does not have")
})
