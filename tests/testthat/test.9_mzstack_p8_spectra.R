context("mzStack: spectra read and written through Spectra")

.lib_file <- function(f) system.file("extdata", "tests", "mzstack", "library",
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

.sorted <- function(x) x[order(as.numeric(names(x)))]

test_that("peaks read through MsBackendParquet equal the arrow reader", {
    skip_if_no_mzstack()
    testthat::skip_if_not_installed("MsBackendParquet")
    lib <- tempfile("lib-")
    suppressMessages(convertLibraryToMzstack(.lib_file("mini_library.sqlite"),
                                             lib))
    flags <- c("x_mspurity_snr_pass_flag", "x_mspurity_minnum_pass_flag",
               "x_mspurity_minfrac_pass_flag", "x_mspurity_ra_pass_flag")
    for (p in c(lib, .results_dataset())) {
        m <- .mzs(".mzs_read_manifest")(p)
        ids <- .mzs(".mzs_read_spectra_meta")(p, m, "id")$spectrum_id_
        expect_gt(length(ids), 0)
        a <- .mzs(".mzs_read_peaks_spectra")(p, ids, flags)
        b <- .mzs(".mzs_read_peaks_arrow")(p, m, ids, flags)
        expect_identical(.sorted(a), .sorted(b))
    }
})

test_that("a Spectra library converts like the MSP file it was read from", {
    skip_if_no_mzstack()
    a <- tempfile("lib-")
    suppressMessages(convertLibraryToMzstack(.lib_file("mini_mona.msp"), a,
                                             format = "msp", dialect = "mona"))
    b <- tempfile("lib-")
    suppressMessages(convertLibraryToMzstack(.mona_spectra(), b))
    read <- function(p) {
        m <- .manifest(p)
        as.data.frame(dplyr::bind_rows(lapply(m$runs, function(r) as.data.frame(
            arrow::read_parquet(file.path(p, r$path, "part-0.parquet")))[, c(
                "id", "ms_level", "scan_polarity", "selected_ion_mz", "mz",
                "intensity", "x_mspurity_precursor_type",
                "x_mspurity_inchikey", "x_mspurity_formula")])))
    }
    sa <- read(a)
    sb <- read(b)
    k <- match(sa$id, sb$id)
    expect_false(anyNA(k))
    for (col in setdiff(names(sa), "id"))
        expect_identical(sb[[col]][k], sa[[col]], info = col)
    v <- validateMzstack(b)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("matching against a Spectra library scores as against the MSP file", {
    skip_if_no_mzstack()
    q <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite",
                     package = "msPurity")
    qd <- tempfile("q-")
    suppressMessages(convertSqliteToMzstack(q, qd))
    a <- tempfile("lib-")
    suppressMessages(convertLibraryToMzstack(.lib_file("mini_mona.msp"), a,
                                             format = "msp", dialect = "mona"))
    b <- tempfile("lib-")
    suppressMessages(convertLibraryToMzstack(.mona_spectra(), b))
    ra <- suppressMessages(spectralMatching(qd, a, q_spectraTypes = "av_all",
                                            cores = 1, format = "mzstack"))
    rb <- suppressMessages(spectralMatching(qd, b, q_spectraTypes = "av_all",
                                            cores = 1, format = "mzstack"))
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
    skip_if_no_mzstack()
    sp <- .mona_spectra()
    expect_error(convertLibraryToMzstack(sp, tempfile(), format = "msp"),
                 "Spectra object")
    expect_error(convertLibraryToMzstack(.lib_file("mini_mona.msp"),
                                         tempfile(), format = "spectra"),
                 "needs a Spectra object")
    expect_error(convertLibraryToMzstack(sp, tempfile(),
                                         mapping = c(not_a_field = "name")),
                 "does not have")
})
