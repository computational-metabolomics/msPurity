.alignment_peaks <- function(mz, ra, w) {
    data.frame(mz = mz, ra = ra, w = w)
}

.alignment <- .mzs("align")

test_that("intensity ties select eligible peaks and preserve later matches", {
    q <- .alignment_peaks(c(100.1234, 150.0987), c(100, 100), c(11, 12))
    l <- .alignment_peaks(c(150.0987, 100.1234), c(100, 100), c(21, 22))
    x <- .alignment(q, l, q_ppmProd = 10, l_ppmProd = 10)
    expect_identical(x$q, c(11, 12))
    expect_identical(x$l, c(22, 21))
    expect_identical(x$mcount, 2L)
    expect_identical(x$total, 2L)
    expect_identical(x$percMatchAll, 1)
})

test_that("ppm fallback cannot select an equal-distance ineligible peak", {
    # Deliberately wide tolerance gives equal ppm errors with only one overlap.
    q <- .alignment_peaks(150.1875, 50, 1)
    l <- .alignment_peaks(c(300.3750, 100.1250), c(20, 20), c(2, 3))
    x <- .alignment(q, l, q_ppmProd = 400000, l_ppmProd = 0)
    expect_identical(x$q, c(1, 0))
    expect_identical(x$l, c(3, 2))
    expect_identical(x$mcount, 1L)
})

test_that("eligible ties keep input order and consumed peaks are not reused", {
    q <- .alignment_peaks(rep(100.1234, 3), rep(100, 3), 11:13)
    l <- .alignment_peaks(c(100.1234, 100.1234), c(100, 100), c(21, 22))
    x <- .alignment(q, l)
    expect_identical(x$q, c(11, 12, 13))
    expect_identical(x$l, c(21, 22, 0))
    expect_identical(x$mcount, 2L)
    expect_identical(x$total, 3L)
})

test_that("intensity selection and nearest-ppm fallback keep their priority", {
    q <- .alignment_peaks(100.1234, 100, 11)
    l <- .alignment_peaks(c(100.1235, 100.1236), c(90, 99), c(21, 22))
    expect_identical(.alignment(q, l)$l, c(22, 21))
    l$ra <- c(50, 50)
    expect_identical(.alignment(q, l)$l, c(21, 22))
    l$ra <- c(50, 90)
    # A difference exactly equal to the intensity threshold uses ppm fallback.
    expect_identical(.alignment(q, l, raDiffThres = 10)$l, c(21, 22))
    l$ra <- c(NA_real_, 99)
    expect_identical(.alignment(q, l)$l, c(22, 21))
    expect_identical(.alignment(q, l, raDiffThres = NA_real_)$l, c(21, 22))
})

test_that("rounded ppm boundaries do not admit outside ties", {
    q <- .alignment_peaks(100.1234, 100, 11)
    for (mz in list(c(100.1223987658, 100.1223987660),
            c(100.1244012342, 100.1244012340))) {
        l <- .alignment_peaks(mz, c(100, 100), c(21, 22))
        x <- .alignment(q, l, q_ppmProd = 10, l_ppmProd = 0)
        case <- paste("reference m/z:",
            paste(format(mz, digits = 12), collapse = ", "))
        expect_identical(x$l, c(22, 21), info = case)
        expect_identical(x$mcount, 1L, info = case)
    }
})

test_that("non-overlapping peaks remain unmatched", {
    q <- .alignment_peaks(100.1234, 100, 11)
    l <- .alignment_peaks(50.0617, 100, 21)
    x <- .alignment(q, l)
    expect_identical(x$q, c(11, 0))
    expect_identical(x$l, c(0, 21))
    expect_identical(x$mcount, 0L)
    expect_identical(x$percMatchAll, 0)
})

test_that("the deprecated matcher restricts ties to its own eligible peaks", {
    match_peaks <- .mzs("matchi")
    q <- cbind(mz = 100.1234, ra = 100, type = 1, w = 11)
    l <- cbind(mz = c(50.0617, 100.1234), ra = c(100, 100), type = 2,
        w = c(21, 22))
    x <- match_peaks(l, q, ppmdiff = matrix(c(1e6, 0), 1),
        idiff = matrix(c(0, 0), 1), ppm_tol_prod = 10, match_alg = "dpc")
    expect_equal(x, c(22 / sqrt(21^2 + 22^2), 0.5, 1), tolerance = 1e-14)
})

test_that("the deprecated matcher keeps strict ppm boundaries and exhaustion", {
    match_peaks <- .mzs("matchi")
    q <- cbind(mz = c(100.1234, 150.0987), ra = c(100, 100), type = 1,
        w = c(11, 12))
    l <- cbind(mz = c(100.1234, 50.0617), ra = c(100, 100), type = 2,
        w = c(21, 22))
    x <- match_peaks(l, q, ppmdiff = rbind(c(10, 100), c(0, 100)),
        idiff = matrix(0, 2, 2), ppm_tol_prod = 10, match_alg = "dpc")
    expect_equal(x, c(12 * 21 / sqrt((11^2 + 12^2) * (21^2 + 22^2)),
            1 / 3, 1), tolerance = 1e-14)
    x <- match_peaks(l, q, ppmdiff = rbind(c(0, 100), c(0, 100)),
        idiff = matrix(0, 2, 2), ppm_tol_prod = 10, match_alg = "dpc")
    expect_equal(x, c(11 * 21 / sqrt((11^2 + 12^2) * (21^2 + 22^2)),
            1 / 3, 1), tolerance = 1e-14)
})

.alignment_library <- function(mz) {
    path <- tempfile(fileext = ".sqlite")
    original <- system.file("extdata", "tests", "library",
        "mini_library.sqlite", package = "msPurity")
    stopifnot(file.copy(original, path))
    con <- DBI::dbConnect(RSQLite::SQLite(), path)
    on.exit(DBI::dbDisconnect(con))
    meta <- DBI::dbReadTable(con, "library_spectra_meta")
    pid <- min(meta$id)
    DBI::dbExecute(con, "DELETE FROM library_spectra_annotation")
    DBI::dbExecute(con, "DELETE FROM library_spectra")
    DBI::dbExecute(con, "DELETE FROM library_spectra_meta WHERE id != ?",
        params = list(pid))
    DBI::dbExecute(con,
        "UPDATE library_spectra_meta SET precursor_mz = 200.2468, name = ?",
        params = list("Synthetic peak-alignment fixture"))
    peaks <- data.frame(id = seq_along(mz), mz = mz, i = 10, other = "",
        library_spectra_meta_id = pid)
    DBI::dbWriteTable(con, "library_spectra", peaks, append = TRUE)
    path
}

.alignment_query <- function() {
    path <- tempfile(fileext = ".sqlite")
    original <- system.file("extdata", "tests", "db",
        "createDatabase_example.sqlite", package = "msPurity")
    stopifnot(file.copy(original, path))
    con <- DBI::dbConnect(RSQLite::SQLite(), path)
    on.exit(DBI::dbDisconnect(con))
    # Use one existing all-member average; preserve its feature associations.
    meta <- DBI::dbGetQuery(con, paste("SELECT pid FROM s_peak_meta",
            "WHERE grpid = '89' AND spectrum_type = 'all'"))
    stopifnot(nrow(meta) == 1L)
    peak <- DBI::dbGetQuery(con,
        "SELECT * FROM s_peaks WHERE pid = ? AND type = 'all' LIMIT 1",
        params = list(meta$pid))
    stopifnot(nrow(peak) == 1L)
    peak$mz <- 100.1234
    peak$i <- 10
    peak$ra <- 100
    DBI::dbExecute(con, "DELETE FROM s_peaks WHERE pid = ?",
        params = list(meta$pid))
    DBI::dbWriteTable(con, "s_peaks", peak, append = TRUE)
    DBI::dbExecute(con,
        "UPDATE s_peak_meta SET precursor_mz = 200.2468 WHERE pid = ?",
        params = list(meta$pid))
    path
}

.alignment_match <- function(q, l, format = "sqlite") {
    suppressMessages(spectralMatching(q, l, format = format,
            q_spectraTypes = "av_all", q_xcmsGroups = 89, l_spectraTypes = NA,
            q_spectraFilter = FALSE, l_spectraFilter = FALSE,
            q_pol = NA, l_pol = NA, q_ppmProd = 10, l_ppmProd = 10,
            cores = 1))$matchedResults
}

test_that("SQLite matching scores the eligible peak independently of Arrow", {
    q <- .alignment_query()
    l <- .alignment_library(c(50.0617, 100.1234))
    a <- .alignment_match(q, l)
    expect_identical(nrow(a), 1L)
    # Equal intensities and the 2:1 mass ratio give 4:1 score weights.
    expected <- c(dpc = 4 / sqrt(17), rdpc = 4 / sqrt(17),
        cdpc = 8 / (3 * sqrt(17)), mcount = 1, allcount = 2, mpercent = 0.5)
    for (s in names(expected))
        expect_equal(as.numeric(a[[s]]), unname(expected[s]),
            tolerance = 1e-14, info = s)
})

test_that("SQLite and Parquet routes score the same eligible peaks", {
    skip_if_no_arrow()
    q <- .alignment_query()
    l <- .alignment_library(c(50.0617, 100.1234))
    qp <- tempfile("query-")
    lp <- tempfile("reference-")
    .without_identity_warning(convertSqliteToParquet(q, qp))
    convertLibraryToParquet(l, lp)
    a <- .alignment_match(q, l)
    b <- .alignment_match(qp, lp, "parquet")
    expect_identical(nrow(a), 1L)
    expect_identical(nrow(b), 1L)
    for (s in c("dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent"))
        expect_identical(as.numeric(a[[s]]), as.numeric(b[[s]]), info = s)
})
