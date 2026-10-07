context("SQLite output: frozen, and deprecated with a warning")

# The SQLite writers must keep producing exactly what msPurity 1.37.3 wrote.
# The golden dumps were generated from master before the mzStack backend was
# added; see fixtures/make-sqlite-golden.R.

for (nm in names(.golden_writers())) local({
    nm <- nm
    test_that(paste("SQLite output is unchanged:", nm), {
        skip_if_no_mspuritydata()
        skip_if_not_installed("xcms")
        suppressPackageStartupMessages(library(xcms))
        golden <- file.path(.golden_dir(), paste0(nm, ".rds"))
        skip_if_not(file.exists(golden), "no golden dump")
        expect_identical(.golden_run(nm), readRDS(golden))
    })
})

# Each SQLite-writing entry point raises exactly one deprecatedWarning per
# call, naming the alternative.
.deprecation_count <- function(writer) {
    td <- tempfile("dep-")
    dir.create(td)
    res <- .with_deprecations(suppressMessages(
        utils::capture.output(writer(td))))
    res
}

test_that("legacy SQLite writers warn exactly once", {
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    suppressPackageStartupMessages(library(xcms))
    w <- .golden_writers()
    for (nm in names(w)) {
        r <- .deprecation_count(w[[nm]])
        expect_identical(r$n, 1L, info = nm)
        expect_match(r$messages, "deprecated", info = nm)
    }
})

test_that("spectral_matching() warns once, even when it builds its database", {
    skip_if_no_mspuritydata()
    q <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite",
                     package = "msPurity")
    r <- .with_deprecations(tryCatch(
        suppressMessages(spectral_matching(q, library_db_pth = NA)),
        error = function(e) NULL))
    expect_identical(r$n, 1L)
})

test_that("purityX() writes no SQLite and so does not warn without saveEIC", {
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    suppressPackageStartupMessages(library(xcms))
    r <- .with_deprecations(suppressMessages(purityX(
        .golden_xcms("msms_only_xset_OLD.rds"), plotP = FALSE,
        xgroups = c(1, 2))))
    expect_identical(r$n, 0L)
})
