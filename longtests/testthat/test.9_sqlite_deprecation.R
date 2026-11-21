context("SQLite output regressions and deprecation warnings")

# Historical dumps retain their original values except verified raw pid
# corrections. Other tables, peak arrays, flags and row ordering stay unchanged.

for (nm in names(.golden_writers())) local({
    nm <- nm
    test_that(paste("SQLite output matches corrected references:", nm), {
        skip_if_no_mspuritydata()
        skip_if_not_installed("xcms")
        suppressPackageStartupMessages(library(xcms))
        golden <- file.path(.golden_dir(), paste0(nm, ".rds"))
        expect_true(file.exists(golden), info = nm)
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
