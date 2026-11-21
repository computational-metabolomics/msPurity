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
