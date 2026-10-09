# Regenerate the golden SQLite dumps.
#
# Run from the package root against an installed msPurity whose SQLite
# writers are the reference, i.e. before any change to them:
#
#   R CMD INSTALL .
#   Rscript tests/testthat/fixtures/make-sqlite-golden.R
#
# The dumps were first generated from master at v1.37.3, before the Parquet
# backend was added. Regenerating them accepts the current output as correct,
# so only do it for an intended change to the SQLite schema or values.

suppressPackageStartupMessages({
    library(msPurity)
    library(xcms)
    library(testthat)
})
source(file.path("tests", "testthat", "helper-sqlite-golden.R"))

out <- file.path("tests", "testthat", "fixtures", "sqlite-golden")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (nm in names(.golden_writers())) {
    message("Dumping ", nm)
    d <- .golden_run(nm)
    saveRDS(d, file.path(out, paste0(nm, ".rds")), compress = "xz")
}
