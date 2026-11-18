# Regenerate the Parquet references.
#
# Run from the package root with the code whose output is the reference:
#
#   Rscript tests/testthat/fixtures/make-parquet-golden.R [name ...]
#
# The references were first generated from the SQLite-staged writer (master
# at 4706bac), before createDatabase(format = "parquet") wrote directly.
# Regenerating them accepts the current output as correct, so only do it for
# an intended change in results. The v1.0 dataset is the 1.0.0 results
# layout, kept to check that datasets written then stay readable and valid.

suppressPackageStartupMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(xcms)
    library(testthat)
})
for (h in c("helper-sqlite-golden.R", "helper-parquet.R",
            "helper-msn-standards.R", "helper-parquet-golden.R"))
    source(file.path("tests", "testthat", h))

out <- file.path("tests", "testthat", "fixtures", "parquet-golden")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writers <- .parquet_golden_writers()
names <- commandArgs(trailingOnly = TRUE)
if (!length(names)) names <- names(writers)
for (nm in setdiff(names, "v1.0")) {
    message(nm)
    saveRDS(.parquet_golden_run(nm, writers), file.path(out, paste0(nm, ".rds")),
            compress = "xz")
}

# Only with `v1.0` named explicitly: the dataset must stay as written then.
if ("v1.0" %in% commandArgs(trailingOnly = TRUE)) {
    message("v1.0")
    td <- tempfile()
    dir.create(td)
    w <- .d02_ms2_workflow()
    p <- createDatabase(w$pa, w$xe, outDir = td, dbName = "d02_ms2.parquet",
                        format = "parquet")
    dest <- file.path("tests", "testthat", "fixtures", "parquet-v1.0")
    unlink(dest, recursive = TRUE)
    dir.create(dest)
    file.copy(p, dest, recursive = TRUE)
}
