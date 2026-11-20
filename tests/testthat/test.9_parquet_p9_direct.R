context("Parquet: direct writer")

test_that("createDatabase(format = 'parquet') keeps all it wrote via SQLite", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    suppressPackageStartupMessages(library(xcms))
    writers <- .parquet_golden_writers()
    for (nm in names(writers)) {
        ref <- readRDS(file.path(.parquet_golden_dir(), paste0(nm, ".rds")))
        got <- .parquet_golden_run(nm, writers)
        differences <- .parquet_golden_contains(got, ref, computed = nm == "d02_ms2")
        expect_true(isTRUE(differences),
                    info = paste(nm, paste(differences, collapse = "; ")))
    }
})
