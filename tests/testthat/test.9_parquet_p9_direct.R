context("Parquet: direct writer")

.bridged <- function(expr) {
    old <- options(msPurity.parquet.route = "bridge")
    on.exit(options(old))
    suppressMessages(expr)
}

test_that("the in-memory table set equals the SQLite-staged one", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    suppressPackageStartupMessages(library(xcms))
    src <- .mzs(".mzs_src")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    for (stage in c("2_frag4feature_pa.rds",
                    "9_averageAllFragSpectra_with_filter_pa.rds")) {
        pa <- .mzs(".pa_update")(.golden_pa(stage))
        direct <- src(pa, x, NA, NA)
        bridged <- .bridged(src(pa, x, NA, NA))
        ## The source row names the database by the time it was written.
        direct$source$name <- bridged$source$name
        tabs <- setdiff(names(bridged), "s_peaks")
        expect_identical(direct[tabs], bridged[tabs], info = stage)
        ## Raw scan peaks that carry no flags are counted, not copied.
        keep <- !bridged$s_peaks$type %in% "scan" |
            bridged$s_peaks$sid %in% direct$s_peaks$sid
        raw <- attr(direct$s_peaks, "raw_scan_peaks")
        expect_identical(c(raw %||% 0L), sum(!keep), info = stage)
        expect_equal(direct$s_peaks, bridged$s_peaks[keep, ],
                     check.attributes = FALSE, info = stage)
    }
})

test_that("createDatabase(format = 'parquet') keeps all it wrote via SQLite", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    suppressPackageStartupMessages(library(xcms))
    writers <- .parquet_golden_writers()
    for (nm in names(writers)) {
        ref <- readRDS(file.path(.parquet_golden_dir(), paste0(nm, ".rds")))
        got <- .parquet_golden_run(nm, writers)
        expect_true(isTRUE(.parquet_golden_contains(got, ref)),
                    info = paste(nm, paste(.parquet_golden_contains(got, ref),
                                           collapse = "; ")))
        if (nm == "xcmsnexp")
            ## The bridge, kept for comparison, writes the same.
            expect_identical(.bridged(.parquet_golden_run(nm, writers)), got)
    }
})
