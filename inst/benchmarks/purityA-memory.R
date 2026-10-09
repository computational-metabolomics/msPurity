# purityA workflow: time and memory with and without the legacy slots.
#
# Runs purityA(), frag4feature(), filterFragSpectra(allfrag = TRUE) and the
# three averaging steps on the two LC-MS/MS files of msPurityData, once with
# options(msPurity.legacySlots = TRUE) (the default) and once with FALSE,
# and reports the elapsed time of each step and the in-memory size of the
# final object and of each of its slots.
#
# Usage:
#
#   Rscript inst/benchmarks/purityA-memory.R [outDir]
#
# outDir defaults to a temporary directory; it receives results.csv.

suppressPackageStartupMessages({
    library(msPurity)
    library(xcms)
})

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[1] else tempfile("pa-bench-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

msms <- list.files(system.file("extdata", "lcms", "mzML",
                               package = "msPurityData"),
                   full.names = TRUE, pattern = "MSMS")
xcmsObj <- readRDS(system.file("extdata", "tests", "xcms",
                               "msms_only_xcmsnexp.rds", package = "msPurity"))
xcmsObj@processingData@files <- msms

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

run <- function(legacy) {
    old <- options(msPurity.legacySlots = legacy)
    on.exit(options(old))
    times <- list()
    step <- function(name, expr) {
        t <- system.time(v <- quiet(expr))[["elapsed"]]
        times[[name]] <<- t
        v
    }
    pa <- step("purityA", purityA(msms))
    pa <- step("frag4feature", frag4feature(pa, xcmsObj))
    pa <- step("filterFragSpectra", filterFragSpectra(pa, plim = 0.7, snr = 3,
                                                      allfrag = TRUE))
    pa <- step("averageIntraFragSpectra", averageIntraFragSpectra(pa))
    pa <- step("averageInterFragSpectra", averageInterFragSpectra(pa))
    pa <- step("averageAllFragSpectra", averageAllFragSpectra(pa))
    slots <- c("puritydf", "grped_df", "grped_ms2", "all_frag_scans",
               "av_spectra", "spectra", "fragSpectra", "avSpectra")
    size <- vapply(slots, function(s)
        as.numeric(utils::object.size(methods::slot(pa, s))), numeric(1))
    rbind(
        data.frame(legacySlots = legacy, measure = "seconds",
                   item = names(times), value = unlist(times)),
        data.frame(legacySlots = legacy, measure = "bytes",
                   item = c(slots, "object"),
                   value = c(size, as.numeric(utils::object.size(pa)))))
}

res <- rbind(run(TRUE), run(FALSE))
rownames(res) <- NULL
utils::write.csv(res, file.path(out, "results.csv"), row.names = FALSE)
print(res)
message("Results written to ", normalizePath(out))
