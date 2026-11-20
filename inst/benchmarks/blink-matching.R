# See blink-matching.md for timing and peak RSS commands.
suppressPackageStartupMessages(library(msPurity))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
mode <- match.arg(args[1], c("prepare", "default", "blink"))
dir <- args[2]
lib <- file.path(dir, "library")
queries <- file.path(dir, "queries.rds")
options(msPurity.parquet.batch_spectra = 8L)
if (mode == "prepare") {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    # Four precursor groups and overlapping fragment series with unequal weights.
    i <- seq_len(128L)
    d <- S4Vectors::DataFrame(msLevel = rep(2L, length(i)), polarity = 1L,
        precursorMz = 200 + (i %% 4L), accession = sprintf("entry-%03d", i))
    d$mz <- lapply(i, function(k) seq_len(64L) + 50 + (k %% 5L) * 0.002)
    d$intensity <- lapply(i, function(k) as.numeric(1 + (seq_len(64L) * k) %% 97L))
    spectra <- Spectra::Spectra(d, source = Spectra::MsBackendMemory())
    saveRDS(spectra[1:8], queries)
    convertLibraryToParquet(spectra, lib, version = "blink-benchmark-1")
} else {
    q <- readRDS(queries)
    p <- if (mode == "blink") BlinkMatchingParam() else SpectralMatchingParam()
    elapsed <- system.time(result <- spectralMatching(q, lib, param = p))
    cat(sprintf("method=%s queries=%d references=%d eligible=%d scored=%d elapsed_s=%.3f\n",
        mode, length(q), length(result$reference), sum(result$queryStatus$nEligible),
        sum(result$queryStatus$nScored), elapsed[["elapsed"]]))
}
