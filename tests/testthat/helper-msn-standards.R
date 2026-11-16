.msn_standard_file <- function(compound, kind) {
    identifiers <- c(tryptophan = "6305", kynurenine = "161166")
    system.file("extdata", "msn-standards",
        paste0("msPurity_L_", compound, "_", identifiers[[compound]], "_", kind,
            ".mzML.gz"), package = "msPurity", mustWork = TRUE)
}

.msn_standard_mz <- c(tryptophan = 205.097154101381,
    kynurenine = 209.092068720951)

.msn_standard_counts <- list(
    tryptophan = list(MS1 = c(193L, 0L, 0L, 0L),
        MSn = c(34L, 78L, 44L, 20L)),
    kynurenine = list(MS1 = c(191L, 0L, 0L, 0L),
        MSn = c(35L, 51L, 33L, 16L))
)

.msn_standard_cases <- expand.grid(kind = c("MS1", "MSn"),
    compound = names(.msn_standard_mz), stringsAsFactors = FALSE)

.msn_standard_feature <- local({
    cache <- list()
    function(compound, kind) {
        key <- paste(compound, kind, sep = "_")
        if (!is.null(cache[[key]])) return(cache[[key]])

        xe <- MsExperiment::readMsExperiment(.msn_standard_file(compound, kind),
            BPPARAM = BiocParallel::SerialParam())
        xe <- xcms::findChromPeaks(xe,
            param = xcms::CentWaveParam(ppm = 10, mzdiff = 0.001,
                peakwidth = c(3, 30), snthresh = 10,
                prefilter = c(3, 100), fitgauss = FALSE, noise = 0,
                verboseColumns = FALSE, firstBaselineCheck = TRUE),
            msLevel = 1L, BPPARAM = BiocParallel::SerialParam())
        target <- .msn_standard_mz[[compound]]
        cp <- xcms::chromPeaks(xe)
        selected <- which(abs(cp[, "mz"] - target) <= target * 10 * 1e-6)
        feature <- integer()
        if (length(selected) == 1L) {
            xe <- xcms::groupChromPeaks(xe,
                param = xcms::PeakDensityParam(sampleGroups = 1L,
                    binSize = 0.005, bw = 0.5, minFraction = 0.5,
                    maxFeatures = 50))
            features <- xcms::featureDefinitions(xe)
            feature <- which(vapply(features$peakidx,
                    function(i) selected %in% i, logical(1L)))
        }
        cache[[key]] <<- list(xe = xe, selected = selected, feature = feature)
        cache[[key]]
    }
})
