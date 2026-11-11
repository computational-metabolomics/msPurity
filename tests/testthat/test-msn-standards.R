context("MSn standard fixtures")

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

.msn_standard_feature <- function(compound, kind) {
    xe <- MsExperiment::readMsExperiment(.msn_standard_file(compound, kind))
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
    list(xe = xe, selected = selected, feature = feature)
}

for (i in seq_len(nrow(.msn_standard_cases))) local({
    compound <- .msn_standard_cases$compound[i]
    kind <- .msn_standard_cases$kind[i]

    test_that(paste(compound, kind, "retains expected spectra"), {
        sp <- Spectra::Spectra(.msn_standard_file(compound, kind),
            source = Spectra::MsBackendMzR(),
            BPPARAM = BiocParallel::SerialParam())
        level <- Spectra::msLevel(sp)
        expect_identical(tabulate(level, nbins = 4L),
            .msn_standard_counts[[compound]][[kind]])
        expect_false(anyDuplicated(sp$spectrumId) > 0L)
        expect_true(all(is.finite(Spectra::rtime(sp))))
        peaks <- Spectra::peaksData(sp)
        expect_true(all(vapply(peaks, function(p) {
            all(is.finite(p)) && all(p[, "intensity"] >= 0)
        }, logical(1L))))
        if (kind == "MSn") {
            parent <- match(sp$precScanNum[level > 1L], sp$acquisitionNum)
            expect_false(anyNA(parent))
            expect_identical(level[parent], level[level > 1L] - 1L)
            expect_true(any(lengths(sp[level == 4L]) == 0L))
            expect_true(any(lengths(sp[level == 4L]) > 0L))
        }
    })

    test_that(paste(compound, kind, "target peak is detected and grouped"), {
        skip_if_not_installed("xcms")
        prepared <- .msn_standard_feature(compound, kind)
        target <- .msn_standard_mz[[compound]]
        expect_length(prepared$selected, 1L)
        if (length(prepared$selected) == 1L) {
            cp <- xcms::chromPeaks(prepared$xe)
            expect_gte(cp[prepared$selected, "sn"], 10)
            expect_length(prepared$feature, 1L)
            if (length(prepared$feature) == 1L) {
                features <- xcms::featureDefinitions(prepared$xe)
                expect_lte(abs(features$mzmed[prepared$feature] - target),
                    target * 10 * 1e-6)
            }
        }
    })

    if (kind == "MSn") {
        test_that(paste(compound, "MS2 workflow on MSn data supports purity,",
                "association, filtering and averaging"), {
                skip_if_not_installed("xcms")
                prepared <- .msn_standard_feature(compound, kind)
                expect_length(prepared$selected, 1L)
                expect_length(prepared$feature, 1L)
                if (length(prepared$selected) == 1L &&
                        length(prepared$feature) == 1L) {
                    pa <- purityA(.msn_standard_file(compound, kind),
                        cores = 1L, plotP = FALSE)
                    purity <- purityTable(pa)
                    expect_identical(nrow(purity),
                        .msn_standard_counts[[compound]]$MSn[[2L]])
                    expect_true(all(is.finite(purity$inPurity)))
                    pa <- frag4feature(pa, prepared$xe, ppm = 10,
                        convert2RawRT = FALSE)
                    linked <- any(pa@grped_df$grpid == prepared$feature)
                    expect_true(linked)
                    if (linked) {
                        pa <- filterFragSpectra(pa)
                        pa <- averageAllFragSpectra(pa, cores = 1L)
                        averages <- averagedSpectra(pa, legacy = FALSE)
                        target_average <- averages[
                            averages$grpid %in% as.character(prepared$feature)]
                        expect_length(target_average, 1L)
                        expect_true(all(lengths(target_average) > 0L))
                    }
                }
            })
    }
})
