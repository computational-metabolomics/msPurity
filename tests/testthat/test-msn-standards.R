context("MSn standard fixtures")

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
