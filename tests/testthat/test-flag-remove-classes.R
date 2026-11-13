test_that("XCMSnExp handles the reported metadata and single-sample classes", {
    skip_if_not_installed("xcms")
    skip_if_not_installed("MSnbase")
    skip_if_not_installed("msPurityData")
    suppressPackageStartupMessages(library(xcms))
    x <- readRDS(system.file("extdata", "tests", "purityA",
        "10_input_filterflagremove.rds", package = "msPurity", mustWork = TRUE))
    # MSnbase has no public file-name replacement method for this legacy
    # fixture. Relocate its files so XCMSnExp validity checks can find them.
    x@processingData@files <- file.path(system.file("extdata", "lcms", "mzML",
        package = "msPurityData", mustWork = TRUE),
        c("LCMS_1.mzML", "LCMS_2.mzML", "LCMSMS_1.mzML", "LCMSMS_2.mzML"))
    pd <- MSnbase::pData(x)
    pd$class <- factor(c("blank", "meoh", "sample", "sample"))
    pd$injectionOrder <- seq_len(nrow(pd))
    pd$mode <- "positive"
    MSnbase::pData(x) <- pd
    x <- xcms::groupChromPeaks(x, param = xcms::PeakDensityParam(
        sampleGroups = pd$class, minFraction = 0, bw = 5, binSize = 0.017))
    values <- xcms::featureValues(x, method = "medret", value = "into")

    for (reference in c("blank", "meoh")) {
        case <- paste("reference:", reference)
        result <- suppressMessages(flag_remove(x, ref.class = reference,
            remove_spectra_bool = FALSE))
        tab <- result$grp_peaklist
        expect_equal(nrow(tab), nrow(values), info = case)
        expect_equal(tab$blank_median_I, unname(values[, 1]), info = case)
        expect_equal(tab$meoh_median_I, unname(values[, 2]), info = case)
        expect_true(all(is.na(tab[, c("blank_RSD_I", "meoh_RSD_I",
            "blank_RSD_RT", "meoh_RSD_RT")])), info = case)
        expect_identical(MSnbase::pData(result$xcmsObj), pd, info = case)
    }
})
