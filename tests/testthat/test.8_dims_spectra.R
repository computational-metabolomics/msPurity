context("purityD Spectra slots, MsExperiment input and msfr files")

purityD_ref <- .parity_ref("purityD")

test_that("the DIMS workflow is unchanged with only the Spectra slots filled", {
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    q <- .parity_quiet
    inDF <- Getfiles(.dims_paths(), pattern = ".mzML", check = FALSE)
    pd <- q(averageSpectra(q(purityD(inDF, mzML = TRUE)), snMeth = "median",
                           snthr = 5))
    expect_identical(length(pd@avPeaks), 0L)
    expect_identical(averagedPeaks(pd), purityD_ref$averageSpectra)
    pd <- q(filterp(pd, thr = 5000, rsd = 10))
    expect_identical(averagedPeaks(pd), purityD_ref$filterp)
    pd <- q(subtract(pd))
    expect_identical(averagedPeaks(pd), purityD_ref$subtract)
    pd <- q(dimsPredictPurity(pd))
    expect_identical(averagedPeaks(pd), purityD_ref$dimsPredictPurity)
    expect_identical(getP(pd), purityD_ref$dimsPredictPurity)
    pd <- q(groupPeaks(pd))
    expect_identical(pd@groupedPeaks, purityD_ref$groupPeaks)

    sp <- averagedPeaks(pd, legacy = FALSE)
    expect_s4_class(sp, "Spectra")
    expect_length(sp, 4)
    expect_setequal(unique(sp$stage), c("orig", "processed"))
    expect_true("medianPurity" %in% Spectra::peaksVariables(sp))
})

test_that("purityD() accepts an MsExperiment", {
    q <- .parity_quiet
    inDF <- Getfiles(.dims_paths(), pattern = ".mzML", check = FALSE)
    me <- MsExperiment::readMsExperiment(
        inDF$filepth, sampleData = S4Vectors::DataFrame(inDF))
    pd <- q(purityD(me, mzML = TRUE))
    expect_identical(pd@fileList$name, inDF$name)
    expect_identical(pd@sampleIdx, q(purityD(inDF, mzML = TRUE))@sampleIdx)
    expect_s4_class(pd@experiment, "MsExperiment")
    pd <- q(averageSpectra(pd, snMeth = "median", snthr = 5))
    expect_identical(averagedPeaks(pd), purityD_ref$averageSpectra)

    # Without sample data the file paths, names and types are filled in.
    pd <- q(purityD(MsExperiment::readMsExperiment(inDF$filepth)))
    expect_identical(pd@fileList$filepth, normalizePath(inDF$filepth))
    expect_identical(pd@fileList$sampleType, c("sample", "sample"))
})

test_that("updateObject() adds the Spectra slots to a saved purityD object", {
    q <- .parity_quiet
    inDF <- Getfiles(.dims_paths(), pattern = ".mzML", check = FALSE)
    pd <- q(averageSpectra(q(purityD(inDF, mzML = TRUE)), snMeth = "median",
                           snthr = 5))
    saved <- pd
    attr(saved, "avSpectra") <- NULL
    attr(saved, "experiment") <- NULL
    expect_false(msPurity:::.pd_is_current(saved))
    up <- updateObject(saved)
    expect_true(msPurity:::.pd_is_current(up))
    expect_identical(averagedPeaks(up), saved@avPeaks)
    expect_identical(averagedPeaks(q(filterp(saved, thr = 5000, rsd = 10))),
                     purityD_ref$filterp)
})

test_that("msfr peak lists are read through Spectra", {
    files <- list.files(system.file("extdata", "dims", "msfr-peaks",
                                    package = "msPurityData"),
                        full.names = TRUE)
    for (f in files) {
        sp <- msPurity:::.msp_read_msfr(f)
        csv <- utils::read.csv(f)
        expect_s4_class(sp, "Spectra")
        expect_length(sp, length(unique(csv$scanid)))
        expect_true(all(c("snr", "noise") %in% Spectra::peaksVariables(sp)))
        expect_identical(msPurity:::.msp_msfr_frame(sp), csv)
    }
    expect_error(msPurity:::.msp_read_msfr(
        system.file("extdata", "tests", "external_annotations", "beams.tsv",
                    package = "msPurity")), "msfr peak list")
})
