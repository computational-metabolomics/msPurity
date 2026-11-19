.px_input_fixture <- function() {
    skip_if_not_installed("xcms")
    meta <- S4Vectors::DataFrame(
        msLevel = rep(c(1L, 2L, 1L, 1L), 2),
        rtime = rep(c(1, 2, 3, 4), 2),
        rtime_adjusted = rep(c(11, 12, 13, 14), 2),
        dataOrigin = rep("same-source", 8),
        acquisitionNum = rep(c(101L, 102L, 103L, 104L), 2),
        spectrumId = rep(paste0("scan=", 101:104), 2))
    arrays <- lapply(seq_len(8), function(i) {
        cbind(mz = c(100.1234, 100.3234), intensity = c(i * 10, 5))
    })
    x <- xcms::XcmsExperiment()
    MsExperiment::sampleData(x) <- S4Vectors::DataFrame(
        sample = c("sample_a", "sample_b"), row.names = c("a", "b"))
    MsExperiment::spectra(x) <- .msp_memory(meta, arrays)
    x <- MsExperiment::linkSampleData(x, with = "spectra",
        sampleIndex = rep(1:2, each = 4), withIndex = 1:8)
    cp <- cbind(mz = c(100.1234, 100.1234),
        mzmin = c(100.123, 100.123), mzmax = c(100.124, 100.124),
        rt = c(13, 13), rtmin = c(11, 11), rtmax = c(14, 14),
        into = c(100, 300), sample = 1:2, sn = c(NA_real_, NA_real_))
    rownames(cp) <- c("peak-a", "peak-b")
    xcms::chromPeaks(x) <- cp
    xcms::chromPeakData(x) <- S4Vectors::DataFrame(
        ms_level = c(1L, 1L), is_filled = c(FALSE, FALSE),
        row.names = rownames(cp))
    fd <- data.frame(mzmed = 100.1234, rtmed = 13, row.names = "feature-a")
    fd$peakidx <- list(c(1L, 2L))
    xcms::featureDefinitions(x) <- fd
    x
}

test_that("purityX preserves modern sample links, times and queued peaks", {
    x <- .px_input_fixture()
    sp <- Spectra::filterIntensity(MsExperiment::spectra(x), intensity = 10)
    MsExperiment::spectra(x) <- sp
    prepared <- .px_prepare(x, singleFile = 2L)

    expect_identical(prepared$targets$sample, 2)
    expect_identical(prepared$targets$peak_id, "peak-b")
    expect_identical(prepared$targets$scan_positions[[1L]], c(1L, 3L, 4L))
    expect_identical(prepared$tracking$sourcePosition, 5:8)
    expect_equal(prepared$tracking$rt_corrected, c(11, 12, 13, 14))
    expect_equal(prepared$tracking$retentionTime, c(1, 2, 3, 4))
    expect_equal(prepared$scanpeaks[[2L]][[1L]],
        cbind(mz = 100.1234, intensity = 50))
    expect_length(prepared$exclusions$reason, 0L)
})

test_that("purityX separates feature membership from excluded inputs", {
    x <- .px_input_fixture()
    fd <- xcms::featureDefinitions(x)
    fd$peakidx <- list(c(1L, 1L, 2L))
    xcms::featureDefinitions(x) <- fd
    p <- .px_prepare(x, fileignore = 2L)
    expect_identical(p$targets$c_peak_id, c(1L, 1L))
    expect_equal(p$groups$i, 200)
    expect_identical(p$exclusions$reason, "ignored_sample")

    cd <- xcms::chromPeakData(x)
    cd$is_filled[1L] <- TRUE
    cd$ms_level[2L] <- 2L
    xcms::chromPeakData(x) <- cd
    empty <- .px_prepare(x)
    expect_equal(nrow(empty$targets), 0L)
    expect_identical(empty$exclusions$reason,
        c("filled_peak", "filled_peak", "non_ms1_peak"))
    expect_equal(nrow(empty$groups), 1L)
    expect_length(empty$tracking$scan, 0L)
})

test_that("purityX uses explicit shared sample links and empty selections", {
    x <- .px_input_fixture()
    expect_warning(x <- MsExperiment::linkSampleData(x, with = "spectra",
            sampleIndex = c(rep(1L, 4), rep(2L, 4)), withIndex = rep(1:4, 2)),
        "Overwriting previously present link")
    p <- .px_prepare(x)
    expect_identical(p$tracking$sourcePosition, rep(1:4, 2))
    expect_equal(p$scanpeaks[[1L]], p$scanpeaks[[2L]])
    expect_equal(nrow(.px_prepare(x, xgroups = integer())$groups), 0L)
    expect_error(.px_prepare(x, singleFile = 3), "valid integer positions")
    expect_error(.px_prepare(x, xgroups = 0), "valid integer positions")
})

test_that("shared origins preserve xcms medret feature intensities", {
    x <- .px_input_fixture()
    cp <- xcms::chromPeaks(x)[c(1L, 1L, 2L), , drop = FALSE]
    rownames(cp) <- c("early", "late", "second-sample")
    cp[, "rt"] <- c(11, 14, 13)
    cp[, "into"] <- c(100, 900, 300)
    xcms::chromPeaks(x) <- cp
    fd <- data.frame(mzmed = 100.1234, rtmed = c(11, 12.5),
        row.names = c("stored-centre", "tied-distance"))
    fd$peakidx <- list(1:3, c(2L, 1L, 3L))
    xcms::featureDefinitions(x) <- fd

    # For the tie at 12.5 seconds, the late peak is the first member.
    shared <- unname(.px_prepare(x)$groups$i)
    expect_equal(shared, c(200, 600))
    sp <- MsExperiment::spectra(x)
    Spectra::dataOrigin(sp) <- rep(c("one", "two"), each = 4L)
    MsExperiment::spectra(x) <- sp
    intensities <- xcms::featureValues(x, method = "medret", value = "into")
    expect_equal(shared, unname(apply(intensities, 1L, median)))
})

test_that("purityX retains the second original LC-MS sample identity", {
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    x <- readRDS(system.file("extdata", "tests", "xcms",
            "ms_only_xset_OLD.rds", package = "msPurity"))
    paths <- list.files(system.file("extdata", "lcms", "mzML",
            package = "msPurityData"), pattern = "^LCMS_[12][.]mzML$",
        full.names = TRUE)
    x@filepaths <- paths
    p <- .px_prepare(x, singleFile = 2L)
    expect_true(all(p$targets$sample == 2L))
    expect_true(all(p$tracking$sample == 2L))
    expect_true(all(lengths(p$targets$scan_positions) > 0L))
    expect_identical(p$inputInfo$source_files, paths)
    expect_equal(p$tracking$retentionTime, x@rt$raw[[2L]])

    xn <- readRDS(system.file("extdata", "tests", "xcms",
            "ms_only_xcmsnexp.rds", package = "msPurity"))
    # Relocate the existing fixture without changing its retained spectra.
    xn@processingData@files <- paths
    pn <- .px_prepare(xn, singleFile = 2L)
    retained <- as.integer(which(MSnbase::fromFile(xn) == 2L))
    expect_identical(pn$tracking$sourcePosition, retained)
    expect_equal(pn$tracking$retentionTime,
        as.numeric(xcms::rtime(xn, adjusted = FALSE)[retained]))
    expect_identical(pn$tracking$acquisitionNum,
        as.integer(MSnbase::acquisitionNum(xn)[retained]))
})
