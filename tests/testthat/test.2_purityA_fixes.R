context("purityA bug fixes: averaging after rmp, createMSP max")

.frag_pa <- function() .golden_pa("2_frag4feature_pa.rds")

test_that("averaging after filterFragSpectra(rmp = TRUE) gives the rmp = FALSE results", {
    q <- .parity_quiet
    chain <- function(pa) q(averageAllFragSpectra(q(averageInterFragSpectra(
        q(averageIntraFragSpectra(pa))))))
    a <- chain(q(filterFragSpectra(.frag_pa(), plim = 0.7, snr = 3, rmp = TRUE)))
    b <- chain(q(filterFragSpectra(.frag_pa(), plim = 0.7, snr = 3)))
    av_a <- averagedSpectra(a)
    av_b <- averagedSpectra(b)
    expect_identical(names(av_a), names(av_b))
    # Averaging uses only passing peaks, so removing the others first must
    # not change the averaged spectra.
    for (g in names(av_b))
        for (lvl in c("av_inter", "av_all"))
            expect_identical(.parity_strip(av_a[[g]][[lvl]]),
                             .parity_strip(av_b[[g]][[lvl]]),
                             label = paste(g, lvl))
    # Some features lose every peak, and their averages are NULL.
    expect_true(any(vapply(av_a, function(x) is.null(x$av_all), logical(1))))
})

.msp_peaks <- function(f) {
    l <- readLines(f)
    p <- l[grepl("^[0-9.]+\t[0-9.]+", l)]
    as.numeric(sub("\t.*$", "", p))
}

test_that("createMSP(method = \"max\") writes the scan with the highest precursor intensity", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
                              "9_averageAllFragSpectra_with_filter_pa.rds",
                              package = "msPurity"))
    tab <- purityTable(pa)
    ms2 <- groupedSpectra(pa)
    # In grpid 147 and 212 the old code picked another scan.
    for (g in c("147", "212")) {
        grpd <- pa@grped_df[as.character(pa@grped_df$grpid) == g, ]
        best <- which.max(tab$precursorIntensity[match(grpd$pid, tab$pid)])
        f <- tempfile(fileext = ".msp")
        createMSP(pa, msp_file_pth = f, method = "max",
                  xcms_groupids = as.numeric(g), filter = FALSE)
        expect_equal(.msp_peaks(f), unname(ms2[[g]][[best]][, "mz"]),
                     tolerance = 1e-6, label = g)
    }
})

test_that("createMSP(method = \"max\") handles one or no passing peaks", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
                              "9_averageAllFragSpectra_with_filter_pa.rds",
                              package = "msPurity"))
    # grpid 452: the chosen scan has a single passing peak.
    f <- tempfile(fileext = ".msp")
    createMSP(pa, msp_file_pth = f, method = "max", xcms_groupids = 452)
    expect_length(.msp_peaks(f), 1)
    # Every group, including those with no passing peak, which are skipped.
    f <- tempfile(fileext = ".msp")
    expect_error(createMSP(pa, msp_file_pth = f, method = "max"), NA)
    expect_false(any(grepl("NUM_PEAK: 0$|Num Peaks: 0$", readLines(f))))
})
