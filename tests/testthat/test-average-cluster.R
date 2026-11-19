test_that("mean averaging accepts clusters with or without normalized intensities", {
    peaks <- data.frame(mz = 100, i = c(2, 6), snr = c(1, 2),
        inPurity = 1, ra = 100)
    average <- function(x) msPurity:::averageCluster(x, av = "mean",
        totalScans = 2L, normTIC = FALSE)

    expect_no_warning(raw <- average(peaks))
    expect_identical(raw[c("i", "inorm")], c(i = 4, inorm = NA_real_))

    peaks$inorm <- c(0.2, 0.6)
    expect_no_warning(normalized <- average(peaks))
    expect_equal(normalized[c("i", "inorm")], c(i = 4, inorm = 0.4))
})
