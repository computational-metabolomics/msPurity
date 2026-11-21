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
