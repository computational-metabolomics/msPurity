.blink_spectra <- function(peaks) {
    d <- S4Vectors::DataFrame(msLevel = rep(2L, length(peaks)),
        precursorMz = 200, polarity = 1L, accession = paste0("entry-", seq_along(peaks)))
    d$mz <- lapply(peaks, function(p) as.numeric(unlist(p[[1L]])))
    d$intensity <- lapply(peaks, function(p) as.numeric(unlist(p[[2L]])))
    Spectra::Spectra(d, source = Spectra::MsBackendMemory())
}

test_that("CPU BLINK agrees with the pinned upstream cosine and count", {
    f <- jsonlite::read_json(test_path("fixtures", "blink", "reference.json"))
    q <- .blink_spectra(list(f$query))
    r <- .blink_spectra(f$references)
    for (config in f$cases) {
        p <- do.call(BlinkMatchingParam, config[c("binWidth", "blurTolerance", "intensityPower")])
        res <- spectralMatching(q, r, referenceRole = "library", param = p)
        hit <- res$matchedResults
        expect_identical(hit$referenceIndex, seq_along(r))
        expect_lte(max(abs(hit$blinkCosine - unlist(config$cosine))), 1e-12)
        expect_equal(hit$blinkMatchingIonCount, unlist(config$count), tolerance = 0)
        expect_equal(hit$rank, rank(-unlist(config$cosine), ties.method = "min"))
        expect_identical(res$queryStatus$nScored, 4L)
        expect_identical(res$scoring$referenceRevision, f$revision)
        expect_identical(res$scoring$settings, as.list(p))
        expect_identical(res$scoring$rankScore, "blinkCosine")
        expect_equal(scoreSpectrumPair(q, r, p)$scores,
            unlist(hit[1, c("blinkCosine", "blinkMatchingIonCount")]))
    }
})

test_that("BLINK respects the input queues, batch boundaries and selection", {
    old <- options(msPurity.parquet.batch_spectra = 2L)
    on.exit(options(old))
    same <- list(c(50, 80), c(100, 20))
    other <- list(90, 10)
    q <- .blink_spectra(list(same, same, same))
    q$precursorMz <- c(200, 205, 200)
    q$intensity <- list(c(100, 20), c(100, 20), c(0, 0))
    r <- .blink_spectra(list(other, same, other, same, same))
    r$polarity <- c(1L, 1L, 1L, 1L, 0L)
    r <- Spectra::filterMzRange(r, c(45, 60))
    res <- spectralMatching(q, r, topn = 1L, param = BlinkMatchingParam())
    expect_identical(res$matchedResults$referenceIndex, c(2L, 4L))
    expect_equal(res$matchedResults$rank, c(1, 1))
    expect_equal(res$matchedResults$blinkMatchingIonCount, c(1, 1))
    expect_identical(res$queryStatus$status, c("scored", "no_eligible_reference", "unscorable"))
    expect_equal(res$queryStatus$nEligible, c(4, 0, 4))
    empty <- spectralMatching(q, r[integer()], param = BlinkMatchingParam())
    expect_null(empty$matchedResults)
    expect_identical(empty$scoring$method, "BLINK")
    expect_identical(empty$queryStatus$nEligible, integer(3))
    expect_error(spectralMatching(q, r, param = BlinkMatchingParam(), q_ppmProd = 10), "blurTolerance")
    expect_error(spectralMatching(q, r, param = BlinkMatchingParam(), raW = 0.5), "intensityPower")
})

test_that("unit normalization preserves fractional powers across a large intensity range", {
    query <- .blink_spectra(list(list(c(50, 80), c(1e-300, 1e300))))
    ref <- .blink_spectra(list(list(50, 1)))
    p <- BlinkMatchingParam(intensityPower = 0.001)
    score <- scoreSpectrumPair(query, ref, p)$scores
    weights <- c(1e-300, 1e300) ^ p@intensityPower
    expect_equal(unname(score["blinkCosine"]), weights[1L] / sqrt(sum(weights^2)), tolerance = 1e-12)
    expect_equal(unname(score["blinkMatchingIonCount"]), 1)
})

test_that("BLINK validates its preprocessing and fragment units", {
    for (x in list(0, Inf, NA_real_, c(0.001, 0.002)))
        expect_error(BlinkMatchingParam(binWidth = x), "binWidth")
    expect_error(BlinkMatchingParam(blurTolerance = 0.0009), "blurTolerance")
    expect_error(BlinkMatchingParam(intensityPower = -1), "intensityPower")
    expect_error(BlinkMatchingParam(q_ppmProd = 1), "blurTolerance")
    expect_error(BlinkMatchingParam(l_daProd = 0.1), "blurTolerance")
})
