.search_sp <- function(mz = 100, level = rep(2L, length(mz))) {
    n <- length(mz)
    d <- S4Vectors::DataFrame(msLevel = level, precursorMz = mz,
        rtime = as.numeric(seq_len(n)), acquisitionNum = seq_len(n),
        spectrumId = paste0("scan=", seq_len(n)), polarity = 0L, dataOrigin = "file")
    d$mz <- rep(list(c(25, 50)), n)
    d$intensity <- rep(list(c(2, 10)), n)
    Spectra::Spectra(d, source = Spectra::MsBackendMemory())
}

test_that("collection searches keep positions, selections, queues and tied ranks", {
    sp <- .search_sp(c(100, 100, 101, 100))
    lib <- sp[c(4, 3, 2)]
    lib$accession <- c("last", "other", "second")
    q <- Spectra::filterMzRange(sp[c(2, 1)], mz = c(20, 30))
    lib <- Spectra::filterMzRange(lib, mz = c(20, 30))
    hits <- spectralMatching(q, lib, topn = 1L)
    expect_identical(hits$query, q)
    expect_identical(hits$reference, lib)
    expect_identical(hits$referenceRole, "experimental")
    mr <- hits$matchedResults
    expect_identical(mr$queryIndex, c(1L, 1L, 2L, 2L))
    expect_identical(mr$referenceIndex, c(1L, 3L, 1L, 3L))
    expect_equal(mr$rank, rep(1, 4))
    expect_equal(mr$mcount, rep(1, 4))
    expect_identical(mr$query_native_id, q$spectrumId[mr$queryIndex])
    expect_identical(mr$reference_native_id, lib$spectrumId[mr$referenceIndex])
    expect_true(all(is.na(mr$library_accession)))
    mr <- spectralMatching(q, lib, referenceRole = "library", l_accessions = "second",
        q_rtrange = c(2, 2), l_rtrange = c(2, 2), rttol = 1)$matchedResults
    expect_identical(mr$queryIndex, 1L)
    expect_identical(mr$referenceIndex, 3L)
    expect_identical(mr$library_accession, "second")
    expect_null(spectralMatching(q, lib, l_pol = "positive")$matchedResults)
    expect_null(spectralMatching(q, lib, q_pol = "positive")$matchedResults)
    expect_null(spectralMatching(q[integer()], lib)$matchedResults)
    expect_null(spectralMatching(q, lib[integer()])$matchedResults)
    expect_null(spectralMatching(q, lib, topn = 0L)$matchedResults)
    expect_null(spectralMatching(q, Spectra::filterMzRange(lib, mz = c(0, 1)))$matchedResults)
})

test_that("collection index includes closed asymmetric windows at all supported tolerances", {
    q <- .search_sp(100)
    lib <- .search_sp(c(100.375 + 1e-10, 99.625, 100.375, 99.625 - 1e-10))
    p <- SpectralMatchingParam(q_ppmPrec = 1250, q_daPrec = 0.125, l_ppmPrec = 0,
        l_daPrec = 0.125)
    for (param in list(p, SpectralMatchingParam(l_ppmPrec = 1e6),
                       SpectralMatchingParam(l_ppmPrec = 2e6))) {
        want <- which(vapply(seq_along(lib), function(j)
            scoreSpectrumPair(q, lib, param, referenceIndex = j)$eligible, logical(1)))
        mr <- spectralMatching(q, lib, param = param)$matchedResults
        expect_identical(mr$referenceIndex, want)
    }
    expect_identical(spectralMatching(q, lib, param = p, l_ppmPrec = 0)$matchedResults$referenceIndex,
        2:3)
    expect_null(spectralMatching(q, lib, q_ppmPrec = 0, l_ppmPrec = 0)$matchedResults)
    expect_error(spectralMatching(q, lib, param = p, l_ppmPrec = 1), "conflicts")
})

test_that("ancestor conflicts and unsupported arguments are resolved before peak reads", {
    sp <- .search_sp(c(200, 150, 80, 210, 150, 80), rep(2:4, 2L))
    sp$dataOrigin <- rep(c("a", "b"), each = 3)
    sp$acquisitionNum <- rep(10:12, 2)
    sp$precScanNum <- rep(c(NA_integer_, 10L, 11L), 2)
    sp$spectrumId <- rep(c("root", "middle", "leaf"), 2)
    q <- Spectra::addProcessing(sp[c(3,1,2)], function(x, ...) stop("peaks read"))
    r <- Spectra::addProcessing(sp[c(4,6,5)], function(x, ...) stop("peaks read"))
    expect_null(spectralMatching(q, r)$matchedResults)
    expect_error(spectralMatching(q, r, usePrecursors = FALSE), "usePrecursors")
    expect_error(spectralMatching(q, r, updateDb = TRUE), "updateDb")
    expect_error(spectralMatching(q, r, copyDb = TRUE), "copyDb")
    expect_error(spectralMatching(q, r, q_xcmsGroups = 1), "q_xcmsGroups")
    expect_error(spectralMatching(q, r, outPth = tempfile()), "outPth")
    expect_error(spectralMatching(q, r, format = "sqlite"), "SQLite")
    expect_error(spectralMatching(q, r, l_accessions = "scan=1"), "library")
    same <- spectralMatching(sp[c(6,3,1,4,5,2)], sp)$matchedResults
    expect_identical(same$referenceIndex, c(6L, 3L, 1L, 4L, 5L, 2L))
    expect_identical(same$query_run_id, same$reference_run_id)
})
