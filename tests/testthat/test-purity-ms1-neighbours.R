.neighbour_spectra <- function(levels, parent = NA_integer_,
    empty = integer()) {
    target <- 205.0972
    n <- length(levels)
    d <- S4Vectors::DataFrame(msLevel = levels, scanIndex = seq_len(n),
        acquisitionNum = 100L + seq_len(n), rtime = as.numeric(seq_len(n)),
        precScanNum = rep(parent, length.out = n), precursorMz = target,
        precursorIntensity = 100, isolationWindowTargetMz = target,
        isolationWindowLowerMz = target - 0.5,
        isolationWindowUpperMz = target + 0.5,
        dataOrigin = "neighbour-fixture", centroided = TRUE)
    d$mz <- rep(list(c(target, target + 0.2)), n)
    d$intensity <- rep(list(c(100, 25)), n)
    d$mz[empty] <- list(numeric())
    d$intensity[empty] <- list(numeric())
    Spectra::Spectra(d, source = Spectra::MsBackendMemory())
}

test_that("one usable MS1 supports the nearest single-measurement estimate", {
    fixtures <- list(before = .neighbour_spectra(c(1L, 2L), c(NA, 101L)),
        after = .neighbour_spectra(c(2L, 1L)),
        empty_before = .neighbour_spectra(c(1L, 2L, 1L),
            c(NA, 101L, NA), empty = 1L))
    expected <- list(before = 101L, after = 102L, empty_before = 103L)
    for (name in names(fixtures)) {
        pa <- purityA(fixtures[[name]], isotopes = FALSE)
        expect_equal(purityTable(pa)$inPurity, 0.8, info = name)
        record <- processingParameters(pa, "purityA")[[1L]]$resolved[[1L]]
        expect_identical(record$observations$method, "single_ms1", info = name)
        expect_identical(record$observations$ms1_acquisitionNum,
            list(expected[[name]]), info = name)
    }
})

test_that("preceding-reference mode never substitutes a later MS1", {
    after <- .neighbour_spectra(c(2L, 1L))
    expect_warning(unavailable <- purityA(after, nearest = FALSE,
            isotopes = FALSE), "no usable permitted MS1")
    expect_true(is.na(purityTable(unavailable)$inPurity))
    record <- processingParameters(unavailable, "purityA")[[1L]]$resolved[[1L]]
    expect_identical(record$observations$ms1_acquisitionNum, list(integer()))

    before <- .neighbour_spectra(c(1L, 2L), c(NA, 101L))
    available <- purityA(before, nearest = FALSE, isotopes = FALSE)
    expect_equal(purityTable(available)$inPurity, 0.8)
    before$precScanNum[2L] <- 999L
    expect_warning(unresolved <- purityA(before, nearest = FALSE,
            isotopes = FALSE), "no usable permitted MS1")
    expect_true(is.na(purityTable(unresolved)$inPurity))
})

test_that("absent and empty MS1 measurements retain unavailable MS2 spectra", {
    absent <- .neighbour_spectra(c(2L, 3L, 4L))
    expect_message(pa <- purityA(absent, isotopes = FALSE),
        "No MS1 measurements")
    expect_identical(Spectra::msLevel(pa@spectra), c(2L, 3L, 4L))
    expect_true(all(is.na(purityTable(pa)$inPurity)))
    empty <- .neighbour_spectra(c(1L, 2L, 1L), c(NA, 101L, NA), c(1L, 3L))
    expect_warning(pa <- purityA(empty, isotopes = FALSE),
        "no usable permitted MS1")
    expect_true(all(is.na(purityTable(pa)$inPurity)))
    record <- processingParameters(pa, "purityA")[[1L]]$resolved[[1L]]
    expect_identical(record$observations$method, "no_usable_ms1")
    expect_identical(record$observations$ms1_scanIndex, list(integer()))
})

test_that("two usable neighbours retain the interpolation and measurement identities", {
    sp <- .neighbour_spectra(c(1L, 2L, 1L), c(NA, 101L, NA))
    sp$intensity <- list(c(100, 25), c(100, 25), c(100, 100))
    pa <- purityA(sp, isotopes = FALSE)
    expect_equal(purityTable(pa)$inPurity, (0.8 + 0.5) / 2)
    record <- processingParameters(pa, "purityA")[[1L]]$resolved[[1L]]
    expect_identical(record$observations$method, "linear")
    expect_identical(record$observations$ms1_acquisitionNum,
        list(c(101L, 103L)))
    expect_false(any(startsWith(names(purityTable(pa)), ".purity")))
})
