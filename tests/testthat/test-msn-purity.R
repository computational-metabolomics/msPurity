.msn_purity_fixture <- function() {
    source <- .ancestry_spectra()
    source$precursorMz <- c(205.0972, 150.1234, 75.0567)
    source$isolationWindowLowerMz <- c(204.5, 149.9, 74.9)
    source$isolationWindowUpperMz <- c(205.7, 150.5, 75.4)
    data <- Spectra::spectraData(source)
    data$mz <- list(c(150.1234, 150.4, 160.1),
        c(75.0567, 75.2, 100.1), 50.1234)
    data$intensity <- list(c(80, 20, 900), c(30, 70, 500), 10)
    data$inPurity <- c(0.2, NA, NA)
    Spectra::Spectra(data, source = Spectra::MsBackendMemory())
}

.estimate_msn_fixture <- function(source, ppm = 10, tolerance = 0,
    offsets = NA_real_) {
    ancestry <- msPurity:::.pa_msn_ancestry(
        msPurity:::.pa_msn_spectra_metadata(source),
        "sourceSpectrumId")
    msPurity:::.pa_msn_purity(source, ancestry, source$pid,
        ppm, tolerance, offsets)
}

test_that("diagnostic purity uses stored peaks and each child's window", {
    source <- .msn_purity_fixture()
    queued <- Spectra::filterIntensity(source, intensity = c(50, Inf))
    out <- .estimate_msn_fixture(queued)
    records <- out$parameters$observations
    expect_equal(records$purity, c(0.8, 0.3))
    expect_identical(records$parent_pid, c(7L, 40L))
    expect_identical(records$parent_peak_index, c(1L, 1L))
    expect_identical(records$status, c("assessed", "assessed"))
    expect_equal(out$spectra$msPurityStagePurity, c(0.2, 0.8, 0.3))
    expect_identical(out$spectra$inPurity, source$inPurity)
    # Resetting a copy for assessment must not change the caller's queue.
    expect_identical(lengths(queued), c(2L, 2L, 0L))
    overridden <- .estimate_msn_fixture(source, offsets = c(0.05, 0.05))
    expect_equal(overridden$parameters$observations$purity, c(1, 1))
    expect_identical(overridden$parameters$observations$window_source,
        c("caller", "caller"))
    expect_error(.estimate_msn_fixture(source, offsets = c(NA, 0.5)),
        "Explicit isolation offsets")
    expect_error(.estimate_msn_fixture(source, offsets = c(NaN, NaN)),
        "Explicit isolation offsets")
})

test_that("diagnostics distinguish missing bounds, parents, peaks and ties", {
    cases <- list(
        window = function(x) {
            x$isolationWindowUpperMz[2L] <- NA
            x
        },
        parent = function(x) {
            x$precScanNum[2L] <- 999L
            x
        },
        missing = function(x) {
            x$precursorMz[2L] <- 150.15
            x
        },
        profile = function(x) {
            x$centroided <- c(FALSE, TRUE, TRUE)
            x
        },
        zero = function(x) {
            x$intensity <- list(c(0, 20, 900), c(30, 70, 500), 10)
            x
        },
        invalid = function(x) {
            x$intensity <- list(c(80, -20, 900), c(30, 70, 500), 10)
            x
        },
        empty = function(x) {
            data <- Spectra::spectraData(x)
            data$mz[1L] <- list(numeric())
            data$intensity[1L] <- list(numeric())
            Spectra::Spectra(data, source = Spectra::MsBackendMemory())
        },
        tie = function(x) {
            x$mz <- list(c(150.1234, 150.1234, 160.1),
                c(75.0567, 75.2, 100.1), 50.1234)
            x
        })
    expected <- c(window = "invalid_isolation_window",
        parent = "missing_parent",
        missing = "missing_precursor_peak", profile = "profile_parent",
        zero = "zero_precursor_intensity", invalid = "invalid_parent_peaks",
        empty = "empty_parent",
        tie = "ambiguous_precursor_peak")
    for (name in names(cases)) {
        source <- cases[[name]](.msn_purity_fixture())
        out <- .estimate_msn_fixture(source)$parameters$observations
        expect_identical(out$status[1L], unname(expected[name]), info = name)
        expect_true(is.na(out$purity[1L]), info = name)
    }
    shifted <- cases$missing(.msn_purity_fixture())
    out <- .estimate_msn_fixture(shifted, ppm = 0, tolerance = 0.03)
    expect_equal(out$parameters$observations$purity, c(0.8, 0.3))
})
