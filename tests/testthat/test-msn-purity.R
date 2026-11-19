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

.estimate_msn_fixture <- function(source, pids = source$pid, ...) {
    ancestry <- msPurity:::.pa_msn_ancestry(
        msPurity:::.pa_msn_spectra_metadata(source),
        "sourceSpectrumId")
    msPurity:::.pa_msn_purity(source, ancestry, pids, ...)
}

test_that("parent purity uses each child's window and records its evidence", {
    source <- .msn_purity_fixture()
    out <- .estimate_msn_fixture(source)
    records <- out$parameters$observations
    expect_equal(records$purity, c(0.8, 0.3))
    expect_identical(records$parentPid, c(7L, 40L))
    expect_identical(records$parent_peak_index, c(1L, 1L))
    expect_identical(records$status, c("assessed", "assessed"))
    expect_identical(records$childLineageScanId, c("102", "103"))
    expect_identical(records$parentLineageScanId, c("101", "102"))
    expect_equal(out$spectra$msPurityStagePurity, c(0.2, 0.8, 0.3))
    expect_equal(out$spectra$inPurity, c(0.2, 0.8, 0.3))
    overridden <- .estimate_msn_fixture(source, offsets = c(0.05, 0.05))
    expect_equal(overridden$parameters$observations$purity, c(1, 1))
    expect_identical(overridden$parameters$observations$window_source,
        c("caller", "caller"))
    expect_error(.estimate_msn_fixture(source, offsets = c(NA, 0.5)),
        "Explicit isolation offsets")
    expect_error(.estimate_msn_fixture(source, offsets = c(NaN, NaN)),
        "Explicit isolation offsets")
    selected <- .estimate_msn_fixture(source, pids = 41L)
    expect_equal(selected$spectra$inPurity, c(0.2, NA_real_, 0.3))
    expect_identical(selected$spectra$msPurityPurityStatus,
        c("assessed", "not_assessed", "assessed"))
})

test_that("parent purity evaluates supplied processing", {
    source <- .msn_purity_fixture()
    queued <- Spectra::filterIntensity(source, intensity = c(50, Inf))
    out <- .estimate_msn_fixture(queued)

    # Only one in-window peak survives in each measured parent.
    expect_equal(out$parameters$observations$purity, c(1, 1))
    expect_equal(out$parameters$observations$matched_mz, c(150.1234, 75.2))
    expect_identical(lengths(queued), c(2L, 2L, 0L))
    expect_identical(queued$inPurity, source$inPurity)
})

test_that("selection reuses intensity preference without a mass cutoff", {
    source <- .msn_purity_fixture()
    source$precursorMz[2L] <- 150.15
    shifted <- .estimate_msn_fixture(source)$parameters$observations
    expect_equal(shifted$purity, c(0.8, 0.3))
    expect_equal(shifted$matched_mz[1L], 150.1234)

    source$mz <- list(c(150.1235, 150.2234, 150.4),
        c(75.0567, 75.2, 100.1), 50.1234)
    source$intensity <- list(c(30, 80, 90), c(30, 70, 500), 10)
    closest <- .estimate_msn_fixture(source)$parameters$observations
    source$precursorIntensity <- c(NA_real_, 80, NA_real_)
    recorded <- .estimate_msn_fixture(source)$parameters$observations
    strongest <- .estimate_msn_fixture(source,
        mostIntense = TRUE)$parameters$observations

    expect_equal(closest$purity[1L], 30 / 200)
    expect_equal(recorded$purity[1L], 80 / 200)
    expect_equal(strongest$purity[1L], 90 / 200)
    selected <- c(closest$parent_peak_index[1L], recorded$parent_peak_index[1L],
        strongest$parent_peak_index[1L])
    expect_identical(selected, 1:3)
    # Equal masses cannot identify a unique numerator in pcalc.
    source$mz <- list(c(150.2234, 150.2234, 150.4),
        c(75.0567, 75.2, 100.1), 50.1234)
    duplicated <- .estimate_msn_fixture(source)$parameters$observations
    expect_identical(duplicated$status[1L], "ambiguous_precursor_peak")
    expect_true(is.na(duplicated$purity[1L]))
})

test_that("outside precursors retain native bounds without an override", {
    source <- .msn_purity_fixture()
    source$precursorMz[2L] <- 150.6
    source$isolationWindowTargetMz <- c(NA_real_, 150.2, NA_real_)
    for (most_intense in c(FALSE, TRUE)) {
        out <- .estimate_msn_fixture(source,
            mostIntense = most_intense)$parameters$observations
        expect_identical(out$status[1L], "precursor_outside_window",
            info = paste("mostIntense", most_intense))
        expect_true(is.na(out$purity[1L]),
            info = paste("mostIntense", most_intense))
    }
    overridden <- .estimate_msn_fixture(source,
        offsets = c(0.5, 0.5))$parameters$observations
    expect_equal(overridden$purity[1L], 20 / 100)
    expect_equal(overridden$precursor_mz[1L], 150.6)
    expect_equal(overridden$native_isolation_target_mz[1L], 150.2)
    expect_equal(overridden$native_lower_mz[1L], 149.9)
    expect_equal(overridden$native_upper_mz[1L], 150.5)
    expect_equal(overridden$lower_mz[1L], 150.1)
    expect_equal(overridden$upper_mz[1L], 151.1)
    expect_identical(overridden$window_source[1L], "caller")
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
            x$mz <- list(c(149.8, 150.6, 160.1),
                c(75.0567, 75.2, 100.1), 50.1234)
            x
        },
        boundary = function(x) {
            x$mz <- list(c(149.9, 150.6, 160.1),
                c(75.0567, 75.2, 100.1), 50.1234)
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
        missing = "missing_precursor_peak", boundary = "invalid_estimate",
        profile = "profile_parent",
        zero = "zero_precursor_intensity", invalid = "invalid_parent_peaks",
        empty = "empty_parent",
        tie = "ambiguous_precursor_peak")
    for (name in names(cases)) {
        source <- cases[[name]](.msn_purity_fixture())
        out <- .estimate_msn_fixture(source)$parameters$observations
        expect_identical(out$status[1L], unname(expected[name]), info = name)
        expect_true(is.na(out$purity[1L]), info = name)
    }
})
