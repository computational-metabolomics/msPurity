.window_spectra <- function() {
    target <- 100.1234
    d <- S4Vectors::DataFrame(msLevel = c(1L, 2L, 2L, 1L),
        scanIndex = 1:4, acquisitionNum = 1:4, rtime = 1:4,
        precScanNum = c(NA_integer_, 1L, 1L, NA_integer_),
        precursorMz = c(NA_real_, target, target, NA_real_),
        precursorIntensity = c(NA_real_, 100, 100, NA_real_),
        isolationWindowTargetMz = c(NA_real_, target, target + 0.1, NA_real_),
        isolationWindowLowerMz = c(NA_real_, target - 0.4,
            target - 0.1, NA_real_),
        isolationWindowUpperMz = c(NA_real_, target + 0.4,
            target + 0.4, NA_real_),
        dataOrigin = "varying-windows", centroided = TRUE)
    d$mz <- list(target + c(-0.3, 0, 0.3), c(57.0344, 70.0456),
        c(57.0344, 70.0456), target + c(-0.3, 0, 0.3))
    d$intensity <- list(c(50, 100, 50), c(10, 20), c(10, 20), c(50, 100, 50))
    Spectra::Spectra(d, source = Spectra::MsBackendMemory())
}

test_that("each MS2 uses its own bounds for initial and linear purity", {
    spectra <- .window_spectra()
    initial <- assessPuritySingle(spectra, interpol = "none", isotopes = FALSE)
    interpolated <- purityA(spectra, isotopes = FALSE)
    settings <-
        processingParameters(interpolated, "purityA")[[1L]]$resolved[[1L]]

    expect_equal(initial$aPurity, c(1 / 2, 2 / 3))
    expect_equal(purityTable(interpolated)$inPurity, c(1 / 2, 2 / 3))
    expect_identical(settings$windows$seqNum, 2:3)
    expect_equal(settings$windows$isolationWindowTargetMz,
        c(100.1234, 100.2234))
    expect_equal(settings$windows$lowerOffset, c(0.4, 0.2))
    expect_true(all(is.na(settings$offsets)))
    overridden <- assessPuritySingle(spectra, offsets = c(0.4, 0.4),
        interpol = "none", isotopes = FALSE)
    expect_equal(overridden$aPurity, c(1 / 2, 1 / 2))
    spectra$isolationWindowTargetMz[2L] <- NA_real_
    resolved <- msPurity:::.msp_isolation_settings(spectra[2L])
    expect_true(is.na(resolved$isolationWindowTargetMz))
    expect_equal(resolved$lowerOffset, 0.4)
})

test_that("missing bounds stay unknown while valid events continue", {
    spectra <- .window_spectra()
    spectra$isolationWindowUpperMz[3L] <- NA_real_
    expect_warning(pa <- purityA(spectra, isotopes = FALSE),
        "1 MS2 spectra have missing or invalid isolation bounds")
    expect_equal(purityTable(pa)$inPurity, c(1 / 2, NA_real_))
    windows <- processingParameters(pa, "purityA")[[1L]]$resolved[[1L]]$windows
    expect_identical(windows$status, c("valid", "invalid_isolation_window"))
    # Vendor names cannot supply an acquisition width that was never recorded.
    spectra$instrumentModel <- "Agilent Q-TOF"
    expect_warning(result <- assessPuritySingle(spectra, isotopes = FALSE),
        "missing or invalid isolation bounds")
    expect_equal(result$inPurity, c(1 / 2, NA_real_))
    for (offsets in list(c(NA, 0.4), c(-0.1, 0.4), c(0, 0), 0.4))
        expect_error(assessPuritySingle(spectra, offsets = offsets),
            "Explicit isolation offsets", info = paste(offsets, collapse = ","))
})

test_that("normalization and history use the same event windows", {
    spectra <- .window_spectra()
    # A simple caller function makes the independent expected weights explicit.
    weight <- function(x) 1 - abs(x)
    environment(weight) <- baseenv()
    pa <- purityA(spectra, iwNorm = TRUE, iwNormFun = weight, isotopes = FALSE)
    # Window midpoints are 100.1234 and 100.2734, respectively.
    expect_equal(purityTable(pa)$inPurity,
        c(100 / (35 + 100 + 35), 85 / (85 + 42.5)))
    settings <- processingParameters(pa, "purityA")[[1L]]$resolved[[1L]]
    expect_identical(settings$iwNormFun, weight)
    generated <- purityA(spectra, iwNorm = TRUE, isotopes = FALSE)
    settings <- processingParameters(generated, "purityA")[[1L]]$resolved[[1L]]
    expect_length(settings$iwNormFunctions, 2L)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(settings, path)
    expect_identical(lapply(readRDS(path)$iwNormFunctions, function(f) f(0)),
        lapply(settings$iwNormFunctions, function(f) f(0)))
})

test_that("generated weights use the midpoint of asymmetric windows", {
    target <- 200.1234
    d <- S4Vectors::DataFrame(msLevel = c(1L, 2L, 1L), rtime = 1:3,
        scanIndex = 1:3, acquisitionNum = 1:3,
        precScanNum = c(NA_integer_, 1L, NA_integer_),
        precursorMz = c(NA_real_, target, NA_real_),
        precursorIntensity = c(NA_real_, 1000, NA_real_),
        isolationWindowTargetMz = c(NA_real_, target, NA_real_),
        isolationWindowLowerMz = c(NA_real_, 199.9234, NA_real_),
        isolationWindowUpperMz = c(NA_real_, 200.9234, NA_real_),
        dataOrigin = "asymmetric-window", centroided = TRUE)
    d$mz <- list(c(199.9734, target, 200.7234), c(57.0344, 70.0456),
        c(199.9734, target, 200.7234))
    d$intensity <- list(c(100, 1000, 100), c(10, 20), c(100, 500, 100))
    spectra <- Spectra::Spectra(d, source = Spectra::MsBackendMemory())
    initial <- assessPuritySingle(spectra, interpol = "none", iwNorm = TRUE,
        ilim = 0, isotopes = FALSE)

    # Analytic weights at exact knots of the default, rescaled Gaussian.
    weights <- (exp(-c(-2.7, -1.8, 1.8)^2 / 2) - exp(-9 / 2)) /
        (1 - exp(-9 / 2))
    ratios <- c(1000 * weights[2L] / sum(c(100, 1000, 100) * weights),
        500 * weights[2L] / sum(c(100, 500, 100) * weights))
    expect_equal(initial$aPurity, ratios[1L])
    # Moving or omitting the reported target does not move recorded bounds.
    for (target_label in c(target, 200.4234, NA_real_)) {
        spectra$isolationWindowTargetMz[2L] <- target_label
        pa <- purityA(spectra, iwNorm = TRUE, ilim = 0, isotopes = FALSE)
        settings <- processingParameters(pa, "purityA")[[1L]]$resolved[[1L]]
        expect_equal(purityTable(pa)$inPurity, mean(ratios),
            info = paste("target", target_label))
        expect_equal(settings$iwNormFun(c(-0.45, -0.30, 0.30)), weights,
            info = paste("target", target_label))
    }
})

test_that("symmetric defaults retain the existing Gaussian calculation", {
    spectra <- .window_spectra()[c(1L, 2L, 4L)]
    generated <- purityA(spectra, iwNorm = TRUE, isotopes = FALSE)
    explicit <- purityA(spectra, iwNorm = TRUE, isotopes = FALSE,
        iwNormFun = iwNormGauss(minOff = -0.4, maxOff = 0.4))
    expect_equal(purityTable(generated), purityTable(explicit))
})

test_that("an outside precursor has its own unavailable reason", {
    input <- .window_spectra()
    input$precursorMz[3L] <- 101.1234
    expect_warning(pa <- purityA(input, isotopes = FALSE, mostIntense = TRUE),
        "1 MS2 precursors are outside")
    settings <- processingParameters(pa, "purityA")[[1L]]$resolved[[1L]]
    expect_equal(purityTable(pa)$inPurity, c(0.5, NA_real_))
    expect_identical(settings$windows$status,
        c("valid", "precursor_outside_window"))
    expect_identical(settings$observations$method[2L],
        "precursor_outside_window")
    expect_identical(pa@spectra$msPurityPurityStatus[2L],
        "precursor_outside_window")
})
