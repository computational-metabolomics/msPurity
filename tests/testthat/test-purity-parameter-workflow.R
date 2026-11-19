.purity_workflow_input <- function() {
    data <- S4Vectors::DataFrame(
        pid = c(100L, 7L, 8L, 101L, 40L, 41L, 60L, 102L),
        msLevel = c(1L, 2L, 2L, 1L, 3L, 3L, 4L, 1L),
        scanIndex = c(11L, 13L, 16L, 19L, 22L, 23L, 25L, 30L),
        acquisitionNum = c(101L, 104L, 109L, 112L, 117L, 121L, 127L, 130L),
        rtime = c(11, 13, 16, 19, 22, 23, 25, 30),
        precScanNum = c(NA_integer_, 101L, 101L, NA_integer_,
            104L, 104L, 117L, NA_integer_),
        precursorMz = c(NA_real_, 500, 500, NA_real_, 150, 150, 75, NA_real_),
        precursorIntensity = c(NA_real_, 100, 100, NA_real_,
            100, 100, 100, NA_real_),
        activationType = c(NA, "CID", "HCD", NA, "CID", "HCD", "ETD", NA),
        dataOrigin = "purity-parameter-workflow", centroided = TRUE)
    data$isolationWindowTargetMz <- data$precursorMz
    width <- c(NA_real_, 0.5, 0.5, NA_real_, 1.1, 1.1, 1.1, NA_real_)
    data$isolationWindowLowerMz <- data$precursorMz - width
    data$isolationWindowUpperMz <- data$precursorMz + width
    data$mz <- list(c(500, 500.2, 500.4), c(150, 150.2, 150.4, 151.003355),
        c(150, 150.2, 150.4, 151.003355), c(500, 500.2, 500.4),
        c(75, 75.2, 75.4, 76.003355), 50, 40, c(500, 500.2, 500.4))
    data$intensity <- list(c(100, 20, 5), c(100, 20, 5, 10),
        c(100, 20, 5, 10), c(100, 100, 5), c(100, 20, 5, 10),
        10, 10, c(100, 100, 5))
    Spectra::Spectra(data, source = Spectra::MsBackendMemory())
}

test_that("complete parameters reproduce scalar assessments", {
    input <- .purity_workflow_input()
    original <- Spectra::spectraData(input)
    scalar <- purityA(input, msLevel = 2:4,
        parentReference = "acquisitionNum", ilim = 0, isotopes = FALSE)
    typed <- purityA(input, msLevel = 2:4,
        parentReference = "acquisitionNum",
        param = PurityAParam(ilim = 0, isotopes = FALSE))

    expect_identical(purityTable(typed), purityTable(scalar))
    expect_equal(typed@spectra$inPurity, scalar@spectra$inPurity)
    expect_identical(typed@spectra$msPurityPurityStatus,
        scalar@spectra$msPurityPurityStatus)
    expect_identical(typed@spectra$pid, c(7L, 8L, 40L, 41L, 60L))
    expect_identical(Spectra::spectraData(input), original)
    expect_length(processingHistory(typed, "purityA"), 1L)
})

test_that("relative file paths retain parameter assignments", {
    input <- .purity_workflow_input()
    input$dataOrigin <- file.path(tempdir(), "input.mzML")
    param <- PurityAParam(ilim = 0, isotopes = FALSE)
    expected <- purityA(input, param = param)
    local_mocked_bindings(.msp_read = function(files) {
        expect_identical(files, "input.mzML")
        input
    }, .package = "msPurity")

    actual <- purityA("input.mzML", param = param)

    expect_identical(purityTable(actual), purityTable(expected))
})

test_that("MS2 groups without MS1 measurements do not read assessment peaks", {
    input <- .purity_workflow_input()
    input <- input[Spectra::msLevel(input) == 2L]
    params <- list(
        CID = PurityAParam(select = list(activationType = "CID"), ilim = 0),
        HCD = PurityAParam(select = list(activationType = "HCD"), ilim = 0.1))
    local_mocked_bindings(.msp_peaks = function(sp) {
        stop("Assessment must not read peaks without MS1 measurements.")
    }, .package = "msPurity")

    pa <- purityA(input, param = params)

    expect_identical(purityTable(pa)$seqNum, c(13L, 16L))
    expect_true(all(is.na(purityTable(pa)$inPurity)))
})

test_that("MS2 groups retain both MS1 neighbours and source positions", {
    input <- .purity_workflow_input()
    params <- list(
        CID = PurityAParam(select = list(msLevel = 2L, activationType = "CID"),
            ilim = 0, isotopes = FALSE),
        HCD = PurityAParam(select = list(msLevel = 2L, activationType = "HCD"),
            ilim = 0.1, isotopes = FALSE))
    reads <- list()
    read_peaks <- msPurity:::.msp_peaks
    local_mocked_bindings(.msp_peaks = function(sp) {
        reads[[length(reads) + 1L]] <<- Spectra::scanIndex(sp)
        read_peaks(sp)
    }, .package = "msPurity")

    pa <- purityA(input, param = params)
    table <- purityTable(pa)
    settings <- processingParameters(pa, "purityA")[[1L]]
    observations <- settings$resolved[[1L]]$observations

    # Events at 13 s and 16 s lie between MS1 measurements at 11 s and 19 s.
    expect_equal(table$inPurity, c(
        3 / 4 * (100 / 125) + 1 / 4 * (100 / 205),
        3 / 8 * (100 / 120) + 5 / 8 * (100 / 200)))
    expect_identical(table$seqNum, c(13L, 16L))
    expect_identical(observations$seqNum, c(13L, 16L))
    expect_identical(observations$ms1_acquisitionNum,
        list(c(101L, 112L), c(112L, 101L)))
    expect_identical(settings$resolved[[1L]]$windows$seqNum, c(13L, 16L))
    expect_identical(settings$configuration$assignments$parameter_set_names,
        list("CID", "HCD"))
    expect_identical(reads, list(Spectra::scanIndex(input)))
    expect_identical(Spectra::msLevel(pa@spectra), c(2L, 2L, 3L, 3L, 4L))
})

test_that("higher corrections follow child settings and share parent reads", {
    input <- .purity_workflow_input()
    # A +1.003355 peak at 10% lies within this explicit model's abundance range.
    isotope_model <- matrix(c(1, 1.003355, 10, 1, 0, 1, 75, 1), nrow = 1L)
    params <- list(
        low = PurityAParam(select = list(msLevel = 3L, activationType = "CID"),
            ilim = 0.075, isotopes = FALSE),
        isotope = PurityAParam(
            select = list(msLevel = 3L, activationType = "HCD"),
            ilim = 0, isotopes = TRUE, im = isotope_model),
        both = PurityAParam(select = list(msLevel = 4L, activationType = "ETD"),
            ilim = 0.075, isotopes = TRUE, im = isotope_model))
    reads <- list()
    read_peaks <- msPurity:::.msp_peaks
    local_mocked_bindings(.msp_peaks = function(sp) {
        reads[[length(reads) + 1L]] <<- Spectra::scanIndex(sp)
        read_peaks(sp)
    }, .package = "msPurity")

    pa <- purityA(input, msLevel = 3:4, parentReference = "acquisitionNum",
        param = params)
    settings <- processingParameters(pa, "purityA")[[1L]]
    observations <- settings$msn_purity$observations

    # Exclude the intensity-5 peak, the isotope-10 peak, or both, respectively.
    expect_equal(observations$purity, c(100 / 130, 100 / 125, 100 / 120))
    expect_identical(observations$parentPid, c(7L, 7L, 40L))
    expect_identical(observations$parent_peak_index, c(1L, 1L, 1L))
    expect_identical(settings$configuration$assignments$parameter_set_names,
        list("low", "isotope", "both"))
    expect_identical(reads, list(c(13L, 22L)))
    expect_equal(pa@spectra$inPurity,
        c(NA_real_, NA_real_, 100 / 130, 100 / 125, 100 / 120))
})

test_that("weighted MS2 and unweighted higher purity retain function history", {
    input <- .purity_workflow_input()
    weight <- function(x) ifelse(x > 0, 0.5, 1)
    environment(weight) <- baseenv()
    params <- list(
        weighted = PurityAParam(select = list(msLevel = 2L),
            ilim = 0, isotopes = FALSE, iwNorm = TRUE, iwNormFun = weight),
        higher = PurityAParam(select = list(msLevel = 3:4),
            ilim = 0, isotopes = FALSE))
    pa <- purityA(input, msLevel = 2:4, parentReference = "acquisitionNum",
        param = params)

    expect_equal(pa@spectra$inPurity, c(
        3 / 4 * (100 / 112.5) + 1 / 4 * (100 / 152.5),
        3 / 8 * (100 / 112.5) + 5 / 8 * (100 / 152.5),
        rep(100 / 135, 3L)))
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(pa, path)
    restored <- readRDS(path)
    settings <- processingParameters(restored, "purityA")[[1L]]
    saved <- settings$configuration$parameter_objects$weighted
    expect_identical(as.list(saved)$iwNormFun, weight)
    expect_equal(as.list(saved)$iwNormFun(c(-0.2, 0, 0.2)), c(1, 1, 0.5))
    expect_identical(processingHistory(restored@spectra[5L]),
        processingHistory(pa, "purityA"))
})

test_that("higher-level parameter assignments reject isolation weighting", {
    input <- .purity_workflow_input()
    params <- list(enabled = PurityAParam(iwNorm = TRUE),
        callback = PurityAParam(iwNormFun = function(x) rep(1, length(x))))
    for (name in names(params)) {
        expect_error(purityA(input, msLevel = 3:4,
                parentReference = "acquisitionNum", param = params[[name]]),
            "MS2", info = name)
    }
})

test_that("absent higher targets retain typed assessment observations", {
    input <- .purity_workflow_input()
    input <- input[Spectra::msLevel(input) != 4L]
    pa <- purityA(input, msLevel = 4L, parentReference = "acquisitionNum",
        param = PurityAParam(ilim = 0, isotopes = FALSE))
    record <- processingHistory(pa, "purityA")[[1L]]
    observations <- record$parameters$msn_purity$observations

    expect_s3_class(observations, "data.frame")
    expect_identical(nrow(observations), 0L)
    expect_identical(nrow(record$output$spectra), 0L)
    expect_identical(Spectra::msLevel(pa@spectra), c(2L, 2L, 3L, 3L))
})
