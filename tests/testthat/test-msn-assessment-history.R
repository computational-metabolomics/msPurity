.msn_assessment_history_input <- function() {
    data <- S4Vectors::DataFrame(
        pid = c(100L, 7L, 40L, 41L, 101L),
        msLevel = c(1L, 2L, 3L, 4L, 1L),
        scanIndex = 1:5, acquisitionNum = 11:15, rtime = 1:5,
        spectrumId = paste0("scan=", 11:15),
        precScanNum = c(NA_integer_, 11L, 12L, 13L, NA_integer_),
        precursorMz = c(NA_real_, 500, 150, 75, NA_real_),
        precursorIntensity = c(NA_real_, 80, 80, 30, NA_real_),
        dataOrigin = "assessment-history-memory", centroided = TRUE)
    data$mz <- list(c(500, 500.2), c(150, 150.2), c(75, 75.2),
        50, c(500, 500.2))
    data$intensity <- list(c(80, 20), c(80, 20), c(30, 70),
        10, c(80, 20))
    Spectra::Spectra(data, source = Spectra::MsBackendMemory())
}

test_that("constructor retention preserves supplied fragmentation pids", {
    input <- .msn_assessment_history_input()
    original <- Spectra::spectraData(input)
    pa <- purityA(input, offsets = c(0.5, 0.5),
        interpol = "none", ilim = 0, isotopes = FALSE)

    expect_identical(pa@spectra$pid, c(7L, 40L, 41L))
    expect_identical(purityTable(pa)$pid, 7L)
    expect_equal(purityTable(pa)$aPurity, 0.8)
    expect_identical(Spectra::spectraData(input), original)

    input$pid <- NULL
    generated <- purityA(input, offsets = c(0.5, 0.5),
        interpol = "none", ilim = 0, isotopes = FALSE)
    expect_equal(generated@spectra$pid, 1:3)
    expect_equal(purityTable(generated)$pid, 1L)
})

test_that("retention rejects invalid supplied fragmentation identities", {
    input <- .msn_assessment_history_input()
    input$pid[3L] <- input$pid[2L]
    expect_error(purityA(input, msLevel = 4L,
            parentReference = "acquisitionNum", offsets = c(0.5, 0.5)),
        "pids must be unique finite")
})

test_that("history retains each requested assessment", {
    input <- .msn_assessment_history_input()
    callback <- function(x) rep(1, length(x))
    environment(callback) <- baseenv()
    first <- purityA(input, msLevel = 3L,
        parentReference = "acquisitionNum", offsets = c(0.5, 0.5),
        iwNormFun = callback)
    earlier <- processingHistory(first)
    first_id <- names(earlier)[1L]

    expect_length(earlier, 1L)
    expect_identical(first@spectra$msPurityProcessingId,
        c(NA_character_, first_id, NA_character_))
    expect_identical(earlier[[1L]]$output$spectra$pid, 40L)
    expect_true(all(is.na(purityTable(first)$inPurity)))
    expect_identical(first@spectra$msPurityPurityStatus[1L], "not_assessed")

    second <- purityA(first@spectra, msLevel = 4L,
        parentReference = "acquisitionNum", offsets = c(0.5, 0.5))
    history <- processingHistory(second)
    second_id <- names(history)[2L]

    expect_length(history, 2L)
    expect_identical(history[names(earlier)], earlier)
    expect_identical(second@spectra$pid, c(7L, 40L, 41L))
    expect_identical(second@spectra$msPurityProcessingId,
        c(NA_character_, first_id, second_id))
    expect_equal(second@spectra$inPurity, c(NA_real_, 0.8, 0.3))
    expect_identical(second@spectra$msPurityPurityStatus,
        c("not_assessed", "assessed", "assessed"))
    expect_identical(history[[2L]]$input$previous_steps, first_id)
    expect_identical(history[[2L]]$output$spectra$pid, 41L)
    expect_identical(processingHistory(second@spectra[2L]), earlier)
    expect_identical(processingHistory(second@spectra[3L]), history)
    expect_length(processingHistory(second@spectra[1L]), 0L)

    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(second, path)
    restored <- readRDS(path)
    expect_identical(processingHistory(restored), history)
    expect_identical(processingHistory(restored@spectra[3L]), history)
    expect_identical(processingParameters(restored, "purityA")[[1L]]$iwNormFun,
        callback)
})

test_that("source peak identities survive reassessment after reordering", {
    first <- purityA(.msn_assessment_history_input(), msLevel = 2:4,
        parentReference = "acquisitionNum", offsets = c(0.5, 0.5))
    input <- first@spectra[c(3L, 2L)]
    second <- purityA(input, msLevel = 4L,
        parentReference = "acquisitionNum", offsets = c(0.5, 0.5))

    expect_identical(second@spectra$pid, c(41L, 40L))
    expect_equal(second@spectra$msPurityRawPeakOffset, c(4, 2))
    expect_identical(second@spectra$msPurityRawPeakCount, c(1L, 2L))
    expect_identical(second@spectra$msPurityRawPeakIdBasis,
        c("source_reference", "source_reference"))
    positions <- list(original_positions = list(1L, 1:2),
        stored_positions = list(1L, 1:2))
    expect_equal(msPurity:::.pa_source_peak_ids(second@spectra, positions),
        list(5, 3:4))
    expect_equal(second@spectra$inPurity, c(0.3, 0.8))
})

test_that("purity reassessment retains existing flag references", {
    pa <- purityA(.msn_assessment_history_input(), msLevel = 2:4,
        parentReference = "acquisitionNum", offsets = c(0.5, 0.5))
    pa <- filterFragSpectra(pa, allfrag = TRUE, msLevel = 2:4, snr = 0)
    input <- pa@spectra
    filtering <- processingHistory(input, "filterFragSpectra")
    out <- purityA(input, msLevel = 3L, parentReference = "acquisitionNum",
        offsets = c(0.5, 0.5))
    expect_identical(out@spectra$purity_pass_flag, input$purity_pass_flag)
    expect_identical(out@spectra$msPurityFilterProcessingId,
        input$msPurityFilterProcessingId)
    expect_identical(processingHistory(out@spectra, "filterFragSpectra"),
        filtering)
})

test_that("incomplete source peak provenance is rejected", {
    input <- .msn_assessment_history_input()
    input$msPurityRawPeakOffset <- rep(0, length(input))
    expect_error(purityA(input, msLevel = 4L,
            parentReference = "acquisitionNum", offsets = c(0.5, 0.5)),
        "source peak provenance is incomplete")
})

test_that("assessment with no requested spectra records no output tips", {
    input <- .msn_assessment_history_input()
    input <- input[Spectra::msLevel(input) != 4L]
    pa <- purityA(input, msLevel = 4L,
        parentReference = "acquisitionNum", offsets = c(0.5, 0.5))
    history <- processingHistory(pa)

    expect_length(history, 1L)
    expect_identical(nrow(history[[1L]]$output$spectra), 0L)
    expect_true(all(is.na(pa@spectra$msPurityProcessingId)))
    expect_length(processingHistory(pa@spectra), 0L)
})

test_that("higher-level assessment needs its parent but no MS2 feature root", {
    input <- .msn_assessment_history_input()
    input <- input[Spectra::msLevel(input) %in% 3:4]
    pa <- purityA(input, msLevel = 4L, parentReference = "acquisitionNum",
        offsets = c(0.5, 0.5))
    expect_equal(pa@spectra$inPurity, c(NA_real_, 0.3))
    expect_true(all(is.na(pa@params$msnLineage$scans$rootPid)))
    expect_identical(pa@spectra$msPurityPurityStatus,
        c("not_assessed", "assessed"))

    queued <- Spectra::filterIntensity(input, intensity = c(50, Inf))
    processed <- purityA(queued, msLevel = 4L,
        parentReference = "acquisitionNum", offsets = c(0.5, 0.5))
    expect_equal(processed@spectra$inPurity, c(NA_real_, 1))
    expect_identical(lengths(queued), c(1L, 0L))
})

test_that("MS2-only assessment does not resolve unrelated higher parents", {
    local_mocked_bindings(.pa_msn_preparation = function(...) {
        stop("MS2-only assessment must not resolve higher-level ancestry")
    }, .package = "msPurity")
    input <- .msn_assessment_history_input()
    pa <- purityA(input, offsets = c(0.5, 0.5), isotopes = FALSE)
    expect_equal(purityTable(pa)$inPurity, 0.8)
    expect_null(pa@params$msnLineage)
})
