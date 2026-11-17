.feature_history_fixture <- function(adjusted = FALSE) {
    data <- S4Vectors::DataFrame(msLevel = c(1L, 1L, 2L, 1L),
        rtime = c(4, 6, 7, 8), acquisitionNum = 1:4, centroided = TRUE,
        scanIndex = 1:4,
        precScanNum = c(NA_integer_, NA_integer_, 2L, NA_integer_),
        precursorMz = c(NA_real_, NA_real_, 500.1234, NA_real_),
        precursorIntensity = c(NA_real_, NA_real_, 80, NA_real_),
        dataOrigin = "feature-history-memory")
    data$mz <- list(c(500.1234, 500.3234), c(500.1234, 500.3234), 150.1234,
        c(500.1234, 500.3234))
    data$intensity <- list(c(80, 20), c(80, 20), 30, c(80, 20))
    pa <- purityA(Spectra::Spectra(data, source = Spectra::MsBackendMemory()),
        offsets = c(0.5, 0.5), isotopes = FALSE, interpol = "none")
    peak <- c(mz = 500.1234, mzmin = 500.1233, mzmax = 500.1235,
        rt = 6, rtmin = 4, rtmax = 8, into = 100, maxo = 80, sample = 1)
    group <- c(mzmed = 500.1234, mzmin = 500.1233, mzmax = 500.1235,
        rtmed = 6, rtmin = 4, rtmax = 8, npeaks = 1, sample = 1)
    xe <- methods::new("xcmsSet", peaks = rbind(peak), groups = rbind(group),
        groupidx = list(1L), filepaths = unname(pa@fileList),
        phenoData = data.frame(class = "sample", row.names = "sample"),
        rt = list(raw = list(c(4, 6, 8)), corrected = list(c(4, 6, 8))))
    if (adjusted) {
        xe@rt$corrected <- list(c(5, 7, 9))
        xe@peaks[, c("rt", "rtmin", "rtmax")] <- c(7, 5, 9)
        xe@groups[, c("rtmed", "rtmin", "rtmax")] <- c(7, 5, 9)
        xe@.processHistory <- list(methods::new("XProcessHistory",
                type = "Retention time correction"))
    }
    list(pa = pa, xe = xe)
}

.feature_history_record <- function(pa) {
    tail(processingHistory(pa, "frag4feature"), 1L)[[1L]]
}

test_that("association records actual serial, parallel and grouped workers", {
    skip_if_not_installed("xcms")
    input <- .feature_history_fixture()
    serial <- frag4feature(input$pa, input$xe, convert2RawRT = FALSE)
    input$pa@cores <- 2
    # plyr emits this compiler warning on some R versions, including 4.6.
    parallel <- withCallingHandlers(
        frag4feature(input$pa, input$xe, convert2RawRT = FALSE),
        warning = function(warning) {
            known <- paste0("<anonymous>: ... may be used in an incorrect ",
                "context: ", sQuote(".fun(piece, ...)"))
            if (identical(trimws(conditionMessage(warning)), known))
                invokeRestart("muffleWarning")
        })
    grouped <- frag4feature(input$pa, input$xe, convert2RawRT = FALSE,
        useGroup = TRUE)
    cases <- list(serial = serial, parallel = parallel, grouped = grouped)
    expected <- list(serial = c(1, 1), parallel = c(2, 2), grouped = c(2, 1))
    for (name in names(cases)) {
        record <- .feature_history_record(cases[[name]])
        expect_equal(unname(unlist(record$parameters[
            c("cores", "effective_cores")])),
        expected[[name]], info = name)
        expect_identical(record$parameters$rt_conversion,
            list(preparation = "not_requested", used_for_association = FALSE),
            info = name)
    }
    expect_identical(parallel@grped_df, serial@grped_df)
    reference <- .feature_history_record(serial)$input$xcms
    expect_identical(reference$class, "xcmsSet")
    expect_identical(reference$files, unname(input$pa@fileList))
    expect_identical(reference$association_input, "individual_peaks")
    expect_identical(.feature_history_record(grouped)$input$xcms$
            association_input, "group_bounds")
})

test_that("history distinguishes prepared and consumed RT conversion", {
    skip_if_not_installed("xcms")
    raw <- .feature_history_fixture()
    unadjusted <- frag4feature(raw$pa, raw$xe)
    aligned <- .feature_history_fixture(adjusted = TRUE)
    converted <- frag4feature(aligned$pa, aligned$xe)
    grouped <- frag4feature(aligned$pa, aligned$xe, useGroup = TRUE)
    cases <- list(unadjusted = unadjusted, converted = converted,
        grouped = grouped)
    expected <- list(unadjusted = list(preparation = "no_adjustment",
            used_for_association = FALSE),
        converted = list(preparation = "applied", used_for_association = TRUE),
        grouped = list(preparation = "applied", used_for_association = FALSE))
    for (name in names(cases))
        expect_identical(.feature_history_record(cases[[name]])$parameters$
                rt_conversion, expected[[name]], info = name)
    expect_equal(converted@grped_df$rtmin, 4)
    expect_equal(grouped@grped_df$rtmin_full, 5)
    expect_equal(grouped@grped_df$rtmax_full, 9)
})

test_that("empty feature inputs record no execution or fingerprint", {
    skip_if_not_installed("xcms")
    input <- .feature_history_fixture()
    input$pa@cores <- 2
    input$xe@peaks <- input$xe@peaks[FALSE, , drop = FALSE]
    out <- frag4feature(input$pa, input$xe)
    record <- .feature_history_record(out)
    expect_length(out@fragSpectra, 0L)
    expect_equal(record$parameters$effective_cores, 0)
    expect_identical(record$parameters$rt_conversion,
        list(preparation = "no_chromatographic_peaks",
            used_for_association = FALSE))
    expect_identical(record$input$xcms$association_input, "not_used")
    expect_null(record$input$xcms$fingerprint)
})

test_that("feature fingerprints retain prepared bounds and memberships", {
    skip_if_not_installed("xcms")
    input <- .feature_history_fixture()
    link <- function(xe) frag4feature(input$pa, xe, convert2RawRT = FALSE)
    fingerprint <- function(out) {
        .feature_history_record(out)$input$xcms$fingerprint
    }
    original <- link(input$xe)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(input$xe, path)
    expect_identical(fingerprint(link(readRDS(path))), fingerprint(original))
    saveRDS(original, path)
    expect_identical(fingerprint(readRDS(path)), fingerprint(original))
    expect_identical(fingerprint(original)[c("schema", "serialization_version",
                "algorithm")], list(schema = "frag4feature-input-v1",
            serialization_version = 2L, algorithm = "sha256"))
    expect_match(fingerprint(original)$value, "^[a-f0-9]{64}$")
    changed <- input$xe
    changed@peaks[1L, "rtmin"] <- 3
    expect_false(identical(fingerprint(link(changed)), fingerprint(original)))
    changed <- input$xe
    changed@groupidx <- list(c(1L, 1L))
    expect_false(identical(fingerprint(link(changed)), fingerprint(original)))
    changed@groupidx <- list(c(1L, 1L), 1L)
    ordered <- fingerprint(link(changed))
    changed@groupidx <- rev(changed@groupidx)
    expect_false(identical(fingerprint(link(changed)), ordered))
})
