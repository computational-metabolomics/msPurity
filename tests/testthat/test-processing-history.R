context("Processing history")

.d05_history_fixture <- function() {
    file <- system.file("extdata", "tests", "purityA",
        "2_frag4feature_pa.rds", package = "msPurity")
    pa <- updateObject(readRDS(file))
    # Restrict the workload to one real feature, retaining repeated links.
    pa@grped_df <- pa@grped_df[
        as.character(pa@grped_df$grpid) == "162", , drop = FALSE]
    pa@fragSpectra <- pa@fragSpectra[
        pa@fragSpectra$pid %in% pa@grped_df$pid]
    pa@grped_ms2 <- msPurity:::.pa_grouped_legacy(pa)
    pa
}

.d05_history_quiet <- function(expr) {
    utils::capture.output(result <- suppressMessages(expr))
    result
}

.d05_history_last <- function(pa, step) {
    records <- processingHistory(pa, step)
    records[[length(records)]]
}

.d05_history_average <- function(pa) {
    q <- .d05_history_quiet
    pa <- q(averageIntraFragSpectra(pa))
    pa <- q(averageInterFragSpectra(pa))
    q(averageAllFragSpectra(pa))
}

test_that("historical settings have stable incomplete identities", {
    file <- system.file("extdata", "tests", "purityA",
        "9_averageAllFragSpectra_with_filter_pa.rds", package = "msPurity")
    old <- readRDS(file)
    first <- processingHistory(old)
    expect_gt(length(first), 0)
    expect_identical(updateObject(old)@params$processing$version, 2L)
    expect_identical(processingHistory(old), first)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(old, path)
    expect_identical(processingHistory(readRDS(path)), first)
    for (record in first) {
        expect_identical(record$status, "imported")
        expect_false(record$complete)
        expect_true(is.na(record$time$started))
        expect_true(is.na(record$time$finished))
        expect_true(is.na(record$time$elapsed_seconds))
        expect_true(is.na(record$software$version))
        expect_identical(record$input$previous_steps, character())
    }
    expect_length(processingHistory(old, character()), 0)
    expect_error(processingHistory(old, "unknown-step"), "step must be")
})

test_that("new calls append distinct records and execution settings", {
    pa <- .d05_history_quiet(filterFragSpectra(.d05_history_fixture(),
            plim = 0.7, cores = 2))
    first <- processingHistory(pa)
    record <- .d05_history_last(pa, "filterFragSpectra")
    expect_identical(names(formals(processingHistory)), c("object", "step"))
    expect_identical(names(formals(processingParameters)), c("object", "step"))
    expect_identical(record$step, "filterFragSpectra")
    expect_identical(record$step_id, tail(names(first), 1))
    expect_false("stage" %in% names(record))
    expect_identical(record$status, "completed")
    expect_true(record$complete)
    expect_identical(record$parameters$cores, 2)
    expect_identical(record$parameters$effective_cores, 1L)
    expect_identical(record$software$version,
        as.character(utils::packageVersion("msPurity")))
    expect_s3_class(record$time$started, "POSIXct")
    expect_gte(record$time$elapsed_seconds, 0)
    expect_false(any(c("mz", "intensity") %in%
                names(record$input$spectra)))
    expect_false(any(c("mz", "intensity") %in%
                names(record$output$spectra)))
    parameters <- processingParameters(pa, step = "filterFragSpectra")
    expect_identical(parameters[[length(parameters)]], record$parameters)
    pa <- .d05_history_quiet(filterFragSpectra(pa, plim = 0.6))
    second <- processingHistory(pa)
    expect_identical(second[names(first)], first)
    expect_gt(length(second), length(first))
    expect_identical(anyDuplicated(names(second)), 0L)
    expect_length(processingHistory(pa, step = "filterFragSpectra"), 2)
})

test_that("earlier trial histories are rejected instead of migrated", {
    pa <- .d05_history_fixture()
    sp <- groupedSpectra(pa, legacy = FALSE)
    history <- pa@params$processing
    history$version <- 1L
    for (id in names(history$records)) {
        history$records[[id]]$stage <- history$records[[id]]$step
        history$records[[id]]$step <- NULL
    }
    pa@params$processing <- history
    metadata <- msPurity:::.pa_spectra_metadata(sp)
    metadata$msPurityProcessing <- history
    msPurity:::.pa_spectra_metadata(sp) <- metadata
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(pa, path)
    expect_error(processingHistory(readRDS(path)),
        "version 1 is from an earlier development trial")
    expect_error(processingHistory(sp),
        "version 1 is from an earlier development trial")
})

test_that("averaging records requested and actual workers separately", {
    for (actual in c(0, 1)) {
        pa <- .d05_history_fixture()
        pa@cores <- actual
        pa <- .d05_history_quiet(averageIntraFragSpectra(pa, cores = 2))
        record <- .d05_history_last(pa, "averageIntraFragSpectra")
        expect_identical(record$parameters$cores, 2)
        expect_identical(record$parameters$effective_cores, 1L)
        expect_identical(record$parameters$av, "median")
        expect_identical(pa@av_intra_params$cores, 2)
    }
})

test_that("refiltering invalidates dependent results but retains history", {
    q <- .d05_history_quiet
    pa <- .d05_history_average(q(filterFragSpectra(
        .d05_history_fixture(), plim = 0.7)))
    old <- processingHistory(pa)
    expect_gt(length(pa@avSpectra), 0)
    pa <- q(filterFragSpectra(pa, plim = 0.6))
    expect_length(pa@avSpectra, 0)
    expect_length(averagedSpectra(pa), 0)
    expect_length(pa@av_spectra, 0)
    expect_length(pa@av_intra_params, 0)
    expect_length(pa@av_inter_params, 0)
    expect_length(pa@av_all_params, 0)
    expect_identical(processingHistory(pa)[names(old)], old)
    expect_false(any(grepl("average", names(pa@params$processing$active))))
})

test_that("refusal after destructive filtering leaves results unchanged", {
    q <- .d05_history_quiet
    pa <- q(filterFragSpectra(.d05_history_fixture(), plim = 0.7,
            rmp = TRUE))
    pa <- .d05_history_average(pa)
    before <- processingHistory(pa)
    out <- q(filterFragSpectra(pa, plim = 0.6))
    expect_identical(out, pa)
    expect_identical(processingHistory(out), before)
})

test_that("intra-file recomputation retains independent all-member output", {
    q <- .d05_history_quiet
    pa <- .d05_history_average(q(filterFragSpectra(
        .d05_history_fixture(), plim = 0.7)))
    before <- pa@avSpectra[pa@avSpectra$av_level == "av_all"]
    history <- processingHistory(before)
    pa <- q(averageIntraFragSpectra(pa, ppm = 8))
    after <- pa@avSpectra[pa@avSpectra$av_level == "av_all"]
    expect_identical(Spectra::peaksData(after), Spectra::peaksData(before))
    expect_identical(Spectra::spectraData(after), Spectra::spectraData(before))
    expect_identical(processingHistory(after), history)
    expect_false(any(pa@avSpectra$av_level == "av_inter"))
    expect_length(pa@av_inter_params, 0)
})

test_that("history and function-valued settings survive extraction and RDS", {
    q <- .d05_history_quiet
    pa <- .d05_history_fixture()
    callback <- function(x) sqrt(x)
    environment(callback) <- baseenv()
    pa@params$purityA <- list(iwNormFun = callback, ppmInterp = 10)
    pa@params$processing <- NULL
    imported <- processingHistory(pa)
    expect_identical(processingHistory(pa), imported)
    expect_identical(processingParameters(pa, "purityA")[[1]]$iwNormFun,
        callback)
    pa <- q(filterFragSpectra(pa, plim = 0.7))
    sp <- groupedSpectra(pa, legacy = FALSE)
    metadata <- msPurity:::.pa_spectra_metadata(sp)
    expect_true("msPurityProcessing" %in% names(metadata))
    expect_gt(length(processingHistory(sp)), 0)
    expect_length(processingHistory(sp[integer()]), 0)
    one <- Spectra::setBackend(sp[1], Spectra::MsBackendMemory())
    expect_identical(processingHistory(one), processingHistory(sp[1]))
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(pa, path)
    restored <- readRDS(path)
    expect_identical(processingHistory(restored), processingHistory(pa))
    expect_equal(processingParameters(restored, "purityA")[[1]]$iwNormFun(4), 2)
    saveRDS(one, path)
    expect_identical(processingHistory(readRDS(path)), processingHistory(one))
})

test_that("invalid history and unresolved processing references fail", {
    pa <- .d05_history_fixture()
    pa <- .d05_history_quiet(filterFragSpectra(pa, plim = 0.7))
    broken <- pa
    broken@params$processing$active$filterFragSpectra <- "absent"
    expect_error(processingHistory(broken), "invalid step")
    broken <- pa
    id <- names(broken@params$processing$records)[1]
    broken@params$processing$records[[id]]$input$previous_steps <- "later"
    expect_error(processingHistory(broken), "unknown or later")
    sp <- groupedSpectra(pa, legacy = FALSE)
    sp$msPurityProcessingId[1] <- "absent"
    expect_error(processingHistory(sp), "absent from its history")
    expect_length(processingHistory(Spectra::Spectra(
        Spectra::MsBackendMemory())), 0)
})


test_that("constructor calls record measured input references and callbacks", {
    d <- S4Vectors::DataFrame(msLevel = c(1L, 2L, 1L),
        rtime = c(1, 2, 3), scanIndex = 1:3, acquisitionNum = 1:3,
        precScanNum = c(NA_integer_, 1L, NA_integer_),
        precursorMz = c(NA_real_, 100.1, NA_real_),
        dataOrigin = rep("history-memory", 3), centroided = TRUE)
    d$mz <- list(100.1, 70, 100.1)
    d$intensity <- list(100, 20, 90)
    sp <- Spectra::Spectra(d, source = Spectra::MsBackendMemory())
    callback <- function(x) x
    environment(callback) <- baseenv()
    pa <- .d05_history_quiet(purityA(sp, offsets = c(0.5, 0.5),
            interpol = "none", isotopes = FALSE, iwNormFun = callback))
    record <- .d05_history_last(pa, "purityA")
    expect_identical(record$status, "completed")
    expect_identical(record$parameters$iwNormFun, callback)
    expect_identical(record$parameters$effective_cores, 1L)
    expect_identical(nrow(record$input$spectra), 3L)
    expect_identical(nrow(record$output$spectra), 1L)
    expect_length(processingHistory(pa), 1)
    expect_identical(processingHistory(purityTable(pa, legacy = FALSE)),
        processingHistory(pa))
})

test_that("relinking invalidates filtered and averaged state", {
    skip_if_not_installed("xcms")
    loadNamespace("xcms")
    path <- system.file("extdata", "tests", "purityA",
        "2_frag4feature_pa.rds", package = "msPurity")
    raw <- updateObject(readRDS(path))
    original <- raw@fragSpectra
    local_mocked_bindings(
        .pa_raw_scans = function(pa, pids) {
            rows <- match(pids, original$pid)
            stopifnot(!anyNA(rows))
            original[rows]
        },
        .xcms_files = function(x) unname(as.character(raw@fileList)),
        .package = "msPurity")
    x <- readRDS(system.file("extdata", "tests", "xcms",
            "msms_only_xcmsnexp.rds", package = "msPurity"))
    q <- .d05_history_quiet
    pa <- .d05_history_average(q(filterFragSpectra(raw, plim = 0.7)))
    before <- processingHistory(pa)
    pa <- q(frag4feature(pa, x))
    expect_length(pa@filter_frag_params, 0)
    expect_identical(pa@all_frag_scans, methods::new("data.frame"))
    expect_length(pa@avSpectra, 0)
    expect_length(averagedSpectra(pa), 0)
    expect_false("pass_flag" %in% Spectra::peaksVariables(pa@fragSpectra))
    expect_false("purity_pass_flag" %in% names(purityTable(pa)))
    expect_identical(processingHistory(pa)[names(before)], before)
    expect_identical(.d05_history_last(pa, "frag4feature")$status, "completed")
    expect_no_error(q(filterFragSpectra(pa, plim = 0.7)))
})
