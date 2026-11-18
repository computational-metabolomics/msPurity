context("Independent MSn product-quality filtering")

.msn_filter_fixture <- function() {
    data <- S4Vectors::DataFrame(msLevel = c(1L, 2L, 3L, 4L, 3L, 4L),
        scanIndex = 1:6, acquisitionNum = 100:105, rtime = as.numeric(1:6),
        precursorMz = c(NA, 500, 100, 75, 500, 100),
        precScanNum = c(NA_integer_, 100:104), polarity = 1L,
        dataOrigin = "synthetic", spectrumId = paste0("scan=", 100:105))
    data$mz <- list(500, c(100, 150), c(75, 120), c(50, 60),
        c(20, 40), numeric())
    data$intensity <- list(1000, c(20, 10), c(100, 5), c(1, 50),
        c(30, 3), numeric())
    raw <- Spectra::Spectra(data, source = Spectra::MsBackendMemory())
    purity <- data.frame(pid = 1L, fileid = 1L, seqNum = 2L,
        acquisitionNum = 101L, retentionTime = 2, precursorMZ = 500,
        precursorScanNum = 100L, inPurity = 0.1)
    pa <- methods::new("purityA", fileList = c(synthetic = "synthetic"),
        cores = 1)
    pa@spectra <- msPurity:::.pa_msn_retained_spectra(raw, rep(1L, 6),
        purity)
    pa@grped_df <- data.frame(grpid = 1L, pid = 1:3, sample = 1L,
        cid = 1:3, msLevel = 2:4, inPurity = c(0.1, NA, NA),
        precurMtchMZ = c(500, 100, 75), retentionTime = as.numeric(2:4))
    peaks <- as.list(Spectra::peaksData(pa@spectra))
    linked <- list(`1` = peaks[1:3])
    pa@fragSpectra <- msPurity:::.pa_frag_from_grouped(pa, linked)
    msPurity:::.pa_sync_legacy(pa, puritydf = purity, grped_ms2 = linked)
}

.msn_filter_call <- function(pa, typed, msLevel = 2:4, ...) {
    if (typed) filterFragSpectra(pa, ilim = 10, snr = 0,
        param = FilterFragSpectraParam(select = list(msLevel = 2L),
            plim = 0.8, ilim = 10, snr = 0), msLevel = msLevel, ...)
    else filterFragSpectra(pa, ilim = 10, snr = 0, msLevel = msLevel, ...)
}

.msn_assessed_filter_fixture <- function() {
    pa <- .msn_filter_fixture()
    pa@spectra$msPurityStagePurity <- c(0.1, 0.8, 0.9, NA, NA)
    pa@spectra$msPurityPurityStatus <- c(rep("assessed", 3L),
        "missing_parent", "not_assessed")
    pa@spectra$activationType <- c("CID", "HCD", "CID", "CID", "HCD")
    for (name in c("msPurityStagePurity", "msPurityPurityStatus",
                "activationType"))
        pa@fragSpectra[[name]] <- pa@spectra[[name]][1:3]
    pa
}

test_that("purity access and updates remain MS2 only with retained MSn", {
    pa <- .msn_filter_fixture()
    expect_identical(Spectra::msLevel(purityTable(pa, legacy = FALSE)), 2L)
    expect_identical(purityTable(pa)$pid, 1L)
    native <- pa@spectra$spectrumId
    changed <- purityTable(pa)
    changed$inPurity <- 0.5
    changed$purity_pass_flag <- FALSE
    pa@spectra <- msPurity:::.pa_set_purity_vars(pa@spectra, changed)
    expect_identical(pa@spectra$pid, 1:5)
    expect_identical(pa@spectra$spectrumId, native)
    expect_identical(pa@spectra$inPurity, c(0.5, rep(NA_real_, 4)))
    expect_true(all(is.na(pa@spectra$purity_pass_flag[-1])))
})

test_that("MS2 purity failure does not reject higher product-quality peaks", {
    old <- options(msPurity.legacySlots = TRUE)
    on.exit(options(old))
    for (legacy in c(TRUE, FALSE)) for (typed in c(TRUE, FALSE)) {
        options(msPurity.legacySlots = legacy)
        case <- paste("legacy:", legacy, "typed:", typed)
        out <- .msn_filter_call(.msn_filter_fixture(), typed)
        spectra <- groupedSpectra(out, legacy = FALSE)
        peaks <- as.list(Spectra::peaksData(spectra, columns = "pass_flag"))
        expect_identical(spectra$purity_pass_flag, c(FALSE, TRUE, TRUE),
            info = case)
        expect_identical(spectra$msPurityPurityStatus,
            c("assessed", "not_assessed", "not_assessed"), info = case)
        expect_identical(lapply(peaks, function(m) unname(as.matrix(m))),
            lapply(list(c(0, 0), c(1, 0), c(0, 1)), matrix, ncol = 1),
            info = case)
        groups <- groupedSpectra(out, legacy = TRUE)[["1"]]
        expect_true(all(groups[[2]][, "purity_pass_flag"] == 1),
            info = case)
        expect_identical(unname(groups[[2]][, "pass_flag"]), c(1, 0),
            info = case)
    }
})

test_that("removal retains higher-level peaks at their original positions", {
    for (typed in c(TRUE, FALSE)) {
        out <- .msn_filter_call(.msn_filter_fixture(), typed, rmp = TRUE)
        spectra <- groupedSpectra(out, legacy = FALSE)
        expect_identical(lengths(Spectra::mz(spectra)), c(0L, 1L, 1L),
            info = paste("typed:", typed))
        peaks <- as.list(Spectra::peaksData(spectra,
                columns = "original_peak_index"))
        expect_identical(unname(peaks[[2]][, 1]), 1)
        expect_identical(unname(peaks[[3]][, 1]), 2)
    }
})

test_that("allfrag includes unlinked and empty higher spectra", {
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    pa <- .msn_filter_fixture()
    pa@fragSpectra <- Spectra::addProcessing(pa@fragSpectra, function(x) {
        stop("allfrag must use source peaks, not reprocess linked peak arrays")
    })
    for (typed in c(TRUE, FALSE)) {
        case <- paste("typed:", typed)
        out <- .msn_filter_call(pa, typed, allfrag = TRUE, rmp = TRUE)
        spectra <- allFragSpectra(out, legacy = FALSE)
        expect_identical(spectra$pid, 1:5, info = case)
        expect_identical(lengths(Spectra::mz(spectra)),
            c(0L, 1L, 1L, 1L, 0L), info = case)
        expect_identical(spectra$purity_pass_flag,
            c(FALSE, rep(TRUE, 4)), info = case)
        id <- tail(names(processingHistory(out, "filterFragSpectra")), 1L)
        expect_identical(spectra$msPurityFilterProcessingId,
            rep(id, 5L), info = case)
        table <- allFragSpectra(out, legacy = TRUE)
        expect_identical(table$pid, 2:4, info = case)
        expect_identical(table$sid, c(3L, 6L, 7L), info = case)
        expect_false(anyNA(table$scan), info = case)
        expect_true(all(table$purity_pass_flag == 1), info = case)
        expect_equal(nrow(purityTable(out)), 1L, info = case)
    }
})

test_that("an empty linked higher spectrum remains empty after filtering", {
    pa <- .msn_filter_fixture()
    empty <- pa@grped_df[3, , drop = FALSE]
    empty$pid <- 5L
    empty$cid <- 4L
    empty$retentionTime <- 6
    pa@grped_df <- rbind(pa@grped_df, empty)
    peaks <- as.list(Spectra::peaksData(pa@spectra))
    pa@fragSpectra <- msPurity:::.pa_frag_from_grouped(pa,
        list(`1` = peaks[pa@grped_df$pid]))
    for (typed in c(TRUE, FALSE)) {
        case <- paste("typed:", typed)
        out <- .msn_filter_call(pa, typed)
        spectra <- groupedSpectra(out, legacy = FALSE)
        expect_identical(spectra$pid, c(1:3, 5L), info = case)
        expect_identical(lengths(Spectra::mz(spectra)), c(2L, 2L, 2L, 0L),
            info = case)
        expect_true(spectra$purity_pass_flag[4], info = case)
        expect_identical(length(groupedSpectra(out,
                    legacy = TRUE)[["1"]][[4]]), 0L, info = case)
    }
})

test_that("existing scalar averaging cannot silently mix MSn levels", {
    pa <- .msn_filter_fixture()
    expect_error(averageAllFragSpectra(pa),
        "Averaging MS3/MS4 feature links is not supported")
})

test_that("processed higher levels cannot enter the MS2-only writers", {
    higher <- .msn_filter_fixture()@spectra
    higher <- higher[Spectra::msLevel(higher) > 2L]
    for (field in c("fragSpectra", "avSpectra")) {
        pa <- methods::new("purityA")
        methods::slot(pa, field) <- higher
        expect_error(create_database(pa, NULL),
            "Processed MS3/MS4 database export is not supported", info = field)
        expect_error(createDatabase(pa, NULL, format = "sqlite"),
            "Processed MS3/MS4 database export is not supported", info = field)
    }
    ## The Parquet writer takes filtered MS3/MS4 spectra, not averages,
    ## which are not produced yet.
    pa <- methods::new("purityA")
    pa@avSpectra <- higher
    expect_error(createDatabase(pa, NULL, format = "parquet"),
        "Averaged MS3/MS4 spectra cannot be exported yet")
    expect_no_error(msPurity:::.pa_require_msn_export(
        methods::new("purityA", fragSpectra = higher)))
})

test_that("empty selections retain empty filtering history", {
    for (allfrag in c(FALSE, TRUE)) {
        out <- filterFragSpectra(methods::new("purityA"), allfrag = allfrag)
        expect_identical(length(out@fragSpectra), 0L,
            info = paste("allfrag:", allfrag))
        expect_equal(nrow(purityTable(out)), 0L)
        history <- processingHistory(out, "filterFragSpectra")
        expect_length(history, 1L)
        settings <- processingParameters(out, "filterFragSpectra")[[1L]]
        expect_identical(settings$purity_filter,
            data.frame(pid = integer(), threshold = numeric(),
                decision = character()))
    }
})

test_that("level selection controls spectra, feature occurrences and history", {
    pa <- .msn_filter_fixture()
    pa@grped_df <- pa@grped_df[c(1L, 2L, 2L, 3L), ]
    cases <- list(
        linked_ms2 = list(all = FALSE, levels = 2L, pids = 1L),
        linked_higher = list(all = FALSE, levels = 3:4, pids = 2:3),
        linked_mixed = list(all = FALSE, levels = c(4, 2, 4), pids = c(1L, 3L)),
        all_higher = list(all = TRUE, levels = 3:4, pids = 2:5),
        all_levels = list(all = TRUE, levels = 2:4, pids = 1:5))
    for (typed in c(FALSE, TRUE)) for (name in names(cases)) {
        case <- cases[[name]]
        label <- paste(name, "typed:", typed)
        out <- .msn_filter_call(pa, typed, allfrag = case$all,
            msLevel = case$levels)
        record <- tail(processingHistory(out, "filterFragSpectra"), 1L)[[1L]]
        expect_identical(out@fragSpectra$pid, case$pids, info = label)
        expect_identical(out@grped_df$pid,
            pa@grped_df$pid[pa@grped_df$pid %in% case$pids], info = label)
        expect_identical(record$input$spectra$pid, case$pids, info = label)
        expect_identical(record$output$spectra$pid, case$pids, info = label)
        expect_identical(record$parameters$msLevel, case$levels, info = label)
        expect_identical(record$parameters$effective_msLevel,
            unique(Spectra::msLevel(out@fragSpectra)), info = label)
        if (case$all && !2L %in% case$levels)
            expect_true(all(is.na(purityTable(out)$purity_pass_flag)),
                info = label)
    }
})

test_that("higher-level filtering clears stale MS2 flag references", {
    all <- filterFragSpectra(.msn_filter_fixture(), allfrag = TRUE,
        msLevel = 2:4, snr = 0)
    previous <- processingHistory(all)
    higher <- filterFragSpectra(all, allfrag = TRUE, msLevel = 3:4, snr = 0)
    raw_ms2 <- purityTable(higher, legacy = FALSE)
    expect_true(is.na(raw_ms2$purity_pass_flag))
    expect_true(is.na(raw_ms2$msPurityFilterProcessingId))
    expect_length(processingHistory(raw_ms2, "filterFragSpectra"), 0L)
    current <- tail(processingHistory(higher, "filterFragSpectra"), 1L)[[1L]]
    expect_identical(higher@fragSpectra$msPurityFilterProcessingId,
        rep(current$step_id, length(higher@fragSpectra)))
    expect_identical(processingHistory(higher)[names(previous)], previous)
})

test_that("default MS2 selection does not read excluded peak arrays", {
    pa <- .msn_filter_fixture()
    pa@fragSpectra <- Spectra::addProcessing(pa@fragSpectra, function(x) {
        if (nrow(x) && x[1L, "mz"] != 100)
            stop("An excluded higher-level spectrum was read")
        x
    })
    for (typed in c(FALSE, TRUE)) {
        args <- if (typed) list(param = FilterFragSpectraParam(snr = 0)) else
            list(snr = 0)
        implicit <- do.call(filterFragSpectra, c(list(pa = pa), args))
        explicit <- do.call(filterFragSpectra,
            c(list(pa = pa, msLevel = 2L), args))
        expect_identical(implicit@fragSpectra$pid, 1L)
        expect_identical(Spectra::peaksData(implicit@fragSpectra),
            Spectra::peaksData(explicit@fragSpectra))
    }
    generic <- methods::getGeneric("filterFragSpectra")
    expect_identical(tail(names(formals(generic)), 2L), c("param", "msLevel"))
})

test_that("parameter selectors do not extend the selected levels", {
    pa <- .msn_filter_fixture()
    out <- filterFragSpectra(pa, allfrag = TRUE,
        param = list(higher = FilterFragSpectraParam(
            select = list(msLevel = 3:4), snr = 0)))
    expect_identical(allFragSpectra(out, legacy = FALSE)$pid, 1L)
    settings <- tail(processingParameters(out, "filterFragSpectra"), 1L)[[1L]]
    expect_identical(settings$configuration$unused_parameter_sets, "higher")
    expect_identical(settings$configuration$assignments$fallback_used, TRUE)
})

test_that("narrowing levels clears stale outputs and retains earlier history", {
    old <- options(msPurity.legacySlots = TRUE)
    on.exit(options(old))
    for (typed in c(FALSE, TRUE)) {
        all <- .msn_filter_call(.msn_filter_fixture(), typed, allfrag = TRUE)
        previous <- processingHistory(all)
        selected <- .msn_filter_call(all, typed, msLevel = 3L)
        expect_identical(selected@grped_df$pid, 2L)
        expect_equal(nrow(selected@all_frag_scans), 0L)
        expect_length(allFragSpectra(selected, legacy = FALSE), 0L)
        expect_identical(processingHistory(selected)[names(previous)], previous)
        empty <- .msn_filter_call(selected, typed, msLevel = 4L)
        expect_length(groupedSpectra(empty, legacy = FALSE), 0L)
        expect_equal(nrow(empty@grped_df), 0L)
        expect_identical(empty@spectra$pid, 1:5)
    }
})

test_that("invalid levels and destructive refiltering preserve inputs", {
    pa <- .msn_filter_fixture()
    before <- serialize(pa, NULL)
    for (levels in list(NULL, integer(), NA, Inf, NaN, TRUE, "2", 1, 5, 2.5))
        expect_error(filterFragSpectra(pa, msLevel = levels), "msLevel",
            info = paste("levels:", paste(levels, collapse = ",")))
    expect_identical(serialize(pa, NULL), before)
    removed <- .msn_filter_call(pa, TRUE, msLevel = 2L, rmp = TRUE)
    expect_message(unchanged <- filterFragSpectra(removed, msLevel = 3:4),
        "previously filtered and removed")
    expect_identical(processingHistory(unchanged), processingHistory(removed))
    expect_identical(unchanged@fragSpectra$pid, 1L)
})
test_that("unassigned higher-level purity thresholds remain disabled", {
    pa <- .msn_filter_fixture()
    pa@spectra$msPurityStagePurity <- c(0.1, 0.01, 0.9, NA, NA)
    pa@spectra$msPurityPurityStatus <- c("assessed", "assessed", "assessed",
        "not_assessed", "not_assessed")
    pa@fragSpectra$msPurityStagePurity <- c(0.1, 0.01, 0.9)
    pa@fragSpectra$msPurityPurityStatus <- rep("assessed", 3L)
    for (typed in c(FALSE, TRUE)) for (all in c(FALSE, TRUE)) {
        label <- paste("typed:", typed, "allfrag:", all)
        baseline <- .msn_filter_call(.msn_filter_fixture(), typed,
            allfrag = all, rmp = TRUE)
        out <- .msn_filter_call(pa, typed, allfrag = all, rmp = TRUE)
        expect_identical(Spectra::peaksData(out@fragSpectra),
            Spectra::peaksData(baseline@fragSpectra), info = label)
        expected <- if (all) c(0.1, 0.01, 0.9, NA, NA) else c(0.1, 0.01, 0.9)
        expect_identical(out@fragSpectra$msPurityStagePurity, expected,
            info = label)
        expect_true(all(out@fragSpectra$purity_pass_flag[-1L]),
            info = label)
    }
})

test_that("explicit thresholds use each level's own purity and peak flags", {
    pa <- .msn_assessed_filter_fixture()
    cases <- list(
        scalar_linked_flags = list(typed = FALSE, all = FALSE, remove = FALSE),
        typed_all_flags = list(typed = TRUE, all = TRUE, remove = FALSE),
        scalar_all_removed = list(typed = FALSE, all = TRUE, remove = TRUE),
        typed_linked_removed = list(typed = TRUE, all = FALSE, remove = TRUE))
    for (name in names(cases)) {
        case <- cases[[name]]
        settings <- if (case$typed)
            list(param = FilterFragSpectraParam(plim = 0.8, ilim = 10,
                    snr = 0)) else list(plim = 0.8, ilim = 10, snr = 0)
        out <- do.call(filterFragSpectra, c(list(pa = pa, msLevel = 2:4,
                    allfrag = case$all, rmp = case$remove), settings))
        spectra <- out@fragSpectra
        ids <- if (case$all) 1:5 else 1:3
        expect_identical(spectra$purity_pass_flag,
            c(FALSE, FALSE, TRUE, FALSE, FALSE)[ids], info = name)
        expect_identical(spectra$msPurityFilterDecision,
            c("failed", "failed", "passed", "missing", "missing")[ids],
            info = name)
        expect_identical(spectra$msPurityPurityStatus,
            pa@spectra$msPurityPurityStatus[ids], info = name)
        if (case$remove) {
            expect_identical(lengths(Spectra::mz(spectra)),
                c(0L, 0L, 1L, 0L, 0L)[ids], info = name)
            peaks <- Spectra::peaksData(spectra,
                columns = "original_peak_index")
            expect_identical(unname(peaks[[3L]][, 1L]), 2, info = name)
        } else {
            peaks <- Spectra::peaksData(spectra, columns = "pass_flag")
            expect_identical(lapply(peaks, function(m) unname(m[, 1L])),
                list(c(0, 0), c(0, 0), c(0, 1), c(0, 0), numeric())[ids],
                info = name)
        }
        expect_identical(out@spectra$msPurityStagePurity,
            pa@spectra$msPurityStagePurity, info = name)
    }
})

test_that("selected thresholds and skipped criteria retain actual history", {
    pa <- .msn_assessed_filter_fixture()
    params <- list(
        hcd_ms3 = FilterFragSpectraParam(
            select = list(msLevel = 3L, activationType = "HCD"),
            plim = 0.7, ilim = 10, snr = 0),
        ms4 = FilterFragSpectraParam(select = list(msLevel = 4L),
            plim = NA, ilim = 10, snr = 0)
    )
    implicit <- filterFragSpectra(pa, msLevel = 2:4, allfrag = TRUE,
        param = params, ilim = 10, snr = 0)
    explicit <- filterFragSpectra(pa, msLevel = 2:4, allfrag = TRUE,
        param = params, plim = 0.95, ilim = 10, snr = 0)
    expect_identical(implicit@fragSpectra$purity_pass_flag,
        c(FALSE, TRUE, TRUE, TRUE, TRUE))
    expect_identical(explicit@fragSpectra$purity_pass_flag,
        c(FALSE, TRUE, TRUE, FALSE, TRUE))
    expect_identical(implicit@fragSpectra$msPurityFilterThreshold,
        c(0.8, 0.7, NA, NA, NA))
    expect_identical(explicit@fragSpectra$msPurityFilterThreshold,
        c(0.95, 0.7, NA, 0.95, NA))
    record <- processingParameters(explicit, "filterFragSpectra")[[1L]]
    expect_identical(record$purity_filter,
        data.frame(pid = 1:5, threshold = c(0.95, 0.7, NA, 0.95, NA),
            decision = c("failed", "passed", "skipped", "missing",
                "skipped")))
    expect_identical(record$configuration$assignments$fallback_used,
        c(TRUE, FALSE, FALSE, TRUE, FALSE))
    effective <- record$configuration$parameter_sets[
        record$configuration$assignments$parameter_id]
    expect_identical(unname(vapply(effective, `[[`, numeric(1), "plim")),
        c(0.95, 0.7, NA, 0.95, NA))

    skipped <- filterFragSpectra(pa, msLevel = 2:4, allfrag = TRUE,
        plim = NA, ilim = 10, snr = 0)
    expect_true(all(skipped@fragSpectra$purity_pass_flag))
    expect_identical(skipped@fragSpectra$msPurityFilterDecision,
        rep("skipped", 5L))
    expect_identical(skipped@fragSpectra$msPurityPurityStatus,
        pa@spectra$msPurityPurityStatus)
})
