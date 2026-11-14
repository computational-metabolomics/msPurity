context("Complete fragment-filter parameters and metadata selection")

.filter_param_defaults <- function() {
    list(ilim = 0, plim = 0.8, ra = 0, snr = 3, snmeth = "median",
        rmp = FALSE, allfrag = FALSE, cores = 1)
}

.filter_param_fixture <- function() {
    tab <- data.frame(pid = 1:4, fileid = 1L, seqNum = 1:4,
        acquisitionNum = 201:204, retentionTime = 10:13,
        precursorMZ = 205, inPurity = c(0.95, 0.75, NA, 0.95))
    sv <- data.frame(pid = 1:4, fileid = 1L, seqNum = 1:4,
        msLevel = 2L, scanIndex = 1:4, acquisitionNum = 201:204,
        rtime = as.numeric(10:13), precursorMz = 205,
        massAnalyzer = "Orbitrap",
        activationType = c("CID", "HCD", NA_character_, "ETD"),
        collisionEnergy = c(20, 30, NA, 20),
        collisionEnergyUnit = c("eV", "NCE", NA_character_, "eV"),
        polarity = c(1L, 1L, 0L, 1L), stringsAsFactors = FALSE)
    peaks <- rep(list(cbind(mz = c(100, 101, 102),
                intensity = c(1, 4, 10), original_peak_index = c(2, 4, 9))), 4)
    sp <- msPurity:::.msp_memory(sv, peaks,
        c("mz", "intensity", "original_peak_index"))
    sp <- msPurity:::.pa_set_purity_vars(sp, tab)
    msPurity:::.pa_spectra_metadata(sp)$msPurityRawPeaksAvailable <- TRUE
    pa <- methods::new("purityA")
    pa@fileList <- c(synthetic = "synthetic")
    pa@spectra <- sp
    pa@fragSpectra <- sp
    pa@puritydf <- tab
    pa@grped_df <- data.frame(grpid = rep(1L, 4), pid = 1:4,
        sample = 1L, inPurity = tab$inPurity)
    pa
}

.filter_param_peaks <- function(pa) {
    Spectra::peaksData(groupedSpectra(pa, legacy = FALSE),
        columns = c("mz", "intensity", "snr", "ra",
            "intensity_pass_flag", "ra_pass_flag", "snr_pass_flag",
            "pass_flag"))
}

.filter_param_options <- function(value, code) {
    old <- options(msPurity.legacySlots = value)
    on.exit(options(old))
    force(code)
}

test_that("filter parameters extend Param with complete fixed defaults", {
    param <- FilterFragSpectraParam()
    expect_s4_class(param, "Param")
    expect_true(methods::validObject(param))
    expect_identical(as.list(param), list(select = list(), ilim = 0,
            plim = 0.8, ra = 0, snr = 3, snmeth = "median"))
    expect_identical(as.list(FilterFragSpectraParam(plim = NA))$plim, NA_real_)
    expect_identical(as.list(FilterFragSpectraParam(snr = 0))$ilim, 0)
    expect_error(FilterFragSpectraParam(combine = "activationType"), "unused")
})

test_that("invalid values and selectors are rejected at construction", {
    expect_error(FilterFragSpectraParam(snr = -1), "nonnegative")
    expect_error(FilterFragSpectraParam(ilim = Inf), "finite")
    expect_error(FilterFragSpectraParam(plim = 2), "between zero and one")
    expect_error(FilterFragSpectraParam(plim = NaN), "nonnegative")
    expect_error(FilterFragSpectraParam(snmeth = "other"), "median")
    expect_error(FilterFragSpectraParam(select = list(other = 2)),
        "unsupported")
    expect_error(FilterFragSpectraParam(select = list(msLevel = "2")), "type")
    expect_error(FilterFragSpectraParam(select = list(msLevel = 1.5)),
        "integer")
    expect_error(FilterFragSpectraParam(select = list(polarity = 2)),
        "polarity")
    expect_error(FilterFragSpectraParam(select = list(activationType = "")),
        "unknown")
    expect_error(FilterFragSpectraParam(select = list(collisionEnergy = NaN)),
        "finite")
    expect_true(methods::validObject(FilterFragSpectraParam(
        select = list(activationType = NA, massAnalyzer = "Orbitrap"))))
})

test_that("parameter inputs require complete named typed objects", {
    pa <- .filter_param_fixture()
    expect_error(filterFragSpectra(pa, param = list(snr = 1)), "requires")
    expect_error(filterFragSpectra(pa, param = list(default = list(snr = 1))),
        "requires")
    expect_error(filterFragSpectra(pa, param = list(FilterFragSpectraParam())),
        "names")
    expect_error(filterFragSpectra(pa, param = list()), "names")
    invalid <- FilterFragSpectraParam()
    invalid@snr <- -1
    expect_error(filterFragSpectra(pa, param = invalid), "nonnegative")
    expect_error(filterFragSpectra(pa,
            param = FilterFragSpectraParam(), cores = 2),
        "cores = 1")
    expect_error(filterFragSpectra(pa, param = FilterFragSpectraParam(),
            allfrag = NA), "allfrag")
})

test_that("selectors use AND between fields OR within fields and explicit NA", {
    sp <- .filter_param_fixture()@fragSpectra
    param <- list(
        collision = FilterFragSpectraParam(select = list(
            activationType = c("CID", "HCD"), polarity = 1L), snr = 0),
        unknown = FilterFragSpectraParam(select = list(activationType = NA),
            snr = 2),
        unused = FilterFragSpectraParam(select = list(msLevel = 3L), snr = 1))
    result <- msPurity:::.pa_filter_resolve(sp, param, .filter_param_defaults())
    expect_identical(vapply(result$parameters, `[[`, numeric(1), "snr"),
        c(0, 0, 2, 3))
    expect_identical(result$assignments$parameter_set_names,
        list("collision", "collision", "unknown", character()))
    expect_identical(result$assignments$fallback_used,
        c(FALSE, FALSE, FALSE, TRUE))
    expect_identical(result$unused_parameter_sets, "unused")
})

test_that("metadata aliases do not silently conflict", {
    source <- .filter_param_fixture()@fragSpectra
    source$activationMethod <- source$activationType
    expect_identical(msPurity:::.pa_filter_metadata(source)$activationType,
        source$activationType)
    source$activationMethod[1] <- "HCD"
    expect_error(msPurity:::.pa_filter_metadata(source), "Conflicting")
    source$activationMethod <- NULL
    source$massAnalyzer <- NULL
    expect_true(all(is.na(msPurity:::.pa_filter_metadata(source)$massAnalyzer)))
    table <- as.data.frame(Spectra::spectraData(source,
            columns = setdiff(Spectra::spectraVariables(source),
                Spectra::peaksVariables(source))))
    expect_identical(msPurity:::.pa_filter_metadata(table),
        msPurity:::.pa_filter_metadata(source))
})

test_that("equivalent overlaps share effective settings and conflicts fail", {
    pa <- .filter_param_fixture()
    param <- list(all = FilterFragSpectraParam(snr = 0),
        cid = FilterFragSpectraParam(select = list(activationType = "CID"),
            snr = 0))
    resolved <- msPurity:::.pa_filter_resolve(pa@fragSpectra, param,
        .filter_param_defaults())
    expect_identical(resolved$assignments$parameter_set_names[[1]],
        c("all", "cid"))
    expect_length(resolved$parameter_sets, 1)
    history <- processingHistory(pa)
    param$cid <- FilterFragSpectraParam(select = list(activationType = "CID"),
        snr = 2)
    expect_error(filterFragSpectra(pa, param = param),
        "Conflicting overlapping")
    expect_identical(processingHistory(pa), history)
    expect_identical(pa@fragSpectra$activationType, c("CID", "HCD", NA, "ETD"))
})

test_that("typed fixed settings equal corresponding scalar scientific output", {
    pa <- .filter_param_fixture()
    # Keep missing purity out of the scalar comparison; its historical NA
    # behavior remains distinct from the parameter route's explicit policy.
    pa@spectra$inPurity[3] <- 0.95
    pa@grped_df$inPurity[3] <- 0.95
    typed <- filterFragSpectra(pa, snr = 100,
        param = FilterFragSpectraParam(plim = 0.7, snr = 0))
    scalar <- filterFragSpectra(pa, plim = 0.7, snr = 0)
    expect_identical(.filter_param_peaks(typed), .filter_param_peaks(scalar))
    expect_identical(typed@grped_df, scalar@grped_df)
    expect_identical(Spectra::msLevel(typed@fragSpectra), rep(2L, 4))
    expect_identical(typed@fragSpectra$activationType,
        c("CID", "HCD", NA_character_, "ETD"))
})

test_that("unmatched spectra use scalar fallback without being discarded", {
    out <- filterFragSpectra(.filter_param_fixture(), plim = 0.7, snr = 3,
        param = FilterFragSpectraParam(select = list(activationType = "CID"),
            snr = 0))
    expect_length(out@fragSpectra, 4)
    expect_identical(out@fragSpectra$msPurityParameterFallback,
        c(FALSE, TRUE, TRUE, TRUE))
    pk <- .filter_param_peaks(out)
    expect_identical(unname(pk[[1]][, "pass_flag"]), rep(1, 3))
    expect_identical(unname(pk[[2]][, "pass_flag"]), rep(0, 3))
    expect_identical(out@fragSpectra$purity_pass_flag,
        c(TRUE, TRUE, FALSE, TRUE))
})

test_that("typed purity skipping is explicit and missing purity fails", {
    param <- FilterFragSpectraParam(snr = 0, plim = NA)
    skipped <- filterFragSpectra(.filter_param_fixture(), param = param)
    assessed <- filterFragSpectra(.filter_param_fixture(),
        param = FilterFragSpectraParam(snr = 0))
    expect_true(all(skipped@fragSpectra$purity_pass_flag))
    expect_identical(skipped@fragSpectra$msPurityPurityStatus,
        rep("skipped", 4))
    expect_false(assessed@fragSpectra$purity_pass_flag[3])
    expect_identical(assessed@fragSpectra$msPurityPurityStatus[3], "missing")
})

test_that("destructive filtering retains peak positions and refuses rerun", {
    pa <- .filter_param_fixture()
    out <- filterFragSpectra(pa, rmp = TRUE,
        param = FilterFragSpectraParam(ilim = 4, plim = NA, snr = 0))
    pk <- Spectra::peaksData(out@fragSpectra,
        columns = c("mz", "intensity", "original_peak_index"))
    expect_true(all(vapply(pk, function(m) identical(unname(m[, "mz"]), 102),
                logical(1))))
    expect_true(all(vapply(pk, function(m) {
        identical(unname(m[, "original_peak_index"]), 9)
    }, logical(1))))
    before <- processingHistory(out)
    expect_message(refused <- filterFragSpectra(out,
            param = FilterFragSpectraParam()), "previously filtered")
    expect_identical(processingHistory(refused), before)
    expect_identical(Spectra::peaksData(refused@fragSpectra),
        Spectra::peaksData(out@fragSpectra))
})

test_that("queued peak filtering retains positions and empty spectra", {
    pa <- .filter_param_fixture()
    pa@fragSpectra <- Spectra::filterMzRange(pa@fragSpectra,
        mz = c(100.5, 102.5))
    out <- filterFragSpectra(pa,
        param = FilterFragSpectraParam(plim = NA, snr = 0))
    pk <- Spectra::peaksData(out@fragSpectra,
        columns = c("mz", "original_peak_index"))
    expect_identical(unname(pk[[1]][, "mz"]), c(101, 102))
    expect_identical(unname(pk[[1]][, "original_peak_index"]), c(4, 9))
    pa@fragSpectra <- Spectra::filterMzRange(pa@fragSpectra, mz = c(200, 201))
    empty <- filterFragSpectra(pa,
        param = FilterFragSpectraParam(plim = NA, snr = 0))
    expect_identical(lengths(Spectra::mz(empty@fragSpectra)),
        rep(0L, 4))
})

test_that("duplicate peak identities require tracking before queued changes", {
    data <- S4Vectors::DataFrame(msLevel = 2L)
    data$mz <- list(c(100, 100, 150))
    data$intensity <- list(c(20, 20, 50))
    source <- Spectra::Spectra(data, source = Spectra::MsBackendMemory())
    expect_false(msPurity:::.pa_spectra_has_queue(source))
    expect_identical(msPurity:::.pa_filter_source_peaks(source)$
            original_positions[[1L]], 1:3)

    queued <- Spectra::addProcessing(source,
        function(x, ...) x[c(2L, 1L, 3L), , drop = FALSE])
    expect_true(msPurity:::.pa_spectra_has_queue(queued))
    expect_error(msPurity:::.pa_filter_source_peaks(queued),
        "ambiguous for duplicate peaks")

    tracked <- msPurity:::.msp_memory(as.data.frame(data[, "msLevel",
                drop = FALSE]),
        list(cbind(mz = c(100, 100, 150), intensity = c(20, 20, 50),
                original_peak_index = 1:3)),
        c("mz", "intensity", "original_peak_index"))
    tracked <- Spectra::addProcessing(tracked,
        function(x, ...) x[c(2L, 1L, 3L), , drop = FALSE])
    expect_identical(msPurity:::.pa_filter_source_peaks(tracked)$
            original_positions[[1L]], c(2, 1, 3))
    retained <- Spectra::filterMzRange(tracked, mz = c(90, 110))
    expect_identical(msPurity:::.pa_filter_source_peaks(retained)$
            original_positions[[1L]], c(2, 1))
    expect_error(msPurity:::.pa_spectra_has_queue(
        methods::new("purityA")), "Unsupported Spectra")
})

test_that("allfrag includes unlinked spectra with their own metadata", {
    pa <- .filter_param_fixture()
    pa@grped_df <- pa@grped_df[1:3, , drop = FALSE]
    out <- filterFragSpectra(pa, allfrag = TRUE,
        param = FilterFragSpectraParam(plim = NA, snr = 0))
    expect_length(allFragSpectra(out, legacy = FALSE), 4)
    expect_identical(out@fragSpectra$acquisitionNum, 201:204)
    expect_identical(out@fragSpectra$collisionEnergyUnit,
        c("eV", "NCE", NA_character_, "eV"))
    expect_true(all(out@grped_df$purity_pass_flag))
    expect_length(processingHistory(out, "filterFragSpectra"), 1)
})

test_that("queued selection cannot hide invalid stored peak identities", {
    source <- msPurity:::.msp_memory(data.frame(msLevel = 2L),
        list(cbind(mz = c(100, 150), intensity = c(20, 50),
                original_peak_index = c(1, 1))),
        c("mz", "intensity", "original_peak_index"))
    selected <- Spectra::filterMzRange(source, mz = c(149, 151))
    expect_error(msPurity:::.pa_filter_source_peaks(selected),
        "Invalid original peak positions before processing")
})

test_that("history retains typed settings and actual serial workers", {
    param <- FilterFragSpectraParam(plim = NA, snr = 0)
    out <- filterFragSpectra(.filter_param_fixture(), param = param)
    record <- processingParameters(out, "filterFragSpectra")[[1]]
    expect_identical(record$configuration$parameter_objects$default, param)
    expect_identical(record$cores, 1)
    expect_identical(record$effective_cores, 1L)
    expect_false(any(c("rule_ids", "parameter_rule_ids", "unused_rules") %in%
                names(record$configuration)))
})

test_that("legacy slot setting does not change parameter output", {
    param <- FilterFragSpectraParam(plim = NA, snr = 0)
    on <- .filter_param_options(TRUE,
        filterFragSpectra(.filter_param_fixture(), param = param))
    off <- .filter_param_options(FALSE,
        filterFragSpectra(.filter_param_fixture(), param = param))
    expect_identical(.filter_param_peaks(on), .filter_param_peaks(off))
    expect_length(off@grped_ms2, 0)
    expect_identical(on@fragSpectra$msPurityParameterSetNames,
        off@fragSpectra$msPurityParameterSetNames)
})

test_that("empty feature selection retains a processing record", {
    pa <- .filter_param_fixture()
    pa@grped_df <- pa@grped_df[FALSE, , drop = FALSE]
    out <- filterFragSpectra(pa, param = FilterFragSpectraParam())
    expect_length(out@fragSpectra, 0)
    expect_length(processingHistory(out, "filterFragSpectra"), 1)
})

test_that("repeated feature links do not duplicate filtering assignments", {
    pa <- .filter_param_fixture()
    extra <- pa@grped_df[1, , drop = FALSE]
    extra$grpid <- 2L
    pa@grped_df <- rbind(pa@grped_df, extra)
    out <- filterFragSpectra(pa,
        param = FilterFragSpectraParam(plim = NA, snr = 0))
    record <- processingParameters(out, "filterFragSpectra")[[1]]
    expect_length(out@fragSpectra, 4)
    expect_equal(nrow(out@grped_df), 5)
    expect_equal(nrow(record$configuration$assignments), 4)
    expect_true(all(out@grped_df$purity_pass_flag))
})

test_that("scalar generic positions and actual serial workers are retained", {
    pa <- .filter_param_fixture()
    positional <- filterFragSpectra(pa, 0, 0.8, 0, 0, 2, FALSE,
        "median", FALSE)
    named <- filterFragSpectra(pa, ilim = 0, plim = 0.8, ra = 0,
        snr = 0, cores = 2, rmp = FALSE, snmeth = "median", allfrag = FALSE)
    expect_identical(.filter_param_peaks(positional),
        .filter_param_peaks(named))
    record <- processingParameters(positional, "filterFragSpectra")[[1]]
    expect_identical(record$cores, 2)
    expect_identical(record$effective_cores, 1L)
    expect_false("cores" %in% names(positional@filter_frag_params))
    expect_identical(names(formals(methods::getMethod(
        "filterFragSpectra", "purityA"))),
    names(formals(methods::getGeneric("filterFragSpectra"))))
})


test_that("default empty objects support both typed filtering selections", {
    for (allfrag in c(FALSE, TRUE)) {
        out <- filterFragSpectra(methods::new("purityA"),
            param = FilterFragSpectraParam(), allfrag = allfrag)
        expect_length(out@spectra, 0)
        expect_length(out@fragSpectra, 0)
        records <- processingHistory(out, "filterFragSpectra")
        expect_length(records, 1)
        expect_identical(records[[1]]$status, "completed")
        expect_equal(nrow(records[[1]]$output$spectra), 0)
        expect_equal(nrow(records[[1]]$parameters$configuration$assignments), 0)
    }
})
