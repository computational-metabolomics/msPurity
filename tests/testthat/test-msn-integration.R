test_that("association level and parameter options follow existing arguments", {
    existing <- c("pa", "xcmsObj", "ppm", "plim", "intense", "convert2RawRT",
        "useGroup", "createDb", "outDir", "dbName", "grpPeaklist", "use_group",
        "out_dir", "create_db", "grp_peaklist", "db_name", "xset")
    generic <- methods::getGeneric("frag4feature")
    method <- methods::selectMethod("frag4feature", "purityA")
    expected <- c(existing, "msLevel", "parentReference", "param")
    expect_identical(names(formals(generic)), expected)
    expect_identical(names(formals(method)), expected)
    expect_identical(formals(method)$msLevel, 2L)
    expect_null(formals(method)$param)
    expect_identical(formals(method)$parentReference, "auto")
    expect_false("linkMSnSpectra" %in% getNamespaceExports("msPurity"))
    pa <- methods::new("purityA")
    expect_error(frag4feature(pa, NULL, includeMSn = TRUE), "unused argument")
    invalid_levels <- list(empty = integer(), missing = NA_integer_,
        fractional = 2.5, unsupported = c(2L, 5L), character = "4")
    for (name in names(invalid_levels))
        expect_error(frag4feature(pa, NULL, msLevel = invalid_levels[[name]]),
            "msLevel", info = name)
    expect_error(frag4feature(pa, NULL, parentReference = "invalid"),
        "arg.*should be one of")
    expect_error(frag4feature(pa, NULL, msLevel = 2:4, createDb = TRUE),
        "deprecated createDb route")
})

test_that("stored ancestry follows acquisition metadata rather than flags", {
    source <- .ancestry_spectra()
    metadata <- msPurity:::.pa_msn_spectra_metadata(source)
    preparation <- msPurity:::.pa_msn_ancestry(metadata, "sourceSpectrumId")
    preparation$source_metadata <- metadata
    pa <- methods::new("purityA")
    pa@spectra <- source
    pa@params$msnLineage <- preparation
    expect_true(msPurity:::.pa_msn_validate_preparation(pa))
    expect_identical(msPurity:::.pa_msn_preparation(pa), preparation)
    expect_error(msPurity:::.pa_msn_preparation(pa, "acquisitionNum"),
        "conflicts with the prepared parents")
    pa@spectra$purity_pass_flag <- c(TRUE, NA, NA)
    expect_true(msPurity:::.pa_msn_validate_preparation(pa))
    reordered <- pa
    reordered@spectra <- pa@spectra[c(3L, 1L, 2L)]
    expect_true(msPurity:::.pa_msn_validate_preparation(reordered))
    changed <- pa
    changed@spectra$precursorMz[2L] <- 151
    expect_error(msPurity:::.pa_msn_validate_preparation(changed),
        "no longer matches acquisition metadata")
    changed <- pa
    changed@spectra$precScanNum[3L] <- 101L
    expect_error(msPurity:::.pa_msn_validate_preparation(changed),
        "no longer matches acquisition metadata")
    changed <- pa
    changed@spectra <- changed@spectra[1:2]
    expect_error(msPurity:::.pa_msn_validate_preparation(changed),
        "no longer matches the source spectra")
})

.msn_association_fixture <- local({
    cache <- list()
    function(compound) {
        if (is.null(cache[[compound]])) {
            prepared <- .msn_standard_feature(compound, "MSn")
            stopifnot(length(prepared$selected) == 1L,
                length(prepared$feature) == 1L)
            prepared$pa <- purityA(.msn_standard_file(compound, "MSn"),
                cores = 1L, plotP = FALSE, msLevel = 2:4,
                parentReference = "acquisitionNum")
            cache[[compound]] <<- prepared
        }
        cache[[compound]]
    }
})

for (compound in c("tryptophan", "kynurenine")) local({
    name <- compound
    test_that(paste(name, "links descendants before MS2 purity exclusion"), {
        skip_if_not_installed("xcms")
        fixture <- .msn_association_fixture(name)
        purity_history <- processingHistory(fixture$pa, "purityA")
        local_mocked_bindings(.pa_msn_purity = function(...) {
            stop("Feature association must not calculate purity")
        }, .package = "msPurity")
        ordinary <- frag4feature(fixture$pa, fixture$xe, ppm = 10,
            convert2RawRT = FALSE)
        expect_true(all(Spectra::msLevel(ordinary@fragSpectra) == 2L))
        expect_identical(ordinary@params$msnLineage,
            fixture$pa@params$msnLineage)
        full <- frag4feature(fixture$pa, fixture$xe, ppm = 10,
            msLevel = 2:4, parentReference = "acquisitionNum",
            convert2RawRT = FALSE)
        selected <- full@grped_df[
            full@grped_df$grpid %in% fixture$feature, , drop = FALSE]
        expected <- if (name == "tryptophan") c(17L, 11L, 6L) else
            c(34L, 21L, 12L)
        expect_identical(tabulate(selected$msLevel - 1L, nbins = 3L), expected)
        source_rows <- match(selected$pid, fixture$pa@spectra$pid)
        expect_equal(selected$inPurity,
            fixture$pa@spectra$inPurity[source_rows])
        expected <- full@grped_df[full@grped_df$msLevel > 2L, , drop = FALSE]
        rejected <- frag4feature(fixture$pa, fixture$xe, ppm = 10, plim = 1,
            msLevel = 2:4, parentReference = "acquisitionNum",
            convert2RawRT = FALSE)
        actual <- rejected@grped_df
        rownames(actual) <- rownames(expected) <- NULL
        expect_identical(actual, expected)
        expect_true(all(Spectra::msLevel(rejected@fragSpectra) > 2L))
        history <- processingHistory(rejected, "frag4feature")
        parameters <- history[[length(history)]]$parameters
        expect_identical(parameters$msLevel, 2:4)
        expect_identical(parameters$resolvedParentReference,
            "acquisitionNum")
        relinked <- frag4feature(rejected, fixture$xe, ppm = 10,
            convert2RawRT = FALSE)
        expect_identical(relinked@params$msnLineage,
            fixture$pa@params$msnLineage)
        expect_true(all(Spectra::msLevel(relinked@fragSpectra) == 2L))
        expect_identical(processingHistory(relinked, "purityA"), purity_history)
        expect_identical(relinked@spectra$inPurity, fixture$pa@spectra$inPurity)
        expect_gt(length(processingHistory(relinked)),
            length(processingHistory(rejected)))
    })
})

test_that("group-width association retains source acquisition times", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    linked <- frag4feature(fixture$pa, fixture$xe, ppm = 10,
        msLevel = 2:4, parentReference = "acquisitionNum",
        convert2RawRT = FALSE, useGroup = TRUE)
    children <- linked@grped_df[linked@grped_df$msLevel > 2L, , drop = FALSE]
    source <- fixture$pa@spectra
    rows <- match(children$pid, source$pid)
    expect_gt(nrow(children), 0L)
    expect_false(anyNA(rows))
    expect_equal(children$retentionTime, Spectra::rtime(source)[rows])
    expect_equal(children$inPurity, source$inPurity[rows])
    expect_identical(linked@f4f_link_type, "group")
})

test_that("empty associations retain diagnostics without reading peaks", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    source <- .ancestry_spectra()[2:3]
    source$dataOrigin <- rep(unname(fixture$pa@fileList)[1L], 2L)
    source$fileid <- rep(1L, 2L)
    source$seqNum <- source$scanIndex
    pa <- methods::new("purityA", fileList = fixture$pa@fileList, cores = 1)
    pa@spectra <- source
    local_mocked_bindings(.pa_raw_scans = function(...) {
        stop("Empty associations must not request peak arrays")
    }, .package = "msPurity")
    higher_only <- frag4feature(pa, fixture$xe, msLevel = 2:4,
        parentReference = "sourceSpectrumId", convert2RawRT = FALSE)
    unmatched_param <- Frag4FeatureParam(ppm = 0,
        parentReference = "acquisitionNum", convert2RawRT = FALSE)
    unmatched <- frag4feature(fixture$pa, fixture$xe, msLevel = 2:4,
        param = unmatched_param)
    cases <- list(higher_only = higher_only, unmatched = unmatched)
    for (name in names(cases)) {
        out <- cases[[name]]
        expect_identical(nrow(out@grped_df), 0L, info = name)
        expect_true(all(c("pid", "grpid", "sample", "inPurity") %in%
                    names(out@grped_df)), info = name)
        expect_identical(length(out@fragSpectra), 0L, info = name)
        expect_identical(length(groupedSpectra(out)), 0L, info = name)
        record <- tail(processingHistory(out, "frag4feature"), 1L)[[1L]]
        expect_identical(record$status, "completed", info = name)
        expect_identical(nrow(record$output$spectra), 0L, info = name)
    }
    unused <- tail(processingHistory(higher_only, "frag4feature"), 1L)[[1L]]
    expect_equal(unused$parameters$effective_cores, 0)
    expect_identical(unused$input$xcms$association_input, "not_used")
    expect_null(unused$input$xcms$fingerprint)
    expect_identical(nrow(higher_only@params$msnLineage$scans), 2L)
    expect_true(all(is.na(higher_only@params$msnLineage$scans$root_pid)))
    expect_true(any(unmatched@params$msnLineage$scans$path_status == "valid"))
    settings <- tail(processingParameters(unmatched, "frag4feature"), 1L)[[1L]]
    expect_identical(settings$ppm, 0)
    expect_identical(settings$param, unmatched_param)
    expect_identical(settings$effective_msLevel, integer())
})

test_that("unavailable MS2 purity cannot create missing feature-link rows", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    pa <- fixture$pa
    pa@spectra$inPurity[Spectra::msLevel(pa@spectra) == 2L] <- NA_real_
    pa@puritydf$inPurity[] <- NA_real_
    ordinary <- frag4feature(pa, fixture$xe, ppm = 10, plim = 0.7,
        convert2RawRT = FALSE)
    expect_equal(nrow(ordinary@grped_df), 0L)
    expect_length(ordinary@fragSpectra, 0L)
    higher <- frag4feature(pa, fixture$xe, ppm = 10, plim = 0.7,
        msLevel = 2:4, parentReference = "acquisitionNum",
        convert2RawRT = FALSE)
    expect_gt(nrow(higher@grped_df), 0L)
    expect_false(anyNA(higher@grped_df$pid))
    expect_true(all(higher@grped_df$msLevel > 2L))
})

test_that("unassessed MS2 roots support association by recorded precursor", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    pa <- purityA(.msn_standard_file("tryptophan", "MSn"), msLevel = 4L,
        parentReference = "acquisitionNum")
    history <- processingHistory(pa, "purityA")
    expect_true(all(is.na(purityTable(pa)$aMz)))
    local_mocked_bindings(.pa_msn_purity = function(...) {
        stop("Feature association must not calculate missing purity")
    }, .package = "msPurity")

    # Without an assessed MS2 target, use its recorded precursor mass.
    linked <- frag4feature(pa, fixture$xe, ppm = 10, intense = FALSE,
        msLevel = 2:4, convert2RawRT = FALSE)
    children <- linked@grped_df[linked@grped_df$msLevel == 4L, , drop = FALSE]
    rows <- match(children$pid, pa@spectra$pid)

    expect_gt(nrow(children), 0L)
    expect_false(anyNA(rows))
    expect_equal(children$inPurity, pa@spectra$inPurity[rows])
    expect_identical(processingHistory(linked, "purityA"), history)
})

test_that("typed association agrees with scalars and retains typed history", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    param <- Frag4FeatureParam(ppm = 10, plim = 0.5, intense = FALSE,
        convert2RawRT = FALSE, useGroup = TRUE)
    scalar <- frag4feature(fixture$pa, fixture$xe, ppm = 10, plim = 0.5,
        intense = FALSE, convert2RawRT = FALSE, useGroup = TRUE, msLevel = 2:4)
    typed <- frag4feature(fixture$pa, fixture$xe, msLevel = 2:4, param = param)
    expect_message(mixed <- frag4feature(fixture$pa, fixture$xe,
            ppm = 10L, plim = 0.5, intense = FALSE, convert2RawRT = FALSE,
            useGroup = TRUE, use_group = TRUE, parentReference = "auto",
            msLevel = 2:4, param = param), "use_group is deprecated")

    expect_gt(nrow(typed@grped_df), 0L)
    expect_identical(typed@grped_df, scalar@grped_df)
    expect_identical(mixed@grped_df, scalar@grped_df)
    expect_identical(Spectra::peaksData(typed@fragSpectra),
        Spectra::peaksData(scalar@fragSpectra))
    record <- processingParameters(typed, "frag4feature")[[1L]]
    expect_identical(record$param, param)
    expect_identical(record[names(as.list(param))], as.list(param))
    expect_identical(record$msLevel, 2:4)
    expect_identical(record$effective_msLevel, 2:4)
    expect_identical(processingHistory(typed, "purityA"),
        processingHistory(fixture$pa, "purityA"))

    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(typed, path)
    expect_identical(processingHistory(readRDS(path)), processingHistory(typed))
})

test_that("MS4-only links retain repeated occurrences and source parents", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    xe <- fixture$xe
    definitions <- xcms::featureDefinitions(xe)
    repeated <- definitions[rep(fixture$feature, 3L), , drop = FALSE]
    rownames(repeated) <- make.unique(rownames(repeated))
    xcms::featureDefinitions(xe) <- repeated
    local_mocked_bindings(.pa_msn_ancestry = function(...) {
        stop("Automatic association must reuse prepared parents")
    }, .pa_msn_purity = function(...) {
        stop("Feature association must not calculate purity")
    }, .package = "msPurity")
    param <- Frag4FeatureParam(ppm = 10, plim = 1, convert2RawRT = FALSE)
    out <- frag4feature(fixture$pa, xe, msLevel = 4L, param = param)

    expect_identical(nrow(out@grped_df), 18L)
    expected_groups <- factor(rep(as.character(1:3), 6L),
        levels = as.character(1:3))
    expect_identical(out@grped_df$grpid, expected_groups)
    expect_identical(out@grped_df$pid, rep(out@fragSpectra$pid, each = 3L))
    expect_identical(Spectra::msLevel(out@fragSpectra), rep(4L, 6L))
    expect_identical(out@spectra$pid, fixture$pa@spectra$pid)
    expect_identical(Spectra::peaksData(out@spectra),
        Spectra::peaksData(fixture$pa@spectra))
    expect_identical(out@params$msnLineage, fixture$pa@params$msnLineage)
    record <- processingParameters(out, "frag4feature")[[1L]]
    expect_identical(record$msLevel, 4L)
    expect_identical(record$effective_msLevel, 4L)
    expect_identical(record$resolvedParentReference, "acquisitionNum")
})
