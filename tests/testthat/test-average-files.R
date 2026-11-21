# Two files, one feature. File 1: MS2 pid 1 (precursor 300; m/z 100, 150)
# and pid 2 (300.0005; m/z 100, 150), both HCD 30; MS3 pid 5 from pid 1
# (precursor 150) and pid 6 from pid 2 (150.0006), both CID 35. File 2: MS2
# pid 3 (m/z 100, 150) and pid 4 (m/z 200), HCD 30; MS3 pid 7 (150.0004) and
# pid 8 (precursor 100) from pid 3, CID 35.
.average_files_fixture <- function(
    ms3_precursors = c(150, 150.0006, 150.0004, 100)) {
    data <- S4Vectors::DataFrame(
        msLevel = rep(c(1L, 2L, 2L, 3L, 3L), 2L),
        acquisitionNum = rep(1:5, 2L), scanIndex = rep(1:5, 2L),
        rtime = rep(10 * (1:5), 2L),
        precScanNum = c(NA, 1L, 1L, 2L, 3L, NA, 1L, 1L, 2L, 2L),
        precursorMz = c(NA, 300, 300.0005, ms3_precursors[1:2],
            NA, 300.0008, 300.0008, ms3_precursors[3:4]),
        activationType = rep(c(NA, "HCD", "HCD", "CID", "CID"), 2L),
        collisionEnergy = rep(c(NA, 30, 30, 35, 35), 2L),
        polarity = 1L, massAnalyzer = "Orbitrap",
        dataOrigin = rep(c("file1", "file2"), each = 5L),
        spectrumId = rep(paste0("scan=", 1:5), 2L))
    data$mz <- list(500, c(100, 150), c(100, 150), c(60, 80), c(60, 80),
        500, c(100, 150), 200, c(60, 80), 50)
    data$intensity <- list(1000, c(1000, 500), c(800, 400), c(50, 100),
        c(40, 120), 1000, c(600, 300), 90, c(30, 60), 10)
    raw <- Spectra::Spectra(data, source = Spectra::MsBackendMemory())
    purity <- data.frame(pid = 1:4, fileid = c(1L, 1L, 2L, 2L),
        seqNum = c(2L, 3L, 2L, 3L), acquisitionNum = c(2L, 3L, 2L, 3L),
        retentionTime = c(20, 30, 20, 30),
        precursorMZ = c(300, 300.0005, 300.0008, 300.0008),
        precursorScanNum = 1L, inPurity = 0.9)
    pa <- methods::new("purityA", fileList = c(file1 = "file1",
        file2 = "file2"), cores = 1)
    pa@spectra <- msPurity:::.pa_msn_retained_spectra(raw,
        rep(1:2, each = 5L), purity)
    metadata <- msPurity:::.pa_msn_spectra_metadata(pa@spectra)
    lineage <- msPurity:::.pa_msn_ancestry(metadata, "acquisitionNum")
    lineage$source_metadata <- metadata
    pa@params$msnLineage <- lineage
    roots <- data.frame(grpid = 1L, pid = 1:4, sample = c(1L, 1L, 2L, 2L),
        cid = 1:4, inPurity = 0.9, precurMtchMZ = purity$precursorMZ,
        retentionTime = purity$retentionTime)
    pa@grped_df <- msPurity:::.pa_msn_feature_links(roots, lineage)
    peaks <- as.list(Spectra::peaksData(pa@spectra))
    linked <- split(peaks[match(pa@grped_df$pid, pa@spectra$pid)],
        as.character(pa@grped_df$grpid))
    pa@fragSpectra <- msPurity:::.pa_frag_from_grouped(pa, linked)
    msPurity:::.pa_sync_legacy(pa, puritydf = purity, grped_ms2 = linked)
}

# "<file or -> <sorted selected pids>" for each average.
.average_members <- function(averages) {
    pids <- vapply(.average_selected_links(averages), function(x)
        paste(sort(x$pid), collapse = ","), character(1))
    stats::setNames(seq_along(averages), paste(ifelse(is.na(averages$sample),
        "-", averages$sample), pids))
}

.average_peak_matrix <- function(average) {
    peaks <- Spectra::peaksData(average,
        columns = Spectra::peaksVariables(average))[[1L]]
    unname(as.matrix(peaks[order(peaks[, "mz"]), c("mz", "intensity",
        "count", "total"), drop = FALSE]))
}

.average_param <- function(...)
    AverageFragSpectraParam(minfrac = 0, sumi = FALSE, ...)

test_that("configured averaging honours cores without changing peaks or provenance", {
    pa <- updateObject(.average_files_fixture())
    pa@cores <- 2
    param <- .average_param(select = list(msLevel = 2L))
    intra <- averageIntraFragSpectra(pa, param = param, cores = 1)
    inputs <- list(All = pa, Intra = pa, Inter = intra)
    make_cluster <- parallel::makeCluster
    for (level in names(inputs)) {
        step <- paste0("average", level, "FragSpectra")
        args <- list(inputs[[level]], param = param, blockSize = 3L)
        serial <- do.call(step, c(args, list(cores = 1)))
        workers <- integer()
        parallel <- with_mocked_bindings(do.call(step, args),
            makeCluster = function(spec, ...) {
                workers <<- c(workers, spec)
                make_cluster(spec, ...)
            }, .package = "parallel")
        expect_true(length(workers) > 0L && all(workers == 2L), info = step)
        a <- averagedSpectra(serial, legacy = FALSE)
        b <- averagedSpectra(parallel, legacy = FALSE)
        expect_identical(Spectra::peaksData(a, columns = Spectra::peaksVariables(a)),
            Spectra::peaksData(b, columns = Spectra::peaksVariables(b)), info = step)
        expect_identical(a$spectrumId, b$spectrumId, info = step)
        for (id in a$spectrumId) {
            expect_identical(averageMembers(a, id), averageMembers(b, id), info = step)
            expect_identical(averagePeakMembers(a, id, depth = "original"),
                averagePeakMembers(b, id, depth = "original"), info = step)
        }
        recorded <- tail(processingParameters(parallel, step), 1L)[[1L]]
        expect_identical(recorded[c("cores", "effective_cores", "blockSize")],
            list(cores = 2, effective_cores = 2L, blockSize = 3L), info = step)
        expect_identical(tail(processingParameters(serial, step), 1L)[[1L]]$effective_cores,
            1L, info = step)
        expect_identical(serial@cores, 2, info = step)
        expect_identical(parallel@cores, 2, info = step)
        expect_error(do.call(step, c(args, list(cores = 0))),
            "cores.*positive integer", info = step)
    }
})

test_that("configured execution records direct and empty work without starting workers", {
    pa <- .average_files_fixture()
    param <- .average_param(select = list(msLevel = 2L))
    direct <- with_mocked_bindings(averageAllFragSpectra(pa, param = param, cores = 4),
        makeCluster = function(...) stop("Unexpected workers"), .package = "parallel")
    settings <- tail(processingParameters(direct, "averageAllFragSpectra"), 1L)[[1L]]
    expect_identical(settings[c("cores", "effective_cores")],
        list(cores = 4, effective_cores = 1L))
    empty <- averageAllFragSpectra(pa, param = .average_param(select = list(msLevel = 4L)),
        cores = 2)
    expect_identical(tail(processingParameters(empty,
        "averageAllFragSpectra"), 1L)[[1L]]$effective_cores, 0L)
})

test_that("intra-file averages keep files apart", {
    pa <- .average_files_fixture()
    pa@fragSpectra$dataOrigin <- paste0("/experiment/", pa@fragSpectra$dataOrigin)
    out <- averageIntraFragSpectra(pa, param = .average_param())
    averages <- averagedSpectra(out, legacy = FALSE)
    members <- .average_members(averages)

    # MS3 pid 5 and 6 differ by 4 ppm (within 5 ppm); pid 7 and 8 have other
    # precursors or another file.
    expect_true(all(averages$av_level == "av_intra"))
    expect_setequal(names(members), c("1 1,2", "1 5,6", "2 3,4", "2 7",
        "2 8"))
    expect_type(averages$spectrumId, "character")
    expect_false("msPurityKey" %in% Spectra::spectraVariables(averages))
    expect_equal(anyDuplicated(averages$spectrumId), 0L)
    # File 1 MS2: medians of 1000/800 and 500/400; file 2 MS2: pid 3 and 4
    # have no common peak.
    expect_equal(.average_peak_matrix(averages[members[["1 1,2"]]]),
        cbind(c(100, 150), c(900, 450), 2, 2))
    expect_equal(.average_peak_matrix(averages[members[["2 3,4"]]]),
        cbind(c(100, 150, 200), c(600, 300, 90), 1, 2))
    expect_equal(.average_peak_matrix(averages[members[["1 5,6"]]]),
        cbind(c(60, 80), c(45, 110), 2, 2))
    first <- averages[members[["1 1,2"]]]
    second <- averages[members[["2 3,4"]]]
    # Only pid 4 supplies the queued peak. The stored average used both scans.
    queued <- Spectra::filterMzRange(second, mz = c(190, 210))
    result <- spectralMatching(first, queued)
    expect_identical(result$matchedResults$sharedMeasurement, FALSE)
    expect_identical(result$matchedResults$sharedSourceFile, FALSE)
    expect_equal(result$matchedResults$mcount, 0)
    sources <- subset(result$sourceMeasurements, side == "reference")
    expect_setequal(sources$scanNumber, c("2", "3"))
})

test_that("native spectrum IDs retain their measured file and feature-link identities", {
    pa <- .average_files_fixture()
    rows <- pa@grped_df[pa@grped_df$msLevel == 2L, , drop = FALSE]
    # A second chromatographic peak of feature 1 also matches scan 1.
    extra <- rows[1L, , drop = FALSE]
    extra$cid <- 22L
    pa@grped_df <- rbind(rows, extra)
    # Both files use native "9", also the first derived ID in this fixture.
    for (slot in c("spectra", "fragSpectra")) {
        sp <- methods::slot(pa, slot)
        sp$spectrumId[sp$pid %in% c(1L, 3L)] <- "9"
        sp$spectrumId[sp$pid == 2L] <- NA_character_
        methods::slot(pa, slot) <- sp
    }
    # Supplied measurement IDs need not be consecutive integer keys.
    pa@spectra$pid <- 10.5 + pa@spectra$pid * 2
    pa@fragSpectra$pid <- 10.5 + pa@fragSpectra$pid * 2
    pa@grped_df$pid <- 10.5 + pa@grped_df$pid * 2

    out <- averageIntraFragSpectra(pa, param = .average_param())
    intra <- averagedSpectra(out, legacy = FALSE)
    id <- intra$spectrumId[intra$sample == "1"]
    expect_identical(id, "9")
    selected <- averageMembers(out, spectrumId = id)
    expect_identical(unique(selected$merged), id)
    expect_identical(selected$member[selected$pid == 14.5], NA_character_)
    repeated <- selected[selected$pid == 12.5, ]
    expect_setequal(repeated$cid, c(1L, 22L))
    expect_identical(unique(repeated$member), "9")
    expect_equal(length(unique(repeated$linkId)), 2L)
    edges <- averagePeakMembers(out, id, peak = 1L)
    expect_setequal(edges$linkId[edges$pid == 12.5], repeated$linkId)
    both <- averageProvenance(out, intra$spectrumId)$mergeMembers
    expect_setequal(both$fileid[both$member %in% "9"], 1:2)
    expect_setequal(both$pid[both$member %in% "9"], c(12.5, 16.5))
    expect_error(averageMembers(out, paste0("0", id)), "spectrumId")
    expect_error(averageMembers(out, NA_character_), "spectrumId")
    expect_identical(out@spectra$spectrumId, pa@spectra$spectrumId)
})

test_that("averages retain acquisition files after reordered or subset linking", {
    skip_if_no_mspuritydata()
    original <- msPurity:::.pa_update(.golden_pa("1_purityA_pa.rds"))
    xcms <- .golden_xcms("msms_only_xset.rds")
    cases <- list(reordered = c(2L, 1L), subset = 2L)
    param <- .average_param(combine = c("activationType", "collisionEnergy"))

    for (name in names(cases)) {
        files <- cases[[name]]
        pa <- original
        pa@fileList <- original@fileList[files]
        sp <- original@spectra
        sp <- sp[sp$fileid %in% files]
        sp$fileid <- match(sp$fileid, files)
        pa@spectra <- sp
        linked <- frag4feature(pa, xcms)

        # XCMS stays in its original order: acquisition file 1 is sample 2.
        expect_true(nrow(linked@grped_df) > 0L, info = name)
        expect_equal(linked@grped_df$sample,
            files[linked@grped_df$fileid], info = name)
        out <- averageIntraFragSpectra(linked, param = param)
        out <- averageInterFragSpectra(out, param = param)
        averages <- averagedSpectra(out, legacy = FALSE)
        intra <- averages[averages$av_level == "av_intra"]
        expect_true(setequal(intra$sample, as.character(seq_along(files))),
            info = name)
        for (i in seq_along(intra)) {
            contributors <- intra$msPurityContributors[[i]]
            file <- unique(contributors$fileid)
            case <- paste(name, intra$spectrumId[i])
            expect_equal(length(file), 1L, info = case)
            expect_identical(intra$sample[i], as.character(file), info = case)
            inputs <- averageMembers(intra, intra$spectrumId[i])
            expect_true(all(contributors$pid %in% inputs$pid), info = case)
        }
        for (i in seq_along(averages)) {
            contributors <- averages$msPurityContributors[[i]]
            expect_identical(contributors$sample,
                as.character(files[contributors$fileid]),
                info = paste(name, averages$spectrumId[i]))
        }
    }
})

test_that("inter-file averages combine compatible intra-file averages", {
    pa <- averageIntraFragSpectra(.average_files_fixture(),
        param = .average_param())
    out <- averageInterFragSpectra(pa, param = .average_param())
    averages <- averagedSpectra(out, legacy = FALSE)
    intra <- averages[averages$av_level == "av_intra"]
    inter <- averages[averages$av_level == "av_inter"]
    members <- .average_members(inter)

    expect_false(any(c("msPuritySelectedContributors", "msPurityAverageInputs")
        %in% Spectra::spectraVariables(averages)))
    expect_length(intra, 5L)
    # MS3: {5, 6} (150, 150.0006) and {7} (150.0004) agree pairwise; {8} has
    # precursor 100.
    expect_setequal(names(members), c("- 1,2,3,4", "- 5,6,7", "- 8"))
    # Each intra-file average counts once: medians of 900/600, 450/300; m/z
    # 200 from file 2 only.
    expect_equal(.average_peak_matrix(inter[members[["- 1,2,3,4"]]]),
        cbind(c(100, 150, 200), c(750, 375, 90), c(2, 2, 1), 2))
    expect_equal(.average_peak_matrix(inter[members[["- 5,6,7"]]]),
        cbind(c(60, 80), c(37.5, 85), 2, 2))
    ms2 <- inter[members[["- 1,2,3,4"]]]
    expect_identical(.average_input_ids(ms2)[[1L]],
        intra$spectrumId[Spectra::msLevel(intra) == 2L])
    expect_identical(sort(ms2$msPurityContributors[[1L]]$pid), 1:4)
    ms3 <- inter[members[["- 5,6,7"]]]
    expect_equal(ms3$msPurityPrecursorPath[[1L]],
        data.frame(msLevel = 3L, precursorMz = 150.0004))

    record <- processingHistory(ms2)[[ms2$msPurityProcessingId]]
    expect_identical(record$step, "averageInterFragSpectra")
    expect_setequal(record$output$spectra$spectrumId, inter$spectrumId)
    expect_setequal(record$input$spectra$spectrumId, intra$spectrumId)
    expect_false("msPurityKey" %in% names(record$output$spectra))
    expect_null(record$parameters$configuration$selected_feature_links)
})

test_that("inter-file links that contribute no retained peak are not used", {
    ms2 <- list(msLevel = 2L)
    pa <- averageIntraFragSpectra(.average_files_fixture(),
        param = .average_param(select = ms2))
    # m/z 200 (pid 4, file 2) appears in one input only and is removed.
    out <- averageInterFragSpectra(pa, rmp = TRUE,
        param = .average_param(select = ms2, minnum = 2))
    inter <- averagedSpectra(out, legacy = FALSE)
    inter <- inter[inter$av_level == "av_inter"]

    expect_equal(.average_peak_matrix(inter), cbind(c(100, 150), c(750, 375),
        2, 2))
    expect_identical(sort(.average_selected_links(inter)[[1L]]$pid), 1:4)
    expect_identical(sort(inter$msPurityContributors[[1L]]$pid), 1:3)

    id <- inter$spectrumId
    inputs <- averagedSpectra(pa, legacy = FALSE)
    file_ids <- inputs$spectrumId[order(as.integer(inputs$sample))]
    direct <- averagePeakMembers(inter, id, peak = 2L, depth = "immediate")
    expect_equal(direct$member, file_ids)
    expect_equal(direct$memberPeak, c(2L, 2L))
    leaf <- averagePeakMembers(out, id, peak = 2L, depth = "original")
    expect_equal(sort(leaf$pid), 1:3)
    expect_equal(leaf$memberPeak, rep(2L, 3L))

    bundle <- averageProvenance(out, id)
    expect_setequal(bundle$averages$spectrumId, c(id, file_ids))
    expect_setequal(bundle$measured$pid, 1:4)
    expect_setequal(bundle$featureLinks$linkId, averageMembers(inter, id,
        depth = "original")$linkId)
    expect_setequal(bundle$fragments$pid, 1:4)
    path <- tempfile()
    on.exit(unlink(path))
    saveRDS(bundle, path)
    snapshot <- readRDS(path)
    expect_identical(averagePeakMembers(snapshot$averages, id, peak = 2L,
        depth = "original"), leaf)
    expect_true(all(unique(snapshot$averages$msPurityProcessingId) %in%
        names(snapshot$processing)))

    flagged <- averageInterFragSpectra(pa, rmp = FALSE,
        param = .average_param(select = ms2, minnum = 2))
    sp <- averagedSpectra(flagged, legacy = FALSE)
    sp <- sp[sp$av_level == "av_inter"]
    edge <- averagePeakMembers(sp, sp$spectrumId, peak = 3L,
        depth = "original")
    expect_equal(edge$pid, 4L)
    expect_equal(edge$memberPeak, 1L)
})

test_that("inter-file peak edges survive reordering feature links", {
    param <- .average_param(select = list(msLevel = 2L))
    pa <- averageIntraFragSpectra(.average_files_fixture(), param = param)
    original <- pa@grped_df
    pa@grped_df <- original[rev(seq_len(nrow(original))), , drop = FALSE]
    out <- averageInterFragSpectra(pa, param = param)
    averages <- averagedSpectra(out, legacy = FALSE)
    final <- averages[averages$av_level == "av_inter"]
    # Only scan 4 has m/z 200; its occurrence keeps its original index.
    edge <- averagePeakMembers(final, final$spectrumId, peak = 3L,
        depth = "original")
    row <- which(original$pid == 4L)
    expect_equal(edge$pid, 4L)
    expect_identical(edge$linkId, original$linkId[row])
    expect_identical(edge$featureLinkIndex, row)
    expect_equal(.average_peak_matrix(final)[3L, ], c(200, 90, 1, 2))
})

test_that("inter-file settings must match the intra-file averages", {
    pa <- .average_files_fixture()
    expect_error(averageInterFragSpectra(pa, param = .average_param()),
        "requires intra-file averages computed with param")
    intra <- averageIntraFragSpectra(pa, param = .average_param())
    expect_error(averageInterFragSpectra(intra,
            param = list(other = .average_param())),
        "has no intra-file parameter set of the same name")
    expect_error(averageInterFragSpectra(intra,
            param = .average_param(combine = "collisionEnergy")),
        "must use the select and combine")
    # File 1 averaged MS3 precursors 150 and 150.0006 (4 ppm) together.
    expect_error(averageInterFragSpectra(intra, param = .average_param(
            ppm = 3)), "differ by more than the inter-file tolerance")
    expect_error(averageInterFragSpectra(intra),
        "at other levels \\(av_intra\\)")
})

test_that("each level replaces its own averages and keeps unrelated ones", {
    pa <- averageAllFragSpectra(.average_files_fixture(),
        param = .average_param())
    all_ids <- averagedSpectra(pa, legacy = FALSE)$spectrumId
    pa <- averageIntraFragSpectra(pa, param = .average_param())
    pa <- averageInterFragSpectra(pa, param = .average_param())
    levels <- function(pa) table(factor(averagedSpectra(pa,
        legacy = FALSE)$av_level, c("av_all", "av_intra", "av_inter")))
    # All-spectrum averages ignore files: MS2 {1-4}, MS3 {5, 6, 7} and {8}.
    expect_equal(as.vector(levels(pa)), c(3, 5, 3))

    again <- averageIntraFragSpectra(pa, param = .average_param())
    expect_equal(as.vector(levels(again)), c(3, 5, 0))
    expect_null(again@params$processing$active$averageInterFragSpectra)
    expect_length(again@av_inter_params, 0L)
    averages <- averagedSpectra(again, legacy = FALSE)
    expect_identical(averages$spectrumId[averages$av_level == "av_all"],
        all_ids)
    old <- averagedSpectra(pa, legacy = FALSE)$spectrumId
    fresh <- averages$spectrumId[averages$av_level == "av_intra"]
    expect_false(any(fresh %in% old))
    expect_setequal(fresh, processingHistory(again)[[
        unique(averages$msPurityProcessingId[averages$av_level == "av_intra"])
    ]]$output$spectra$spectrumId)
    expect_equal(as.vector(levels(averageAllFragSpectra(pa,
        param = .average_param()))), c(3, 5, 3))
})

test_that("averaging without param replaces its own level only", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
        "3_filterFragSpectra_pa.rds", package = "msPurity"))
    intra <- averageIntraFragSpectra(pa,
        param = AverageFragSpectraParam(combine = "collisionEnergy"))
    expect_error(averageAllFragSpectra(intra), "at other levels \\(av_intra\\)")

    scalar <- averageIntraFragSpectra(intra)
    expect_false(any(msPurity:::.pa_av_configured(scalar@avSpectra)))
    expect_type(averagedSpectra(scalar), "list")
    ids <- averagedSpectra(scalar, legacy = FALSE)$spectrumId
    expect_false(anyNA(ids))
    expect_false(any(ids %in% averagedSpectra(intra, legacy = FALSE)$spectrumId))
    # No columns of the replaced averages remain.
    expect_setequal(Spectra::spectraVariables(scalar@avSpectra),
        Spectra::spectraVariables(averageIntraFragSpectra(pa)@avSpectra))

    # A call with param that stored no averages leaves no settings behind.
    empty <- averageAllFragSpectra(pa,
        param = AverageFragSpectraParam(select = list(polarity = 0L)))
    expect_false(is.null(empty@av_all_params$parameter_sets))
    expect_null(averageAllFragSpectra(empty)@av_all_params$parameter_sets)
})

test_that("each file contributes at most one input to an inter-file average", {
    # File 1 MS3 precursors 150 and 150.0012 (8 ppm) are separate intra-file
    # averages at 5 ppm; at 10 ppm both agree with file 2's 150.0004.
    pa <- averageIntraFragSpectra(.average_files_fixture(
        ms3_precursors = c(150, 150.0012, 150.0004, 100)),
        param = .average_param())
    out <- averageInterFragSpectra(pa, param = .average_param(ppm = 10))
    inter <- averagedSpectra(out, legacy = FALSE)
    inter <- inter[inter$av_level == "av_inter"]
    expect_setequal(names(.average_members(inter)), c("- 1,2,3,4", "- 5,7",
        "- 6", "- 8"))
})

test_that("intra-file averages without passing peaks are not inputs", {
    # With minnum = 2, file 2's averages have no passing peak: its MS2 peaks
    # come from one spectrum each, and pid 7 and 8 are single spectra.
    for (remove in c(FALSE, TRUE)) {
        case <- paste("intra rmp =", remove)
        pa <- averageIntraFragSpectra(.average_files_fixture(), rmp = remove,
            param = .average_param(minnum = 2))
        intra <- averagedSpectra(pa, legacy = FALSE)
        expect_false("msPurityAverageStatus" %in%
            Spectra::spectraVariables(intra), info = case)
        out <- averageInterFragSpectra(pa, param = .average_param())
        inter <- averagedSpectra(out, legacy = FALSE)
        inter <- inter[inter$av_level == "av_inter"]
        members <- .average_members(inter)

        expect_true(setequal(names(members), c("- 1,2", "- 5,6")), info = case)
        expect_identical(.average_input_ids(inter)[[members[["- 1,2"]]]],
            intra$spectrumId[intra$sample == "1" & Spectra::msLevel(intra) == 2L],
            info = case)
    }
})

test_that("inter-file averaging ignores obsolete output status metadata", {
    pa <- averageIntraFragSpectra(.average_files_fixture(),
        param = .average_param())
    # A previous object may carry a stale label; actual peaks determine use.
    pa@avSpectra$msPurityAverageStatus <- rep("no_output_peaks_retained",
        length(pa@avSpectra))

    out <- averageInterFragSpectra(pa, param = .average_param())
    averages <- averagedSpectra(out, legacy = FALSE)
    inter <- averages[averages$av_level == "av_inter"]
    expect_setequal(names(.average_members(inter)),
        c("- 1,2,3,4", "- 5,6,7", "- 8"))
    expect_false("msPurityAverageStatus" %in%
        Spectra::spectraVariables(averages))
})

test_that("intra-file averages without eligible input peaks are not inputs", {
    filtered <- filterFragSpectra(.average_files_fixture(), ilim = 1e6,
        plim = NA, snr = 0, msLevel = 2:3)
    pa <- averageIntraFragSpectra(filtered, param = .average_param())
    intra <- averagedSpectra(pa, legacy = FALSE)
    expect_gt(length(intra), 0L)
    expect_true(all(lengths(Spectra::mz(intra)) == 0L))
    expect_true(all(vapply(.average_selected_links(intra), nrow,
        integer(1)) > 0L))
    # Cleanup also applies when the next call creates no new average.
    pa@avSpectra$msPurityAverageStatus <- rep("produced", length(intra))

    out <- averageInterFragSpectra(pa, param = .average_param())
    averages <- averagedSpectra(out, legacy = FALSE)
    expect_false(any(averages$av_level == "av_inter"))
    expect_false("msPurityAverageStatus" %in%
        Spectra::spectraVariables(averages))
    expect_identical(.average_selected_links(averages),
        .average_selected_links(intra))
    history <- processingParameters(out, "averageInterFragSpectra")
    configuration <- history[[length(history)]]$configuration
    expect_null(configuration$selected_feature_links)
    expect_null(configuration$inputs)
})

test_that("inter-file compatibility uses only passing intra-file clusters", {
    pa <- .average_files_fixture()
    pa@fragSpectra$dataOrigin <- paste0("/experiment/", pa@fragSpectra$dataOrigin)
    row <- match(6L, pa@fragSpectra$pid)
    mz <- Spectra::mz(pa@fragSpectra)
    intensity <- Spectra::intensity(pa@fragSpectra)
    mz[[row]] <- c(120, 130)
    intensity[[row]] <- c(40, 40)
    pa@fragSpectra$mz <- mz
    pa@fragSpectra$intensity <- intensity
    ms3 <- list(msLevel = 3L)

    for (remove in c(FALSE, TRUE)) {
        case <- paste("intra rmp =", remove)
        intra <- averageIntraFragSpectra(pa, rmp = remove,
            param = .average_param(select = ms3, ppm = 5, snr = 1.1))
        stored <- averagedSpectra(intra, legacy = FALSE)
        file1 <- stored[stored$sample == "1"]
        # pid 6 has SNR 1 at both peaks, so none passes. Keeping its failed
        # peaks must not make the inter-file 3 ppm path check see its 4 ppm
        # precursor offset. pid 5 and 7 contribute only m/z 80 (100 and 60).
        expect_identical(.average_selected_links(file1)[[1L]]$pid,
            c(5L, 6L), info = case)
        expect_identical(file1$msPurityContributors[[1L]]$pid,
            if (remove) 5L else c(5L, 6L), info = case)

        out <- averageInterFragSpectra(intra,
            param = .average_param(select = ms3, ppm = 3))
        inter <- averagedSpectra(out, legacy = FALSE)
        inter <- inter[inter$av_level == "av_inter"]
        expect_equal(length(inter), 1L, info = case)
        expect_equal(.average_peak_matrix(inter),
            matrix(c(80, 80, 2, 2), nrow = 1L), info = case)
        # pid 6 was selected upstream. Its peaks dropped out before inter.
        expect_setequal(.average_selected_links(inter)[[1L]]$pid, c(5L, 6L, 7L))
        expect_identical(inter$msPurityContributors[[1L]]$pid,
            c(5L, 7L), info = case)
        hit <- spectralMatching(inter, inter)
        leaves <- subset(hit$sourceMeasurements, side == "query")
        used <- pa@fragSpectra[match(c(5L, 7L), pa@fragSpectra$pid)]
        expect_setequal(paste(basename(leaves$dataOrigin), leaves$scanNumber),
            paste(basename(used$dataOrigin), used$acquisitionNum))
        expect_identical(hit$matchedResults$sharedMeasurement, TRUE)
        expect_equal(unname(inter$msPurityClusterContributors[[1L]]),
            list(c(5L, 7L)), info = case)
        expect_equal(inter$msPurityPrecursorPath[[1L]],
            data.frame(msLevel = 3L, precursorMz = 150.0002), info = case)
        edges <- averagePeakMembers(inter, inter$spectrumId, peak = 1L)
        # With failed peaks kept, m/z 80 is stored at position 2 in both
        # intra spectra. Removing their failed peaks moves it to position 1.
        expect_equal(edges$memberPeak, rep(if (remove) 1L else 2L, 2), info = case)
        expect_identical(.average_input_ids(inter)[[1L]],
            stored$spectrumId[c(which(stored$sample == "1"),
                which(vapply(.average_selected_links(stored), function(x)
                    identical(as.integer(x$pid), 7L), logical(1))))], info = case)
    }
})
