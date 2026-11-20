# One file with two features; scan positions follow acquisition order.
# Feature 1: MS2 pid 1 and 2 (HCD 30) and 3 (CID 35); MS3 (CID 35) pid 6 from
# pid 1, pid 7 and 8 from pid 2, and pid 9 from the CID MS2 pid 3; MS4
# (HCD 30) pid 10 from pid 6, pid 11 and 12 from pid 7, pid 13 from pid 8 and
# pid 14 from pid 9. Feature 2: MS2 pid 4 and 5 (HCD 30).
.average_msn_fixture <- function(ms3_precursors = c(150, 150.0006, 100),
    ms4_precursors = c(80, 80.0002, 80.002, 80, 80),
    ms2_activation = c("HCD", "HCD", "CID"), reverse_links = FALSE,
    duplicate_root = FALSE, feature2_mz = list(200, c(200, 250))) {
    data <- S4Vectors::DataFrame(
        msLevel = c(1L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 2L, 2L, rep(4L, 5L)),
        acquisitionNum = 1:15, scanIndex = 1:15, rtime = 10 * (1:15),
        precScanNum = c(NA, 1L, 1L, 1L, 2L, 3L, 3L, 4L, 1L, 1L, 5L, 6L, 6L,
            7L, 8L),
        precursorMz = c(NA, 300, 300.001, 300, ms3_precursors, 150, 400, 400,
            ms4_precursors),
        activationType = c(NA, ms2_activation, rep("CID", 4L), "HCD", "HCD",
            rep("HCD", 5L)),
        collisionEnergy = c(NA, 30, 30, 35, rep(35, 4L), 30, 30, rep(30, 5L)),
        polarity = 1L, massAnalyzer = "Orbitrap", dataOrigin = "synthetic",
        spectrumId = paste0("scan=", 1:15))
    data$mz <- c(list(500, c(100, 150, 200), c(100, 150), c(110, 120),
        c(60, 80), c(60, 80), 50, 80), feature2_mz, rep(list(40), 5L))
    data$intensity <- list(1000, c(1000, 500, 50), c(800, 400), c(10, 20),
        c(50, 100), c(40, 120), 10, 90, 100, c(50, 20), 10, 20, 30, 40, 50)
    raw <- Spectra::Spectra(data, source = Spectra::MsBackendMemory())
    purity <- data.frame(pid = 1:5, fileid = 1L, seqNum = c(2L, 3L, 4L, 9L,
            10L), acquisitionNum = c(2L, 3L, 4L, 9L, 10L),
        retentionTime = c(20, 30, 40, 90, 100),
        precursorMZ = c(300, 300.001, 300, 400, 400), precursorScanNum = 1L,
        inPurity = 0.9)
    pa <- methods::new("purityA", fileList = c(synthetic = "synthetic"),
        cores = 1)
    pa@spectra <- msPurity:::.pa_msn_retained_spectra(raw, rep(1L, 15L),
        purity)
    metadata <- msPurity:::.pa_msn_spectra_metadata(pa@spectra)
    lineage <- msPurity:::.pa_msn_ancestry(metadata, "acquisitionNum")
    lineage$source_metadata <- metadata
    pa@params$msnLineage <- lineage
    roots <- data.frame(grpid = c(1L, 1L, 1L, 2L, 2L), pid = 1:5, sample = 1L,
        cid = 1:5, inPurity = 0.9, precurMtchMZ = purity$precursorMZ,
        retentionTime = purity$retentionTime)
    if (duplicate_root) roots <- rbind(roots, roots[1L, ])
    links <- msPurity:::.pa_msn_feature_links(roots, lineage)
    # Reversing every feature link, MS3/MS4 included, tests order independence.
    if (reverse_links) links <- links[rev(seq_len(nrow(links))), , drop = FALSE]
    pa@grped_df <- links[order(links$grpid), , drop = FALSE]
    peaks <- as.list(Spectra::peaksData(pa@spectra))
    linked <- split(peaks[match(pa@grped_df$pid, pa@spectra$pid)],
        as.character(pa@grped_df$grpid))
    pa@fragSpectra <- msPurity:::.pa_frag_from_grouped(pa, linked)
    msPurity:::.pa_sync_legacy(pa, puritydf = purity, grped_ms2 = linked)
}

.average_by_pid <- function(averages) {
    pids <- vapply(.average_selected_links(averages), function(x)
        paste(sort(x$pid), collapse = ","), character(1))
    stats::setNames(seq_along(averages), paste(averages$grpid, pids))
}

.average_peaks <- function(average) {
    peaks <- Spectra::peaksData(average,
        columns = Spectra::peaksVariables(average))[[1L]]
    as.matrix(peaks[order(peaks[, "mz"]), c("mz", "intensity", "count",
        "total"), drop = FALSE])
}

test_that("combining energy and activation preserves other acquisition axes", {
    cases <- list(massAnalyzer = c("Orbitrap", "ion trap"),
        scanFilterAnalyzer = c("FTMS", "ITMS"),
        collisionEnergyUnit = c("eV", "NCE"), polarity = c(1L, 0L))
    for (axis in names(cases)) {
        pa <- .average_msn_fixture()
        pa@grped_df <- pa@grped_df[pa@grped_df$pid %in% 1:2, , drop = FALSE]
        values <- rep(cases[[axis]][1L], length(pa@fragSpectra))
        values[pa@fragSpectra$pid == 2L] <- cases[[axis]][2L]
        pa@fragSpectra[[axis]] <- values
        out <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
            minfrac = 0, combine = c("activationType", "collisionEnergy")))
        s <- averagedSpectra(out, legacy = FALSE)

        expect_equal(length(s), 2L, info = axis)
        expect_equal(sort(s[[axis]]), sort(cases[[axis]]), info = axis)
        expect_identical(sort(vapply(.average_selected_links(s), nrow,
            integer(1))), c(1L, 1L), info = axis)
    }
})

test_that("provenance queries require the recorded processing activity", {
    out <- averageAllFragSpectra(.average_msn_fixture(),
        param = AverageFragSpectraParam(minfrac = 0))
    one <- averagedSpectra(out, legacy = FALSE)[1L]
    id <- one$spectrumId
    baseline <- averageMembers(one, id)
    expect_gt(nrow(baseline), 0L)
    for (value in list(missing = NA_character_, foreign = "unrelated",
            removed = NULL)) {
        broken <- one
        broken$msPurityProcessingId <- value
        expect_error(averageMembers(broken, id), "processing activity")
    }
    selected <- Spectra::selectSpectraVariables(one,
        setdiff(Spectra::spectraVariables(one), "av_level"))
    expect_identical(averageMembers(selected, id), baseline)
    selected$msPurityProcessingId <- NULL
    expect_error(averageMembers(selected, id), "spectrumId")
})

test_that("peak queries reject processing that changes their recorded positions", {
    out <- averageAllFragSpectra(.average_msn_fixture(),
        param = AverageFragSpectraParam(minfrac = 0))
    one <- averagedSpectra(out, legacy = FALSE)[1L]
    id <- one$spectrumId
    expect_equal(Spectra::mz(one)[[1L]], c(100, 150, 200))
    filtered <- Spectra::filterMzRange(one, mz = c(140, 210))
    expect_error(averagePeakMembers(filtered, id, peak = 1L), "unmodified")
    applied <- Spectra::applyProcessing(filtered)
    expect_error(averagePeakMembers(applied, id), "unmodified")
    expect_identical(averageMembers(applied, id), averageMembers(one, id))
    reordered <- Spectra::addProcessing(one,
        function(x) x[rev(seq_len(nrow(x))), , drop = FALSE])
    expect_error(averagePeakMembers(reordered, id), "unmodified")
    expect_equal(averagePeakMembers(one, id, peak = 2L)$memberPeak,
        c(2L, 2L))
})

test_that("parameter objects hold complete validated settings", {
    defaults <- as.list(AverageFragSpectraParam())
    scalar <- formals(methods::getGeneric("averageAllFragSpectra"))
    for (name in c("minfrac", "minnum", "ppm", "snr", "ra", "av", "sumi"))
        expect_equal(defaults[[name]], scalar[[name]], info = name)
    expect_setequal(names(defaults), c("select", "combine", "minfrac",
        "minnum", "ppm", "tolerance", "snr", "ra", "av", "sumi"))
    expect_identical(defaults$select, list())
    expect_identical(defaults$combine, character())
    expect_identical(defaults$tolerance, 0)
    expect_identical(AverageFragSpectraParam(ppm = 0, tolerance = 0.3)@ppm, 0)

    invalid <- list(
        combine = list(combine = "msLevel"),
        duplicate_combine = list(combine = c("collisionEnergy",
            "collisionEnergy")),
        minfrac = list(minfrac = 1.5), minnum = list(minnum = -1),
        ppm = list(ppm = -1), no_tolerance = list(ppm = 0, tolerance = 0),
        tolerance = list(tolerance = -0.1), snr = list(snr = NA_real_),
        av = list(av = "sum"), sumi = list(sumi = NA),
        select = list(select = list(analyzer = "Orbitrap")))
    for (name in names(invalid))
        expect_error(do.call(AverageFragSpectraParam, invalid[[name]]),
            info = name)

    pa <- .average_msn_fixture()
    expect_error(averageAllFragSpectra(pa,
            param = list(AverageFragSpectraParam())),
        "unique, nonempty names")
    expect_error(averageAllFragSpectra(pa,
            param = list(a = FilterFragSpectraParam())),
        "requires AverageFragSpectraParam")
})

test_that("explicit settings must agree with the parameter object", {
    pa <- .average_msn_fixture()
    param <- AverageFragSpectraParam(minfrac = 0)
    expect_error(averageAllFragSpectra(pa, minfrac = 0.5, param = param),
        "Explicit minfrac conflicts with param\\$minfrac")
    expect_error(averageAllFragSpectra(pa, sumi = FALSE,
            param = list(a = param, b = AverageFragSpectraParam(sumi = FALSE))),
        "Explicit sumi conflicts with param\\$a\\$sumi")
    expect_silent(averageAllFragSpectra(pa, minfrac = 0, param = param))
    expect_error(averageAllFragSpectra(pa, cores = 1.5, param = param),
        "cores.*positive integer")
    expect_error(averageAllFragSpectra(pa, rmp = NA, param = param),
        "rmp must be TRUE or FALSE")
})

test_that("one uniform MS2 group reproduces the scalar average", {
    pa <- .average_msn_fixture()
    ms2 <- pa
    ms2@grped_df <- pa@grped_df[pa@grped_df$msLevel == 2L, , drop = FALSE]
    scalar <- averagedSpectra(averageAllFragSpectra(ms2, minfrac = 0,
        rmp = TRUE), legacy = FALSE)
    typed <- averagedSpectra(averageAllFragSpectra(ms2, rmp = TRUE,
        param = AverageFragSpectraParam(minfrac = 0,
            combine = c("activationType", "collisionEnergy"))),
        legacy = FALSE)

    expect_identical(length(typed), length(scalar))
    for (i in seq_along(scalar)) {
        peaks <- function(x) Spectra::peaksData(x,
            columns = Spectra::peaksVariables(scalar))[[1L]]
        expect_equal(peaks(typed[i]), peaks(scalar[i]), info = i)
        expect_identical(typed$msPurityContributors[[i]],
            scalar$msPurityContributors[[i]], info = i)
        expect_identical(typed$msPurityClusterContributors[[i]],
            scalar$msPurityClusterContributors[[i]], info = i)
        expect_equal(typed$precursorMz[i], scalar$precursorMz[i], info = i)
        expect_equal(Spectra::rtime(typed)[i], Spectra::rtime(scalar)[i],
            info = i)
    }
})

test_that("conditions and precursor paths define separate averages", {
    pa <- .average_msn_fixture()
    out <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
        minfrac = 0, sumi = FALSE))
    averages <- averagedSpectra(out, legacy = FALSE)
    groups <- .average_by_pid(averages)

    # HCD and CID MS2 stay apart; MS3 pid 6 and 7 differ by 4 ppm, within the
    # default 5 ppm; pid 8 has another precursor; pid 9 descends from a CID
    # rather than HCD MS2. MS4 groups are checked in their own test.
    # The median of two observed peaks equals their mean.
    expect_setequal(names(groups), c("1 1,2", "1 3", "1 6,7", "1 8", "1 9",
        "1 10,11", "1 12", "1 13", "1 14", "2 4,5"))
    expect_false("msPurityKey" %in% Spectra::spectraVariables(averages))
    expect_equal(anyDuplicated(averages$spectrumId), 0L)
    expect_true(all(grepl("^[1-9][0-9]*$", averages$spectrumId)))
    expect_identical(pa@spectra$spectrumId, out@spectra$spectrumId)
    expect_identical(Spectra::msLevel(averages)[groups[c("1 1,2", "1 6,7",
            "2 4,5")]], c(2L, 3L, 2L))

    # Means of the observed peaks: (1000 + 800) / 2, (500 + 400) / 2, 50.
    expect_equal(unname(.average_peaks(averages[groups[["1 1,2"]]])),
        cbind(c(100, 150, 200), c(900, 450, 50), c(2, 2, 1), c(2, 2, 2)))
    # MS3: (50 + 40) / 2 at m/z 60 and (100 + 120) / 2 at m/z 80.
    ms3 <- averages[groups[["1 6,7"]]]
    expect_equal(unname(.average_peaks(ms3)),
        cbind(c(60, 80), c(45, 110), c(2, 2), c(2, 2)))
    expect_equal(ms3$precursorMz, 150.0003)
    expect_equal(ms3$msPurityPrecursorPath[[1L]],
        data.frame(msLevel = 3L, precursorMz = 150.0003))
    ancestors <- ms3$msPurityAncestorConditions[[1L]]
    expect_identical(ancestors$activationType, "HCD")
    expect_equal(ancestors$collisionEnergy, 30)
})

test_that("MS4 averages compare every stage after MS2", {
    pa <- .average_msn_fixture()
    out <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
        minfrac = 0, sumi = FALSE))
    averages <- averagedSpectra(out, legacy = FALSE)
    groups <- .average_by_pid(averages)

    # Stage precursors (MS3, MS4) at the default 5 ppm, allowing about
    # 0.00075 at m/z 150 and 0.0004 at m/z 80:
    # pid 10 (150, 80) and pid 11 (150.0006, 80.0002) agree at both stages;
    # pid 12 (150.0006, 80.002) differs at MS4 only; pid 13 (100, 80) differs
    # at MS3 only; pid 14 (150, 80) descends from the CID MS2 pid 3.
    expect_true(all(c("1 10,11", "1 12", "1 13", "1 14") %in% names(groups)))
    ms4 <- averages[groups[["1 10,11"]]]
    expect_identical(Spectra::msLevel(ms4), 4L)
    # Median of 10 and 20 at m/z 40, from two spectra.
    expect_equal(unname(.average_peaks(ms4)), cbind(40, 15, 2, 2))
    expect_equal(ms4$precursorMz, 80.0001)
    expect_equal(ms4$msPurityPrecursorPath[[1L]],
        data.frame(msLevel = 3:4, precursorMz = c(150.0003, 80.0001)))
    ancestors <- ms4$msPurityAncestorConditions[[1L]]
    expect_identical(ancestors$msLevel, 2:3)
    expect_identical(ancestors$activationType, c("HCD", "CID"))
    expect_equal(ancestors$collisionEnergy, c(30, 35))
    bundle <- averageProvenance(out, ms4$spectrumId)
    expect_setequal(bundle$measured$pid, c(1L, 2L, 6L, 7L, 10L, 11L))
    expect_setequal(bundle$featureLinks$pid, c(10L, 11L))
    row <- match(11L, bundle$lineage$pid)
    expect_equal(bundle$lineage$path[[row]], c(2L, 7L, 11L))
})

test_that("combine relaxes only the averaged stage", {
    pa <- .average_msn_fixture()
    out <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
        minfrac = 0, combine = c("activationType", "collisionEnergy")))
    groups <- .average_by_pid(averagedSpectra(out, legacy = FALSE))

    # The MS2 scans combine; MS3 pid 9 and MS4 pid 14 keep their CID rather
    # than HCD MS2 parent.
    expect_setequal(names(groups), c("1 1,2,3", "1 6,7", "1 8", "1 9",
        "1 10,11", "1 12", "1 13", "1 14", "2 4,5"))
})

test_that("activation names and accessions of one method group together", {
    pa <- .average_msn_fixture(ms2_activation = c("HCD",
        "beam-type collision-induced dissociation", "CID"))
    out <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
        select = list(activationType = "HCD", msLevel = 2L), minfrac = 0))
    averages <- averagedSpectra(out, legacy = FALSE)
    groups <- .average_by_pid(averages)
    expect_true("1 1,2" %in% names(groups))
    hcd <- averages[groups[["1 1,2"]]]
    expect_identical(hcd$activationType, "HCD")
    expect_identical(hcd$msPurityAcquisitionStatus[[1L]]$activationType,
        "uniform")

    # MS3 pid 6 and 7 descend from pid 1 and pid 2: the same parent method.
    ms3 <- .average_by_pid(averagedSpectra(averageAllFragSpectra(pa,
        param = AverageFragSpectraParam(select = list(msLevel = 3L),
            minfrac = 0)), legacy = FALSE))
    expect_true("1 6,7" %in% names(ms3))
})

test_that("a repeated feature link is averaged once per link", {
    pa <- .average_msn_fixture(duplicate_root = TRUE)
    out <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
        minfrac = 0, sumi = FALSE))
    averages <- averagedSpectra(out, legacy = FALSE)
    groups <- .average_by_pid(averages)

    # pid 1 is linked twice, so are its descendants pid 6 and pid 10.
    expect_true(all(c("1 1,1,2", "1 6,6,7", "1 10,10,11") %in% names(groups)))
    # m/z 100: peaks from pid 1 (twice, 1000) and pid 2 (800); median 1000.
    peaks <- .average_peaks(averages[groups[["1 1,1,2"]]])
    expect_equal(unname(peaks[peaks[, "mz"] == 100, ]), c(100, 1000, 3, 3))
})

test_that("MS3 precursors use the set's tolerances without chaining", {
    # pid 6, 7 and 8: MS3 precursors 150, 150.0006 and 150.0012, all from
    # HCD 30 MS2 parents. Neighbours differ by 0.0006 (4 ppm), the ends by
    # 0.0012 (8 ppm).
    chain <- c(150, 150.0006, 150.0012)
    groups <- function(..., reverse = FALSE) {
        pa <- .average_msn_fixture(ms3_precursors = chain,
            reverse_links = reverse)
        names(.average_by_pid(averagedSpectra(averageAllFragSpectra(pa,
            param = AverageFragSpectraParam(minfrac = 0, ...)),
            legacy = FALSE)))
    }
    for (reverse in c(FALSE, TRUE))
        expect_true(all(c("1 6,7", "1 8") %in% groups(ppm = 5,
            reverse = reverse)), info = paste("reverse links:", reverse))
    expect_true("1 6,7,8" %in% groups(ppm = 10))
    # Allowed differences: 3 ppm = 0.00045; 0.0002 + 3 ppm = 0.00065;
    # 0.0008 + 3 ppm = 0.00125; 0.0013 alone.
    expect_true(all(c("1 6", "1 7", "1 8") %in% groups(ppm = 3)))
    expect_true(all(c("1 6,7", "1 8") %in% groups(ppm = 3,
        tolerance = 0.0002)))
    expect_true("1 6,7,8" %in% groups(ppm = 3, tolerance = 0.0008))
    expect_true("1 6,7,8" %in% groups(ppm = 0, tolerance = 0.0013))
    # The MS2 precursors 300 and 300.001 differ by 3.3 ppm but are not
    # compared: frag4feature() linked both to the feature.
    expect_true("1 1,2" %in% groups(ppm = 3))
})

test_that("MS3 grouping respects ppm-only and absolute precursor boundaries", {
    # Compared with pid 6, pid 7 is just inside 5 ppm and pid 8 just outside.
    pa <- .average_msn_fixture(ms3_precursors = c(100, 100.0005, 100.000500008))
    ppm <- AverageFragSpectraParam(select = list(msLevel = 3L), minfrac = 0,
        ppm = 5, tolerance = 0)
    out <- averageAllFragSpectra(pa, param = ppm)
    groups <- .average_by_pid(averagedSpectra(out, legacy = FALSE))
    expect_setequal(names(groups), c("1 6,7", "1 8", "1 9"))

    # At the decimal absolute cutoff, subtraction rounds slightly upward.
    absolute <- AverageFragSpectraParam(select = list(msLevel = 3L), minfrac = 0,
        ppm = 0, tolerance = 0.000500008)
    out <- averageAllFragSpectra(pa, param = absolute)
    groups <- .average_by_pid(averagedSpectra(out, legacy = FALSE))
    expect_setequal(names(groups), c("1 6,7,8", "1 9"))
})

test_that("a precursor must agree with every member of an average", {
    order <- function(precursors) {
        pa <- .average_msn_fixture(ms3_precursors = precursors)
        names(.average_by_pid(averagedSpectra(averageAllFragSpectra(pa,
            param = AverageFragSpectraParam(minfrac = 0)), legacy = FALSE)))
    }
    # pid 8 (150.0012) agrees with pid 6 (150.0006) but not with pid 7 (150).
    expect_true(all(c("1 6,7", "1 8") %in% order(c(150.0006, 150, 150.0012))))
    # pid 8 (150.0006) agrees with both earlier averages, {6} and {7}, and
    # joins the first.
    expect_true(all(c("1 6,8", "1 7") %in% order(c(150, 150.0012, 150.0006))))
})

test_that("absolute and relative peak tolerances add", {
    # Feature 2: pid 4 has m/z 200.00 (100), pid 5 has 200.20 (50) and 250.
    pa <- .average_msn_fixture(feature2_mz = list(200, c(200.2, 250)))
    peaks <- function(ppm, tolerance) {
        out <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
            select = list(msLevel = 2L), minfrac = 0, sumi = FALSE,
            ppm = ppm, tolerance = tolerance))
        averages <- averagedSpectra(out, legacy = FALSE)
        unname(.average_peaks(averages[averages$grpid == "2"]))
    }
    apart <- cbind(c(200, 200.2, 250), c(100, 50, 20), 1, 2)
    joined <- cbind(c(200.1, 250), c(75, 20), c(2, 1), 2)
    # Allowed differences: 0.001; 0.15 + 0.001 = 0.151; 0.15 + 0.060 = 0.210;
    # 0.25 alone. Only the last two reach the 0.2 difference.
    expect_equal(peaks(ppm = 5, tolerance = 0), apart)
    expect_equal(peaks(ppm = 5, tolerance = 0.15), apart)
    expect_equal(peaks(ppm = 300, tolerance = 0.15), joined)
    expect_equal(peaks(ppm = 0, tolerance = 0.25), joined)
})

test_that("contributor indices are positions within each feature", {
    pa <- .average_msn_fixture()
    averages <- averagedSpectra(averageAllFragSpectra(pa,
        param = AverageFragSpectraParam(minfrac = 0)), legacy = FALSE)
    feature2 <- averages[averages$grpid == "2"]
    links <- which(pa@grped_df$grpid == 2L)

    contributors <- feature2$msPurityContributors[[1L]]
    expect_identical(contributors$pid, 4:5)
    expect_identical(contributors$contributorIndex, 1:2)
    expect_identical(contributors$featureLinkIndex, links)
    expect_identical(unname(feature2$msPurityClusterContributors[[1L]]),
        list(1:2, 2L))
})

test_that("selected links that contribute no retained peak are not used", {
    pa <- .average_msn_fixture()
    # MS2 of feature 1 combined: m/z 100 and 150 (pid 1 and 2) have two peaks;
    # pid 3 (CID, m/z 110 and 120) and m/z 200 have one and are removed.
    out <- averageAllFragSpectra(pa, rmp = TRUE,
        param = AverageFragSpectraParam(select = list(msLevel = 2L),
            combine = c("activationType", "collisionEnergy"), minnum = 2,
            minfrac = 0))
    average <- averagedSpectra(out, legacy = FALSE)
    average <- average[average$grpid == "1"]

    expect_identical(.average_selected_links(average)[[1L]]$pid, 1:3)
    inputs <- averageMembers(average, average$spectrumId)
    expect_identical(inputs$contributes[order(inputs$pid)], c(TRUE, TRUE, FALSE))
    expect_identical(average$msPurityContributors[[1L]]$pid, 1:2)
    expect_identical(average$activationType, "HCD")
    expect_equal(average$collisionEnergy, 30)
    expect_equal(average$precursorMz, 300.0005)
    expect_equal(Spectra::rtime(average), 25)
    expect_false("msPurityAverageStatus" %in%
        Spectra::spectraVariables(average))
})

test_that("averaged peak flags are exposed through Spectra", {
    pa <- .average_msn_fixture()
    param <- AverageFragSpectraParam(
        select = list(msLevel = 2L, activationType = "HCD"),
        minnum = 2, minfrac = 0, sumi = FALSE)
    flagged <- averagedSpectra(averageAllFragSpectra(pa, rmp = FALSE,
        param = param), legacy = FALSE)
    retained <- averagedSpectra(averageAllFragSpectra(pa, rmp = TRUE,
        param = param), legacy = FALSE)
    flagged <- flagged[flagged$grpid == "1"]
    retained <- retained[retained$grpid == "1"]

    columns <- c("mz", "intensity", "pass_flag")
    expect_equal(unname(as.matrix(Spectra::peaksData(flagged,
        columns = columns)[[1L]])),
        cbind(c(100, 150, 200), c(900, 450, 50), c(1, 1, 0)))
    expect_equal(unname(as.matrix(Spectra::peaksData(retained,
        columns = columns)[[1L]])),
        cbind(c(100, 150), c(900, 450), c(1, 1)))
    expect_false("msPurityAverageStatus" %in%
        Spectra::spectraVariables(flagged))
})

test_that("averages without input or retained peaks keep selected links", {
    pa <- .average_msn_fixture()
    none <- averageAllFragSpectra(pa,
        param = AverageFragSpectraParam(select = list(polarity = 0L)))
    expect_length(averagedSpectra(none, legacy = FALSE), 0L)
    history <- processingParameters(none, "averageAllFragSpectra")
    expect_identical(history[[length(history)]]$configuration$
        unused_parameter_sets, "default")
    expect_null(history[[length(history)]]$configuration$selected_feature_links)

    removed <- averageAllFragSpectra(pa, rmp = TRUE,
        param = AverageFragSpectraParam(select = list(msLevel = 3L),
            minnum = 5))
    averages <- averagedSpectra(removed, legacy = FALSE)
    expect_gt(length(averages), 0L)
    expect_false("msPurityAverageStatus" %in%
        Spectra::spectraVariables(averages))
    expect_true(all(lengths(Spectra::mz(averages)) == 0L))
    expect_true(all(is.na(averages$precursorMz)))
    expect_true(all(vapply(averages$msPurityContributors, nrow,
        integer(1)) == 0L))
    expect_true(all(vapply(.average_selected_links(averages), nrow,
        integer(1)) > 0L))
    record <- processingHistory(averages)[[averages$msPurityProcessingId[1L]]]
    expect_null(record$parameters$configuration$selected_feature_links)

    filtered <- filterFragSpectra(pa, ilim = 1e6, plim = NA, snr = 0,
        msLevel = 2:4)
    unfiltered <- averagedSpectra(averageAllFragSpectra(filtered,
        param = AverageFragSpectraParam()), legacy = FALSE)
    expect_gt(length(unfiltered), 0L)
    expect_true(all(lengths(Spectra::mz(unfiltered)) == 0L))
    expect_false("msPurityAverageStatus" %in%
        Spectra::spectraVariables(unfiltered))
    expect_true(all(vapply(.average_selected_links(unfiltered), nrow,
        integer(1)) > 0L))
})

test_that("named alternatives select and average independently", {
    pa <- .average_msn_fixture(duplicate_root = TRUE)
    out <- averageAllFragSpectra(pa, param = list(
        hcd = AverageFragSpectraParam(select = list(activationType = "HCD")),
        ms3 = AverageFragSpectraParam(select = list(msLevel = 3L),
            combine = "activationType")))
    averages <- averagedSpectra(out, legacy = FALSE)

    expect_identical(unique(averages$msPurityParameterSet), c("hcd", "ms3"))
    expect_true(all(Spectra::msLevel(averages[
        averages$msPurityParameterSet == "ms3"]) == 3L))
    history <- processingParameters(out, "averageAllFragSpectra")
    parameters <- history[[length(history)]]
    expect_identical(names(parameters$configuration$parameter_sets),
        c("hcd", "ms3"))
    expect_identical(parameters$effective_cores, 1L)
    assignments <- parameters$configuration$assignments
    expect_identical(assignments$parameter_set_names[
        assignments$pid == 1L][[1L]], "hcd")
    expect_identical(assignments$parameter_set_names[
        assignments$pid == 3L][[1L]], character())
    expect_null(parameters$configuration$selected_feature_links)

    # The first feature's HCD MS2 average selects the two occurrences of
    # pid 1 separately. Feature-link and within-feature positions are 1,2,4.
    hcd <- averages[averages$grpid == "1" &
        averages$msPurityParameterSet == "hcd" &
        Spectra::msLevel(averages) == 2L]
    expect_length(hcd, 1L)
    selected <- averageMembers(hcd, hcd$spectrumId)
    selected <- selected[order(selected$featureLinkIndex), ]
    expect_identical(selected$featureLinkIndex, c(1L, 2L, 4L))
    expect_identical(selected$contributorIndex, c(1L, 2L, 4L))
    expect_equal(selected$pid, c(1L, 2L, 1L))
    expect_equal(anyDuplicated(selected$linkId), 0L)
    edges <- averagePeakMembers(hcd, hcd$spectrumId, peak = 1L)
    expect_setequal(edges$linkId, selected$linkId)
    expect_false(any(c("msPuritySelectedContributors", "msPurityAverageInputs")
        %in% Spectra::spectraVariables(averages)))
})

test_that("repeated averaging allocates new spectrum IDs and preserves saved selections", {
    pa <- .average_msn_fixture()
    first <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
        select = list(msLevel = 2L, activationType = "HCD")))
    second <- averageAllFragSpectra(first, param = AverageFragSpectraParam(
        select = list(msLevel = 2L, activationType = "CID")))
    before <- averagedSpectra(first, legacy = FALSE)
    before <- before[before$grpid == "1"]
    after <- averagedSpectra(second, legacy = FALSE)

    expect_false(any(after$spectrumId %in% before$spectrumId))
    expect_false(identical(before$msPurityProcessingId, after$msPurityProcessingId))
    expect_equal(averageMembers(before, before$spectrumId)$pid, 1:2)
    expect_equal(averageMembers(second, after$spectrumId)$pid, 3L)
    expect_error(averageMembers(second, before$spectrumId))
})

test_that("independently averaged branches require their own graph", {
    pa <- .average_msn_fixture()
    param <- AverageFragSpectraParam(select = list(msLevel = 2L))
    left <- averagedSpectra(averageIntraFragSpectra(pa, param = param),
        legacy = FALSE)[1L]
    right <- averagedSpectra(averageIntraFragSpectra(pa, param = param),
        legacy = FALSE)[1L]
    expect_identical(left$spectrumId, right$spectrumId)
    expect_false(identical(left$msPurityProcessingId, right$msPurityProcessingId))
    mixed <- c(left, right)
    expect_error(averageMembers(mixed, left$spectrumId),
        "different|absent|ambiguous|multiple")
})

test_that("legacy access and export reject configured averages", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
        "3_filterFragSpectra_pa.rds", package = "msPurity"))
    configured <- averageAllFragSpectra(pa,
        param = AverageFragSpectraParam(combine = "collisionEnergy"))

    expect_error(averagedSpectra(configured), "legacy = FALSE")
    expect_error(averageIntraFragSpectra(configured),
        "cannot keep the averages computed with AverageFragSpectraParam")
    expect_error(averageInterFragSpectra(configured),
        "cannot keep the averages computed with AverageFragSpectraParam")
    out <- tempfile(fileext = ".db")
    expect_error(createDatabase(configured, NULL, outDir = dirname(out),
        dbName = basename(out)), "AverageFragSpectraParam is not supported")
    expect_false(file.exists(out))
    msp <- tempfile(fileext = ".msp")
    on.exit(unlink(msp), add = TRUE)
    expect_error(createMSP(configured, msp, method = "av_all"),
        "legacy = FALSE")
    expect_no_error(createMSP(configured, msp, method = "max"))

    scalar <- averageAllFragSpectra(configured)
    expect_type(averagedSpectra(scalar), "list")
    expect_false(any(msPurity:::.pa_av_configured(scalar@avSpectra)))
})

test_that("Parquet export holds configured averages with their provenance", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    suppressPackageStartupMessages(library(xcms))
    pa <- averageIntraFragSpectra(.golden_pa("3_filterFragSpectra_pa.rds"))
    configured <- averageAllFragSpectra(pa,
        param = AverageFragSpectraParam(combine = "collisionEnergy"))
    xe <- .golden_xcms("msms_only_xcmsnexp.rds")
    out <- tempfile("configured-")
    dir.create(out)
    path <- suppressMessages(createDatabase(configured, xe, outDir = out,
        dbName = "out.parquet", format = "parquet"))
    expect_false(any(validateParquet(path)$level == "MUST"))
    all <- .read_run(path, "av_all")
    expect_identical(sort(all$x_mspurity_spectrum_key),
        sort(configured@avSpectra$spectrumId[
            msPurity:::.pa_av_configured(configured@avSpectra)]))
    expect_true(all(grepl("^msPurity:av_all:[0-9]+:", all$id)))
    ## The intra-file averages computed without param keep their record.
    expect_true(length(grep("^av_intra_", vapply(.manifest(path)$runs,
        `[[`, "", "run_id"))) > 0L)
    mm <- .read_table(path, "merge_member")
    id <- all$spectrum_id_[1L]
    got <- mm[mm$merged_spectrum_id_ == id, ]
    ref <- averageMembers(configured, all$x_mspurity_spectrum_key[1L])
    expect_identical(got$member_rank[order(got$member_rank)], ref$rank)
    expect_identical(got$x_mspurity_contributes[order(got$member_rank)],
        ref$contributes)

    missing <- configured
    missing@params$averageGraph <- NULL
    expect_error(createDatabase(missing, xe, outDir = out, dbName = "m.parquet",
        format = "parquet"), "recorded provenance")
    changed <- configured
    changed@avSpectra <- Spectra::filterIntensity(changed@avSpectra,
        intensity = c(1e6, Inf))
    expect_error(createDatabase(changed, xe, outDir = out, dbName = "c.parquet",
        format = "parquet"), "Peak processing of averaged spectra")
})

test_that("configured averages keep provenance through subsetting and saving", {
    pa <- .average_msn_fixture()
    param <- AverageFragSpectraParam()
    out <- averageAllFragSpectra(pa, param = param)
    out <- averageIntraFragSpectra(out, param = param)
    out <- averageInterFragSpectra(out, param = param)
    averages <- averagedSpectra(out, legacy = FALSE)
    selected <- .average_selected_links(averages)
    expect_false(any(c("msPuritySelectedContributors", "msPurityAverageInputs")
        %in% Spectra::spectraVariables(averages)))
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(out, path)
    restored <- readRDS(path)
    expect_identical(processingHistory(restored), processingHistory(out))
    expect_identical(.average_selected_links(
        averagedSpectra(restored, legacy = FALSE)), selected)

    order <- rev(seq.int(1L, length(averages), by = 2L))
    subset <- averages[order]
    expect_identical(subset$spectrumId, averages$spectrumId[order])
    expect_identical(subset$msPurityAncestorConditions,
        averages$msPurityAncestorConditions[order])
    expect_identical(.average_selected_links(subset), selected[order])
    saveRDS(subset, path)
    restored_subset <- readRDS(path)
    expect_identical(.average_selected_links(restored_subset), selected[order])
    expect_identical(processingHistory(restored_subset),
        processingHistory(subset))
    inter <- subset[subset$av_level == "av_inter"]
    restored_inter <- restored_subset[
        restored_subset$av_level == "av_inter"]
    expect_identical(.average_input_ids(restored_inter),
        .average_input_ids(inter))
})

test_that("historical spectra do not acquire guessed peak edges", {
    historical <- Spectra::Spectra(S4Vectors::DataFrame(
        msLevel = 2L, spectrumId = "average:av_all:1:default:1",
        mz = list(100), intensity = list(1000)),
        source = Spectra::MsBackendMemory())
    historical$msPurityContributors <- list(data.frame(pid = 2L))
    expect_error(averageMembers(historical, historical$spectrumId), "unknown")
    expect_error(averagePeakMembers(historical, historical$spectrumId, peak = 1L),
        "unknown")
})
