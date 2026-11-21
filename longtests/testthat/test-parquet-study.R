context("Raw MS1-MS4 studies converted by mzMLToParquet()")

test_that("both engines store the D02 spectra as MsBackendMzR reads them", {
    skip_if_no_study()
    vars <- c("msLevel", "acquisitionNum", "scanIndex", "precScanNum",
              "precursorMz", "precursorIntensity", "precursorCharge",
              "collisionEnergy", "isolationWindowLowerMz",
              "isolationWindowTargetMz", "isolationWindowUpperMz",
              "polarity", "spectrumId", "centroided", "filterString",
              "dataOrigin")
    for (case in seq_len(nrow(.msn_standard_cases))) {
        f <- .msn_standard_file(.msn_standard_cases$compound[case],
                                .msn_standard_cases$kind[case])
        ref <- Spectra::Spectra(f, source = Spectra::MsBackendMzR())
        for (engine in c("spectra", "mzr")) {
            info <- paste(basename(f), engine)
            sp <- .msn_study_spectra(.msn_study(f, engine))
            expect_identical(unique(sp$run_id),
                             sub("\\.mzML\\.gz$", "", basename(f)), info = info)
            for (v in vars)
                expect_identical(Spectra::spectraData(sp, v)[[v]],
                                 Spectra::spectraData(ref, v)[[v]],
                                 info = paste(info, v))
            expect_equal(Spectra::rtime(sp), Spectra::rtime(ref),
                         tolerance = 1e-12, info = info)
            expect_identical(sp$peaksCount, lengths(sp), info = info)
            expect_equal(as.list(Spectra::peaksData(sp)),
                         as.list(Spectra::peaksData(ref)), info = info)
            ## Each MSn spectrum's parent is the scan precScanNum names.
            msn <- Spectra::msLevel(sp) > 1L
            parent <- match(Spectra::precScanNum(sp)[msn],
                            Spectra::acquisitionNum(sp))
            expect_false(anyNA(parent), info = info)
            expect_true(all(Spectra::msLevel(sp)[parent] ==
                            Spectra::msLevel(sp)[msn] - 1L), info = info)
            ## The activation the mzML records is the method its filter
            ## string names for the spectrum's own stage.
            meta <- msPurity:::.pa_filter_metadata(sp)
            parsed <- msPurity:::.pa_filter_metadata(
                data.frame(msLevel = Spectra::msLevel(ref),
                           filterString = ref$filterString))
            expect_identical(meta$activationType, parsed$activationType,
                             info = info)
            if (any(msn))
                expect_identical(unique(attr(meta, "basis")$activationType[
                    msn]), "recorded", info = info)
        }
    }
})

test_that("the MSn workflow gives the same results on the study", {
    skip_if_no_study()
    skip_if_not_installed("xcms")
    f <- .d02_copy()
    xe <- .msn_standard_feature("tryptophan", "MSn")$xe
    ref <- .msn_workflow(f, xe)
    for (engine in c("spectra", "mzr")) {
        g <- .d02_copy()
        study <- tempfile("study-")
        suppressMessages(MsBackendParquet::mzMLToParquet(
            g, study, engine = engine, verbose = FALSE))
        unlink(g)                    # nothing to fall back on
        pa <- .msn_workflow(.msn_study_spectra(study), xe)
        expect_identical(pa@params$purityA$study$path, normalizePath(study))
        drop <- c("filename")
        expect_equal(purityTable(pa)[setdiff(names(purityTable(pa)), drop)],
                     purityTable(ref)[setdiff(names(purityTable(ref)), drop)],
                     info = engine)
        lineage <- function(x) {
            s <- x@params$msnLineage$scans
            s <- s[order(s$pid), c("pid", "parentPid", "rootPid",
                                   "parentStatus", "pathStatus",
                                   "parentReferenceType")]
            rownames(s) <- NULL
            s
        }
        expect_identical(lineage(pa), lineage(ref), info = engine)
        expect_identical(pa@spectra$activationType, ref@spectra$activationType)
        expect_identical(unique(na.omit(ref@spectra$activationTypeBasis)),
                         "filter_string")
        expect_identical(unique(na.omit(pa@spectra$activationTypeBasis)),
                         "recorded")
        cols <- c("grpid", "pid", "msLevel", "msPurityRootPid",
                  "precurMtchRT", "precurMtchMZ", "inPurity",
                  "msPurityStagePurity", "purity_pass_flag")
        expect_equal(pa@grped_df[cols], ref@grped_df[cols], info = engine)
        expect_identical(Spectra::msLevel(pa@fragSpectra),
                         Spectra::msLevel(ref@fragSpectra))
        pv <- Spectra::peaksVariables(ref@fragSpectra)
        expect_equal(as.list(Spectra::peaksData(pa@fragSpectra, pv)),
                     as.list(Spectra::peaksData(ref@fragSpectra, pv)),
                     info = engine)
        expect_equal(pa@fragSpectra$msPurityEffectiveParameterId,
                     ref@fragSpectra$msPurityEffectiveParameterId)
        expect_equal(allFragSpectra(pa, legacy = TRUE),
                     allFragSpectra(ref, legacy = TRUE), info = engine)
        expect_equal(getScanPeaks(pa), getScanPeaks(ref), info = engine)
    }
})

test_that("Parquet results reference every MS2-MS4 scan and its parents", {
    skip_if_no_study()
    skip_if_not_installed("xcms")
    f <- .d02_copy()
    xe <- .msn_standard_feature("tryptophan", "MSn")$xe
    study <- tempfile("study-")
    suppressMessages(MsBackendParquet::mzMLToParquet(f, study,
                                                     verbose = FALSE))
    pa_mz <- .msn_workflow(f, xe)
    pa_pq <- .msn_workflow(.msn_study_spectra(study), xe)
    out <- tempfile("results-")
    dir.create(out)
    write <- function(pa, name, ...) suppressMessages(createDatabase(
        pa, xe, outDir = out, dbName = name, format = "parquet", ...))
    p_mz <- write(pa_mz, "mz.parquet", study = study)
    expect_message(p_pq <- createDatabase(pa_pq, xe, outDir = out,
                                          dbName = "pq.parquet",
                                          format = "parquet"),
                   "Linking scans to the Parquet dataset purityA\\(\\) read")
    p_ext <- write(pa_pq, "ext.parquet", study = NA)
    for (p in c(p_mz, p_pq, p_ext))
        expect_false(any(validateParquet(p)$level == "MUST"), info = p)

    sc_mz <- .read_table(p_mz, "x_mspurity_scan")
    sc_pq <- .read_table(p_pq, "x_mspurity_scan")
    expect_identical(as.vector(table(sc_pq$ms_level)),
                     as.vector(table(Spectra::msLevel(pa_pq@spectra))))
    expect_identical(unique(sc_pq$link_basis), "dataset")
    expect_identical(unique(sc_mz$link_basis), "native_id")
    same <- c("ms_level", "ms2_source", "ms2_run_id", "ms2_spectrum_id_",
              "ms2_native_id", "parent_source", "parent_run_id",
              "parent_spectrum_id_", "parent_native_id", "root_spectrum_id_",
              "parent_scan_annotation_id_", "root_scan_annotation_id_",
              "path_status", "stage_purity", "in_purity", "activation_type")
    expect_identical(sc_pq[same], sc_mz[same])
    expect_identical(sc_mz$activation_type_basis[sc_mz$ms_level > 1L],
                     rep("filter_string", nrow(sc_mz)))
    expect_identical(unique(sc_pq$activation_type_basis), "recorded")

    # Datasets written before averaging still expose their measured inputs.
    scans <- readProcessedSpectra(p_pq, msLevel = 3L, includeRaw = TRUE)
    expect_length(scans$processed, 0L)
    expect_setequal(scans$raw$study$spectrum_id_, sc_pq$ms2_spectrum_id_[sc_pq$ms_level == 3L])
    pr <- readProcessedProvenance(p_pq, msLevel = 3L)
    expect_true(all(pr$msnFeatureLinks$ms_level == 3L))
    expect_setequal(pr$measured$ms2_spectrum_id_, scans$raw$study$spectrum_id_)

    ## Every scan and every parent is a study spectrum one level down.
    sp <- .msn_study_spectra(study)
    lvl <- Spectra::msLevel(sp)[match(sc_pq$ms2_spectrum_id_,
                                      as.numeric(sp$spectrum_id_))]
    expect_identical(lvl, sc_pq$ms_level)
    plvl <- Spectra::msLevel(sp)[match(sc_pq$parent_spectrum_id_,
                                       as.numeric(sp$spectrum_id_))]
    expect_identical(plvl, sc_pq$ms_level - 1L)
    rlvl <- Spectra::msLevel(sp)[match(sc_pq$root_spectrum_id_[
        sc_pq$ms_level > 2L], as.numeric(sp$spectrum_id_))]
    expect_true(all(rlvl == 2L))

    ## Without a study, scans are external native-id references.
    sc_ext <- .read_table(p_ext, "x_mspurity_scan")
    expect_true(all(is.na(sc_ext$ms2_spectrum_id_)))
    expect_identical(sc_ext$ms2_native_id, sc_pq$ms2_native_id)
    expect_identical(sc_ext$parent_native_id[sc_ext$ms_level > 2L],
                     sc_pq$parent_native_id[sc_pq$ms_level > 2L])

    ## MSn feature links, filter flags and the processing history.
    mf <- .read_table(p_pq, "x_mspurity_msn_feature")
    g <- pa_pq@grped_df
    expect_identical(as.vector(table(mf$ms_level)),
                     as.vector(table(g$msLevel[g$msLevel > 2L])))
    peaks <- .read_table(p_pq, "x_mspurity_scan_peak")
    flagged <- sc_pq$ms_level[match(peaks$scan_annotation_id_,
                                    sc_pq$scan_annotation_id_)]
    expect_true(all(c(2L, 3L, 4L) %in% flagged))
    tp <- .read_table(p_pq, "tool_provenance")
    hist <- tp$x_mspurity_history[tp$step_name == "filterFragSpectra"]
    expect_identical(jsonlite::fromJSON(hist, simplifyVector = FALSE)[[1L]]
                     $step, "filterFragSpectra")
    sf <- .read_table(p_pq, "spectrum_feature")
    expect_true("x_mspurity_rt_map_status" %in% names(sf))

    ## Spectral matching queries the MS2 scans only.
    q <- msPurity:::.mzs_query_spectra(p_mz, MsBackendParquet::readManifest(
        p_mz), list(q_spectraTypes = "scan", q_spectraFilter = FALSE,
                    q_purity = NA, q_pol = NA, q_rtrange = NA,
                    q_xcmsGroups = NA), NULL)
    expect_identical(sort(q$meta$query_spectrum_id_),
                     sort(sc_mz$ms2_spectrum_id_[sc_mz$ms_level == 2L]))

    raw2 <- readProcessedSpectra(p_mz, msLevel = 2L, includeRaw = TRUE)$raw$study
    first <- q$meta$query_spectrum_id_[1L]
    pk <- Spectra::peaksData(raw2[match(first, raw2$spectrum_id_)])[[1L]]
    expect_equal(q$peaks[[1L]]$mz, unname(pk[, "mz"]))
    expect_equal(q$peaks[[1L]]$i, unname(pk[, "intensity"]))
    filtered <- msPurity:::.mzs_query_spectra(p_mz, .manifest(p_mz),
        list(q_spectraTypes = "scan", q_spectraFilter = TRUE, q_purity = NA,
             q_pol = NA, q_rtrange = NA, q_xcmsGroups = NA), NULL)
    scan_id <- sc_mz$scan_annotation_id_[match(first, sc_mz$ms2_spectrum_id_)]
    ann <- .read_table(p_mz, "x_mspurity_scan_peak")
    ann <- ann[ann$scan_annotation_id_ == scan_id, ]
    ann <- ann[order(ann$peak_rank), ]
    pass <- with(ann, purity_pass_flag & intensity_pass_flag & ra_pass_flag & snr_pass_flag)
    expect_equal(filtered$peaks[[1L]]$mz, ann$mz[which(pass)])
    expect_equal(filtered$peaks[[1L]]$i, ann$intensity[which(pass)])

    ## The SQLite writer holds MS2 only and says so.
    expect_error(suppressWarnings(createDatabase(pa_pq, xe, outDir = out,
        dbName = "x.sqlite")), "not supported by the SQLite writer")
})
