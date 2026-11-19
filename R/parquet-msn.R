# MS2-MS4 scans, their acquisition lineage and their feature links, for the
# Parquet results of a live purityA object.
#
# msPurity's SQLite schema holds MS2 scans only, and the table set the
# converter maps (`src`) follows it. Everything above MS2, and what SQLite
# never held for MS2 (lineage, stage purity, activation, retention-time
# mapping details, filter flags of feature-linked scans, processing
# history), is taken from the live object here.

#' The MS2 core of a purityA object: what the SQLite schema can hold.
#'
#' Feature links and fragmentation spectra above MS2 are left out; the
#' spectra slot and the MSn lineage are kept, as they describe the source.
#'
#' @noRd
.pa_ms2_view <- function(pa) {
    g <- pa@grped_df
    if (nrow(g) && "msLevel" %in% names(g))
        pa@grped_df <- g[is.na(g$msLevel) | g$msLevel == 2L, , drop = FALSE]
    if (length(pa@fragSpectra))
        pa@fragSpectra <- pa@fragSpectra[Spectra::msLevel(pa@fragSpectra) ==
                                         2L]
    pa
}

#' The MS3 and MS4 scans a purityA object processed, in the columns of
#' `.mzs_scans()`, in pid order: all of them when purityA() assessed those
#' levels, else those frag4feature() linked or filterFragSpectra() kept.
#' Scans only retained alongside an MS2 assessment stay in the study.
#'
#' @noRd
.mzs_msn_scans <- function(pa) {
    sp <- pa@spectra
    if (!length(sp))
        return(NULL)
    sp <- sp[Spectra::msLevel(sp) > 2L]
    if (!length(sp))
        return(NULL)
    if (!any(pa@params$purityA$msLevel > 2L)) {
        g <- pa@grped_df
        used <- c(if (nrow(g) && "msLevel" %in% names(g))
                      g$pid[!is.na(g$msLevel) & g$msLevel > 2L],
                  if (length(pa@fragSpectra)) pa@fragSpectra$pid[
                      Spectra::msLevel(pa@fragSpectra) > 2L])
        sp <- sp[sp$pid %in% used]
        if (!length(sp))
            return(NULL)
    }
    have <- Spectra::spectraVariables(sp)
    col <- function(v) if (v %in% have) sp[[v]] else rep(NA, length(sp))
    out <- data.frame(
        pid = as.integer(sp$pid), fileid = as.integer(as.character(sp$fileid)),
        seqNum = as.integer(col("seqNum")),
        acquisitionNum = as.integer(col("acquisitionNum")),
        precursorScanNum = as.integer(col("precursorScanNum")),
        precursorMZ = as.numeric(col("precursorMZ")),
        precursorIntensity = as.numeric(col("precursorIntensity")),
        retentionTime = as.numeric(col("retentionTime")),
        spectrum_type = "scan", stringsAsFactors = FALSE)
    out[order(out$pid), , drop = FALSE]
}

#' Per scan: level, native id, dataset key, lineage, stage purity and the
#' acquisition labels msPurity used, from the live object, by pid.
#'
#' @noRd
.mzs_scan_annotations <- function(pa, pids) {
    sp <- pa@spectra
    k <- match(pids, sp$pid)
    have <- Spectra::spectraVariables(sp)
    col <- function(v, na = NA) if (v %in% have) sp[[v]][k] else
        rep(na, length(pids))
    out <- data.frame(
        pid = pids, ms_level = as.integer(Spectra::msLevel(sp)[k]),
        native_id = as.character(col("spectrumId", NA_character_)),
        dataset_spectrum_id_ = as.numeric(col("spectrum_id_", NA_real_)),
        dataset_run_id = as.character(col("run_id", NA_character_)),
        stage_purity = as.numeric(col("msPurityStagePurity", NA_real_)),
        purity_status = as.character(col("msPurityPurityStatus",
                                         NA_character_)),
        activation_type = as.character(col("activationType", NA_character_)),
        activation_type_basis = as.character(col("activationTypeBasis",
                                                 NA_character_)),
        scan_filter_analyzer = as.character(col("scanFilterAnalyzer",
                                                NA_character_)),
        stringsAsFactors = FALSE)
    lin <- pa@params$msnLineage$scans
    j <- if (is.data.frame(lin)) match(pids, lin$pid) else
        rep(NA_integer_, length(pids))
    ## The first of the given names the lineage has (objects saved before
    ## the lineage columns took their current names use the others).
    lcol <- function(v) {
        v <- if (is.data.frame(lin)) intersect(v, names(lin)) else character()
        if (length(v)) lin[[v[1L]]][j] else rep(NA, length(pids))
    }
    out$parent_pid <- as.integer(lcol(c("parentPid", "parent_pid")))
    out$root_pid <- as.integer(lcol(c("rootPid", "root_pid")))
    out$parent_status <- as.character(lcol(c("parentStatus", "parent_status")))
    out$path_status <- as.character(lcol(c("pathStatus", "path_status")))
    out$parent_reference_type <- as.character(lcol("parentReferenceType"))
    out
}

#' Peak flags of the fragmentation spectra `pids`, in the columns
#' `.mzs_scan_tables()` builds from s_peaks, from the live fragSpectra.
#'
#' @noRd
.mzs_live_scan_peaks <- function(pa, pids) {
    fs <- pa@fragSpectra
    if (!length(fs) || !length(pa@filter_frag_params) || !length(pids))
        return(NULL)
    fs <- fs[fs$pid %in% pids]
    if (!length(fs))
        return(NULL)
    pv <- Spectra::peaksVariables(fs)
    cols <- intersect(c("mz", "intensity", "snr", "ra",
                        "intensity_pass_flag", "ra_pass_flag",
                        "snr_pass_flag"), pv)
    pk <- as.list(Spectra::peaksData(fs, columns = cols))
    n <- vapply(pk, nrow, integer(1L))
    if (!sum(n))
        return(NULL)
    v <- do.call(rbind, lapply(pk, as.data.frame))
    flag <- if ("purity_pass_flag" %in% Spectra::spectraVariables(fs))
        as.logical(fs$purity_pass_flag) else rep(NA, length(fs))
    get <- function(c, f) if (c %in% names(v)) f(v[[c]]) else
        f(rep(NA, nrow(v)))
    data.frame(pid = rep(as.integer(fs$pid), n),
               mz = as.numeric(v$mz), i = as.numeric(v$intensity),
               snr = get("snr", as.numeric), ra = get("ra", as.numeric),
               purity_pass_flag = rep(flag, n),
               intensity_pass_flag = get("intensity_pass_flag", as.logical),
               ra_pass_flag = get("ra_pass_flag", as.logical),
               snr_pass_flag = get("snr_pass_flag", as.logical),
               stringsAsFactors = FALSE)
}

#' Retention-time mapping details of each MS2 feature link, which
#' frag4feature() records and the SQLite schema has no place for.
#'
#' @noRd
.MZS_LINK_DETAILS <- c(
    precurMtchRTadj = "precursor_retention_time_adjusted",
    rtMapStatus = "rt_map_status", rtApexDist = "rt_apex_distance",
    matchRank = "match_rank", precurMzSource = "precursor_mz_source")

#' The x_mspurity_msn_feature table: grped_df rows above MS2.
#'
#' @param scans the exported scans (pid, scan_annotation_id_ order).
#'
#' @noRd
.mzs_msn_feature <- function(pa, scans, scan_ids, ftabs) {
    g <- pa@grped_df
    if (!nrow(g) || !"msLevel" %in% names(g))
        return(NULL)
    g <- g[!is.na(g$msLevel) & g$msLevel > 2L, , drop = FALSE]
    if (!nrow(g))
        return(NULL)
    gcol <- function(v, f) if (v %in% names(g)) f(g[[v]]) else
        f(rep(NA, nrow(g)))
    by_peak <- pa@f4f_link_type == "individual" && "cid" %in% names(g)
    k <- match(g$pid, scans$pid)
    if (anyNA(k))
        .mzs_abort("semantic", "An MSn feature link names a scan that is ",
                   "not in the purityA object.")
    root <- gcol("msPurityRootPid", as.integer)
    out <- data.frame(
        scan_annotation_id_ = scan_ids[k],
        root_scan_annotation_id_ = scan_ids[match(root, scans$pid)],
        feature_id_ = match(as.integer(as.character(g$grpid)),
                            ftabs$groups$grpid),
        chromatographic_peak_id_ = if (by_peak)
            match(as.integer(g$cid), ftabs$peaks$cid) else NA_integer_,
        link_mode = if (by_peak) "chromatographic_peak" else "feature_width",
        ms_level = as.integer(g$msLevel),
        precursor_mz_used = gcol("precurMtchMZ", as.numeric),
        precursor_retention_time = gcol("precurMtchRT", as.numeric),
        precursor_retention_time_adjusted = gcol("precurMtchRTadj",
                                                 as.numeric),
        retention_time = gcol("retentionTime", as.numeric),
        rt_map_status = gcol("rtMapStatus", as.character),
        rt_apex_distance = gcol("rtApexDist", as.numeric),
        match_rank = gcol("matchRank", as.integer),
        precursor_mz_source = gcol("precurMzSource", as.character),
        stringsAsFactors = FALSE)
    if (anyNA(out$feature_id_))
        .mzs_abort("semantic", "An MSn feature link names a feature absent ",
                   "from the XCMS object.")
    out <- out[order(out$feature_id_, out$scan_annotation_id_), , drop = FALSE]
    out$msn_feature_id_ <- seq_len(nrow(out))
    rownames(out) <- NULL
    out
}

#' Processing history records by step, as JSON, for tool_provenance.
#' Imported records, which only restate the slots of an older object, are
#' left out.
#'
#' @noRd
.mzs_history_json <- function(pa) {
    records <- pa@params$processing$records
    if (!length(records))
        return(list())
    records <- Filter(function(r) !identical(r$status, "imported"), records)
    steps <- unique(vapply(records, `[[`, "", "step"))
    out <- lapply(steps, function(s) {
        x <- Filter(function(r) identical(r$step, s), records)
        as.character(.mzs_to_json(.mzs_json_safe(unname(x)), pretty = FALSE))
    })
    names(out) <- steps
    out
}

# Acquisition conditions msPurity reads from a study, by Spectra variable.
.MZS_STUDY_CONDITIONS <- c(activationType = "activation method",
                           collisionEnergy = "collision energy",
                           polarity = "polarity",
                           isolationWindowTargetMz = "isolation window")

#' Loss ledger rows for the scans above MS2 and the acquisition conditions
#' a study does not record.
#'
#' @noRd
.mzs_msn_losses <- function(files, scans, pa) {
    out <- list()
    msn <- scans$ms_level > 2L
    if (!is.null(pa) && any(msn) && !length(pa@filter_frag_params))
        out[[length(out) + 1L]] <- list(
            construct = "MS3/MS4 scan peaks",
            disposition = "no_such_concept_in_target",
            affected_rows = sum(msn),
            reason = paste("Raw scans are referenced, not copied (R-009);",
                           "no filterFragSpectra() flags were set."))
    lin <- if (!is.null(pa)) pa@params$msnLineage$scans
    if (is.data.frame(lin) && "path" %in% names(lin))
        out[[length(out) + 1L]] <- list(
            construct = "msnLineage.path",
            disposition = "dropped_by_configuration",
            affected_rows = nrow(lin), target_ref = "x_mspurity_scan",
            reason = paste("The chain of parents of each scan, recomputable",
                           "from parent_scan_annotation_id_."))
    study <- attr(files, "study")
    if (!is.null(study) && nrow(scans)) {
        frag <- .mzs_read_spectra_vars(study, attr(files, "study_manifest"),
            c("msLevel", names(.MZS_STUDY_CONDITIONS)),
            unique(files$run_id))
        frag <- frag[frag$msLevel %in% 2:4, , drop = FALSE]
        for (v in names(.MZS_STUDY_CONDITIONS)) {
            known <- !is.na(frag[[v]])
            if (!nrow(frag) || all(known))
                next
            out[[length(out) + 1L]] <- list(
                construct = paste0("study.", v),
                disposition = "not_present_in_source",
                affected_rows = sum(!known),
                target_ref = "x_mspurity_scan.ms2_spectrum_id_",
                reason = paste0("The study records no ",
                                .MZS_STUDY_CONDITIONS[[v]], " for ",
                                if (any(known)) "some of " else "",
                                "its fragmentation spectra; it is reported ",
                                "as unavailable, not inferred",
                                if (v == "activationType")
                                    paste0(" (activation_type_basis names ",
                                           "any filter-string label msPurity ",
                                           "used)") else "", "."))
        }
    }
    out
}

#' The MS2 rows of an x_mspurity_scan table; a table written before it held
#' other levels has no ms_level and is all MS2.
#'
#' @noRd
.mzs_ms2_scans <- function(sc) {
    if (is.null(sc) || is.null(sc$ms_level))
        return(sc)
    sc[sc$ms_level %in% c(NA, 2L), , drop = FALSE]
}
