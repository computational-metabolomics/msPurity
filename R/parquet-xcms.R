# Features, chromatographic peaks and their links, mapped from the
# createDatabase() tables (c_peaks, c_peak_groups, c_peak_X_c_peak_group,
# c_peak_X_s_peak_meta or c_peak_group_X_s_peak_meta). These are the same for
# every xcms class, so every route gives the same results tables.

# c_peak_groups columns that are not per-sample values. Class-count columns
# are recognised separately.
.MZS_GROUP_FIXED <- c("grpid", "mz", "mzmin", "mzmax", "rt", "rtmin",
                      "rtmax", "npeaks", "peakidx", "ms_level", "grp_name")

#' The per-sample value columns of c_peak_groups.
#'
#' Names vary by route and class setup, so they are taken positionally: the
#' last `n_files` columns outside the fixed set. Any other column must be a
#' class count (an integer no larger than `npeaks`), else the table is refused.
#'
#' @noRd
.mzs_sample_columns <- function(groups, n_files) {
    cand <- setdiff(names(groups), .MZS_GROUP_FIXED)
    if (length(cand) < n_files)
        .mzs_abort("semantic", "c_peak_groups has ", length(cand),
                   " per-sample column(s) (", paste(cand, collapse = ", "),
                   ") for ", n_files, " source file(s) (R-098).")
    samples <- utils::tail(cand, n_files)
    extra <- setdiff(cand, samples)
    np <- suppressWarnings(as.numeric(groups$npeaks))
    for (c in extra) {
        v <- suppressWarnings(as.numeric(groups[[c]]))
        ok <- all(is.na(groups[[c]]) | (!is.na(v) & v == round(v) & v >= 0 &
                                        (is.na(np) | v <= np)))
        if (!ok)
            .mzs_abort("semantic", "c_peak_groups column '", c, "' is ",
                       "neither one of ", n_files, " per-sample value ",
                       "column(s) nor a class count (R-098).")
    }
    list(samples = samples, classes = extra)
}

#' Features, chromatographic peaks, their map, abundances, assays and
#' samples.
#'
#' @noRd
.mzs_feature_tables <- function(src, files) {
    cp <- src$c_peaks[order(src$c_peaks$cid), , drop = FALSE]
    run_of <- files$run_id[match(cp$fileid, files$fileid)]
    if (anyNA(run_of))
        .mzs_abort("semantic", "c_peaks names file(s) absent from fileinfo.")
    within <- stats::ave(seq_len(nrow(cp)), cp$fileid, FUN = seq_along)
    peak <- data.frame(
        chromatographic_peak_id_ = seq_len(nrow(cp)),
        run_id = run_of,
        peak_index_in_run = as.integer(within),
        exp_mass_to_charge = cp$mz,
        mass_to_charge_min = cp$mzmin,
        mass_to_charge_max = cp$mzmax,
        retention_time_in_seconds = cp$rt,
        retention_time_in_seconds_start = cp$rtmin,
        retention_time_in_seconds_end = cp$rtmax,
        signal_to_noise = if ("sn" %in% names(cp)) cp$sn else NA_real_,
        stringsAsFactors = FALSE)
    ## Remaining peak picker columns (into, intb, maxo, ...) pass through.
    for (c in setdiff(names(cp), c("cid", "mz", "mzmin", "mzmax", "rt",
                                   "rtmin", "rtmax", "fileid", "sn")))
        peak[[.mzs_passthrough(c)]] <- cp[[c]]
    if ("sn" %in% names(cp))
        peak$x_mspurity_sn <- cp$sn

    grp <- src$c_peak_groups[order(src$c_peak_groups$grpid), , drop = FALSE]
    feature <- data.frame(
        feature_id_ = seq_len(nrow(grp)),
        exp_mass_to_charge = grp$mz,
        mass_to_charge_min = grp$mzmin,
        mass_to_charge_max = grp$mzmax,
        retention_time_in_seconds = grp$rt,
        retention_time_in_seconds_start = grp$rtmin,
        retention_time_in_seconds_end = grp$rtmax,
        n_chromatographic_peaks = as.integer(round(grp$npeaks)),
        feature_name = grp$grp_name,
        stringsAsFactors = FALSE)
    if ("ms_level" %in% names(grp))
        feature$x_mspurity_ms_level <- suppressWarnings(
            as.integer(grp$ms_level))

    x <- src$c_peak_X_c_peak_group
    cpf <- data.frame(
        chromatographic_peak_id_ = match(x$cid, cp$cid),
        feature_id_ = match(x$grpid, grp$grpid),
        is_representative = x$bestpeak %in% 1L,
        x_mspurity_idi = x$idi)
    if (anyNA(cpf$chromatographic_peak_id_) || anyNA(cpf$feature_id_))
        .mzs_abort("semantic", "c_peak_X_c_peak_group names peaks or ",
                   "features absent from c_peaks or c_peak_groups.")
    cpf$peak_index_in_run <- peak$peak_index_in_run[
        cpf$chromatographic_peak_id_]

    sc <- .mzs_sample_columns(grp, nrow(files))
    abundance <- do.call(rbind, lapply(seq_along(sc$samples), function(j)
        data.frame(feature_id_ = feature$feature_id_, assay_id_ = j,
                   quantity_kind = "integrated_area",
                   value = suppressWarnings(as.numeric(grp[[sc$samples[j]]])),
                   value_source = "medret", stringsAsFactors = FALSE)))

    assay <- data.frame(assay_id_ = seq_len(nrow(files)),
                        assay_name = files$filename,
                        sample_id_ = seq_len(nrow(files)),
                        source = files$source, run_id = files$run_id,
                        stringsAsFactors = FALSE)
    sample <- data.frame(sample_id_ = seq_len(nrow(files)),
                         sample_name = sub("\\.[^.]*$", "", files$filename),
                         sample_class = files$class,
                         stringsAsFactors = FALSE)
    list(tables = list(feature = feature, chromatographic_peak = peak,
                       chromatographic_peak_feature = cpf,
                       abundance = abundance, assay = assay, sample = sample),
         sample_columns = sc, peaks = cp, groups = grp)
}

#' A source column name as `x_mspurity_<name>`, with `_id` rewritten to
#' `_ID` to avoid mzStack's reserved suffix.
#'
#' @noRd
.mzs_passthrough <- function(name) {
    name <- gsub("[^A-Za-z0-9_]+", "_", sub("^_+", "", name))
    name <- gsub("_id(?=_|$)", "_ID", name, perl = TRUE)
    paste0("x_mspurity_", name)
}

#' The spectrum_feature table: MS2 spectra associated with features.
#'
#' Scans linked via a chromatographic peak are run-scoped
#' `chromatographic_peak` links; scans linked via the full feature width
#' (frag4feature(useGroup = TRUE)) are `feature_width` links. Averaged
#' spectra are `manual` links to their feature.
#'
#' @param match_details live route only: data.frame with pid, cid, grpid,
#'     precurMtchMZ, precurMtchPPM, precurMtchRT.
#'
#' @noRd
.mzs_spectrum_feature <- function(src, scans, scan_refs, av, ftabs,
                                  match_details = NULL) {
    link <- .mzs_scan_links(src)
    s <- match(link$pid, scans$pid)
    grp_ids <- ftabs$groups$grpid
    by_peak <- !all(is.na(link$cid))
    sf <- data.frame(
        ms2_source = scan_refs$source[s], ms2_run_id = scan_refs$run_id[s],
        ms2_spectrum_id_ = scan_refs$spectrum_id_[s],
        ms2_native_id = scan_refs$native_id[s],
        feature_id_ = match(link$grpid, grp_ids),
        chromatographic_peak_id_ = if (by_peak)
            match(link$cid, ftabs$peaks$cid) else NA_integer_,
        link_mode = if (by_peak) "chromatographic_peak" else "feature_width",
        run_scoped = by_peak,
        precursor_mz_used = scans$precursorMZ[s],
        precursor_mz_error = NA_real_,
        precursor_retention_time = scans$precursorRT[s],
        precursor_purity = scans$inPurity[s],
        x_mspurity_scan_number = scans$acquisitionNum[s],
        stringsAsFactors = FALSE)
    ## The link's row id in the source database.
    link_id <- if (by_peak)
        src$c_peak_X_s_peak_meta$cXp_id[match(
            paste(link$pid, link$cid),
            paste(src$c_peak_X_s_peak_meta$pid,
                  src$c_peak_X_s_peak_meta$cid))]
    else src$c_peak_group_X_s_peak_meta$gXp_id[match(
        paste(link$pid, link$grpid),
        paste(src$c_peak_group_X_s_peak_meta$pid,
              src$c_peak_group_X_s_peak_meta$grpid))]
    if (!is.null(match_details) && nrow(match_details)) {
        key <- if (by_peak) paste(link$pid, link$cid) else
            paste(link$pid, link$grpid)
        mkey <- if (by_peak) paste(match_details$pid, match_details$cid) else
            paste(match_details$pid, match_details$grpid)
        k <- match(key, mkey)
        have <- !is.na(k)
        sf$precursor_mz_used[have] <- match_details$precurMtchMZ[k[have]]
        sf$precursor_mz_error[have] <- match_details$precurMtchPPM[k[have]]
        sf$precursor_retention_time[have] <-
            match_details$precurMtchRT[k[have]]
        ## precursor_mz_error is float32; keep the exact value too.
        sf$x_mspurity_precursor_mz_error_exact <- NA_real_
        sf$x_mspurity_precursor_mz_error_exact[have] <-
            match_details$precurMtchPPM[k[have]]
    }
    if (any(is.na(sf$feature_id_)))
        .mzs_abort("semantic", "An MS2 link names a feature absent from ",
                   "c_peak_groups.")
    if (nrow(av)) {
        self <- data.frame(
            ms2_source = "self", ms2_run_id = av$run_id,
            ms2_spectrum_id_ = av$spectrum_id_, ms2_native_id = av$native_id,
            feature_id_ = match(av$grpid, grp_ids),
            chromatographic_peak_id_ = NA_integer_,
            link_mode = "manual", run_scoped = av$method == "intra",
            precursor_mz_used = av$precursor_mz,
            precursor_mz_error = NA_real_,
            precursor_retention_time = av$retention_time,
            precursor_purity = av$inPurity,
            stringsAsFactors = FALSE)
        sf <- .mzs_rbind_fill(sf, self)
    }
    sf$association_id_ <- seq_len(nrow(sf))
    attr(sf, "link_id") <- link_id
    attr(sf, "link_column") <- if (by_peak) "cXp_id" else "gXp_id"
    attr(sf, "link_table") <- if (by_peak) "c_peak_X_s_peak_meta"
                              else "c_peak_group_X_s_peak_meta"
    sf
}
