# Averaged spectra in the Parquet results: those computed with
# AverageFragSpectraParam, which the SQLite schema cannot hold, and the
# members and peaks every average was made from, as its recorded provenance
# (pa@params$averageGraph) gives them.

#' Whether averaged spectra no longer hold the peaks averaging stored.
#'
#' @param n_peaks the peak counts the provenance records.
#'
#' @noRd
.mzs_av_peaks_changed <- function(sp, n_peaks) {
    if (!length(sp))
        return(FALSE)
    if (any(.msp_lengths(sp) != n_peaks))
        return(TRUE)
    .pa_spectra_has_queue(sp) &&
        !identical(.msp_peaks(sp), .msp_peaks(Spectra::reset(sp)))
}

# Run metadata of averages computed with AverageFragSpectraParam.
.MZS_AV_CONFIGURED_TYPES <- c(
    collision_energy = "f64", x_mspurity_activation_type = "str",
    x_mspurity_parameter_set = "str", x_mspurity_average_group = "str",
    x_mspurity_spectrum_key = "str", x_mspurity_processing_step = "str",
    x_mspurity_acquisition_configuration = "str",
    x_mspurity_acquisition_values = "str",
    x_mspurity_acquisition_status = "str",
    x_mspurity_precursor_path = "str",
    x_mspurity_ancestor_conditions = "str")

#' The averages computed with AverageFragSpectraParam, in the columns of
#' `.mzs_derived_index()`, and their peaks in the columns of s_peaks.
#'
#' Each is identified as `msPurity:av_<method>:<grpid>:<set>:<ordinal>`,
#' which no average computed without param takes. Its pid is a negative
#' number, a join key no scan or stored spectrum has.
#'
#' @return `NULL`, or list(av, peaks).
#'
#' @noRd
.mzs_live_averages <- function(pa, files, uid, metadata) {
    sp <- pa@avSpectra
    configured <- .pa_av_configured(sp)
    if (!any(configured))
        return(NULL)
    view <- averagedSpectra(pa, legacy = FALSE)
    view <- view[.pa_av_configured(view)]
    n <- length(view)
    vars <- Spectra::spectraVariables(view)
    col <- function(v, na = NA) if (v %in% vars) view[[v]] else rep(na, n)
    json <- function(v) vapply(col(v, list(NULL)), function(x)
        if (is.null(x)) NA_character_
        else as.character(.mzs_to_json(.mzs_json_safe(x), pretty = FALSE)), "")
    method <- sub("^av_", "", as.character(view$av_level))
    if (any(!method %in% c("intra", "inter", "all")))
        .mzs_abort("semantic", "An averaged spectrum has an unknown level.")
    fileid <- ifelse(method == "intra", suppressWarnings(
        as.integer(as.character(col("sample")))), NA_integer_)
    k <- match(fileid, files$fileid)
    intra <- method == "intra"
    if (any(intra & is.na(k)))
        .mzs_abort("semantic", "Averaged spectra name file(s) absent from ",
                   "fileinfo.")
    set <- as.character(col("msPurityParameterSet", NA_character_))
    group <- as.character(col("msPurityAverageGroup", NA_character_))
    ordinal <- suppressWarnings(as.integer(sub("^.*/", "", group)))
    grpid <- as.integer(as.character(view$grpid))
    ## The MS2 precursor of each average: its own, or for MS3/MS4 averages
    ## the first stage of the precursor path.
    paths <- col("msPurityPrecursorPath", list(NULL))
    link_mz <- vapply(seq_len(n), function(i) {
        p <- paths[[i]]
        if (is.data.frame(p) && nrow(p) && any(p$msLevel %in% 2L))
            stats::median(p$precursorMz[p$msLevel %in% 2L])
        else as.numeric(Spectra::precursorMz(view)[i])
    }, numeric(1))
    md <- function(name) if (is.list(metadata) && !is.null(metadata[[name]]))
        rep(as.character(metadata[[name]]), n) else rep(NA_character_, n)
    av <- data.frame(
        pid = -seq_len(n), fileid = fileid, spectrum_type = method,
        precursor_mz = as.numeric(Spectra::precursorMz(view)),
        retention_time = as.numeric(Spectra::rtime(view)),
        inPurity = NA_real_, grpid = grpid, polarity = md("polarity"),
        instrument = md("instrument"), instrument_type = md("instrument_type"),
        method = method,
        run_id = ifelse(intra, paste0("av_intra_", files$run_id[k]),
                        paste0("av_", method)),
        data_origin = ifelse(intra, files$data_origin[k],
                             paste0("mspurity:", uid, ":av_", method)),
        source_run = ifelse(intra, files$run_id[k], NA_character_),
        source_key = ifelse(intra, files$source[k], NA_character_),
        native_id = paste0("msPurity:av_", method, ":", grpid, ":",
                           utils::URLencode(set, reserved = TRUE), ":",
                           ordinal),
        ms_level = as.integer(Spectra::msLevel(view)),
        configured = TRUE, is_null = as.logical(col("is_null", FALSE)),
        link_precursor_mz = link_mz,
        parameter_set = set, average_group = group, ordinal = ordinal,
        spectrum_key = as.character(view$spectrumId),
        processing_step = as.character(col("msPurityProcessingId",
                                           NA_character_)),
        collision_energy = as.numeric(col("collisionEnergy", NA_real_)),
        activation_type = as.character(col("activationType", NA_character_)),
        acquisition_configuration = as.character(
            col("msPurityAcquisitionConfiguration", NA_character_)),
        acquisition_values = json("msPurityAcquisitionValues"),
        acquisition_status = json("msPurityAcquisitionStatus"),
        precursor_path = json("msPurityPrecursorPath"),
        ancestor_conditions = json("msPurityAncestorConditions"),
        stringsAsFactors = FALSE)
    if (anyDuplicated(paste(av$run_id, av$native_id)))
        .mzs_abort("semantic", "Averages computed with AverageFragSpectraParam ",
                   "do not have distinct parameter-set groups.")
    pv <- Spectra::peaksVariables(view)
    pk <- as.list(Spectra::peaksData(view, columns = pv))
    counts <- vapply(pk, nrow, integer(1L))
    peaks <- do.call(rbind, lapply(pk, as.data.frame))
    if (is.null(peaks))
        peaks <- as.data.frame(matrix(numeric(), 0L, length(pv),
                                      dimnames = list(NULL, pv)))
    names(peaks)[names(peaks) == "intensity"] <- "i"
    peaks$pid <- rep(av$pid, counts)
    peaks$type <- rep(av$spectrum_type, counts)
    rownames(peaks) <- NULL
    list(av = av, peaks = peaks)
}

#' Averages from SQLite's tables and from the live object together, in write
#' order, with dataset spectrum ids. Averages computed without param keep
#' the order and ids they have alone.
#'
#' @noRd
.mzs_av_combine <- function(av, live, files) {
    av$ms_level <- rep(2L, nrow(av))
    av$configured <- rep(FALSE, nrow(av))
    av$link_precursor_mz <- av$precursor_mz
    if (is.null(live))
        return(av)
    all <- .mzs_rbind_fill(av, live$av)
    ord <- order(match(all$method, c("intra", "inter", "all")),
                 match(all$source_run, files$run_id),
                 all$grpid, all$ms_level, all$configured,
                 all$parameter_set, all$ordinal, all$pid)
    all <- all[ord, , drop = FALSE]
    all$spectrum_id_ <- seq_len(nrow(all))
    rownames(all) <- NULL
    all
}

#' Each exported average's node in the provenance graph, whether its
#' recorded provenance can be used (complete, peaks unchanged), and the
#' exported rank of each of its stored peaks.
#'
#' Averages computed without param are found by feature, level and file.
#'
#' @param pk the exported peaks (pid, mz, cl), as `.mzs_derived_runs()`
#'     orders them.
#'
#' @noRd
.mzs_av_graph <- function(pa, av, pk) {
    n <- nrow(av)
    none <- list(graph = NULL, node = rep(NA_integer_, n),
                 use = rep(FALSE, n), peak_rank = vector("list", n))
    g <- if (!is.null(pa)) .pa_graph_data(pa)
    sp <- if (!is.null(pa)) pa@avSpectra
    if (is.null(g) || !n || !length(sp))
        return(none)
    keys <- .pa_graph_sp_keys(sp, g)
    configured <- .pa_av_configured(sp)
    sample <- if ("sample" %in% Spectra::spectraVariables(sp))
        as.character(sp$sample) else rep(NA_character_, length(sp))
    level <- as.character(sp$av_level)
    legacy_key <- paste(as.character(sp$grpid), level,
                        ifelse(level == "av_intra", sample, NA))
    live <- av$configured %in% TRUE
    row <- rep(NA_integer_, n)
    row[live] <- match(av$spectrum_key[live], .pa_graph_sp_ids(sp))
    want <- paste(as.character(av$grpid), paste0("av_", av$spectrum_type),
                  ifelse(av$spectrum_type == "intra",
                         as.character(av$fileid), NA))
    cand <- which(!configured)
    row[!live] <- cand[match(want[!live], legacy_key[cand])]
    node <- match(keys[row], g$nodes$key)
    counts <- tabulate(match(pk$pid, av$pid), n)
    use <- !is.na(node) & g$nodes$complete[node] %in% TRUE &
        g$nodes$nPeaks[node] %in% counts
    use[use] <- vapply(which(use), function(i)
        !.mzs_av_peaks_changed(sp[row[i]], g$nodes$nPeaks[node[i]]),
        logical(1))
    peak_rank <- vector("list", n)
    stored <- if (any(use)) as.list(Spectra::peaksData(
        sp[row[use]], columns = intersect(c("mz", "cl"),
                                          Spectra::peaksVariables(sp))))
    for (j in seq_along(stored)) {
        m <- stored[[j]]
        o <- if ("cl" %in% colnames(m)) order(m[, "mz"], m[, "cl"])
             else order(m[, "mz"])
        r <- integer(length(o))
        r[o] <- seq_along(o)
        peak_rank[[which(use)[j]]] <- r
    }
    list(graph = g, node = node, use = use, peak_rank = peak_rank)
}

#' merge_member rows of the averages whose provenance is recorded: every
#' selected immediate member, whether it contributed a stored peak, and its
#' feature link. Measured members are referenced as scans, earlier averages
#' as spectra of the dataset.
#'
#' @noRd
.mzs_graph_merge_member <- function(av, gk, scans, scan_refs) {
    g <- gk$graph
    rows <- lapply(which(gk$use), function(i) {
        m <- .pa_graph_step(g, g$nodes$key[gk$node[i]])
        if (!nrow(m))
            return(NULL)
        measured <- !is.na(m$linkId)
        li <- match(m$linkId, g$links$linkId)
        k <- match(g$links$pid[li], scans$pid)
        d <- match(m$member, g$nodes$key[gk$node])
        if (any(measured & is.na(k)) || any(!measured & is.na(d)))
            .mzs_abort("semantic", "The recorded inputs of averaged spectrum ",
                       av$native_id[i], " are not among the exported scans ",
                       "and averages.")
        ref <- function(scan_col, av_col) ifelse(measured, scan_refs[[scan_col]][k],
                                                 av[[av_col]][d])
        data.frame(
            merged_spectrum_id_ = av$spectrum_id_[i],
            merged_run_id = av$run_id[i], merged_native_id = av$native_id[i],
            member_source = ifelse(measured, scan_refs$source[k], "self"),
            member_run_id = ref("run_id", "run_id"),
            member_spectrum_id_ = as.numeric(ref("spectrum_id_",
                                                 "spectrum_id_")),
            member_native_id = ref("native_id", "native_id"),
            x_mspurity_member_scan_number = ifelse(measured,
                scans$acquisitionNum[k], NA_integer_),
            member_rank = as.integer(m$rank), members_complete = TRUE,
            x_mspurity_contributes = as.logical(m$contributes),
            x_mspurity_feature_link_index = as.integer(
                g$links$featureLinkIndex[li]),
            x_mspurity_contributor_index = as.integer(
                g$links$contributorIndex[li]),
            stringsAsFactors = FALSE)
    })
    out <- do.call(.mzs_rbind_fill, rows)
    if (!is.null(out))
        out$merged_source <- "self"
    out
}

#' The x_mspurity_merge_peak table: the stored peak of each immediate input
#' every averaged peak was made from, for averages whose provenance is
#' recorded. Peaks of averages are given by their exported rank; peaks of
#' measured scans by their position on the recorded basis.
#'
#' @noRd
.mzs_merge_peak <- function(av, gk) {
    g <- gk$graph
    rows <- lapply(which(gk$use), function(i) {
        e <- .pa_graph_step(g, g$nodes$key[gk$node[i]], peaks = TRUE)
        if (!nrow(e))
            return(NULL)
        d <- match(e$member, g$nodes$key[gk$node])
        derived <- is.na(e$linkId)
        member_peak <- as.integer(e$memberPeak)
        for (j in which(derived)) {
            r <- gk$peak_rank[[d[j]]]
            member_peak[j] <- if (is.null(r) || is.na(member_peak[j])) NA
                              else r[member_peak[j]]
        }
        out <- data.frame(
            merged_spectrum_id_ = av$spectrum_id_[i],
            merged_peak_rank = gk$peak_rank[[i]][e$mergedPeak],
            member_rank = as.integer(e$rank),
            member_peak_rank = member_peak,
            basis = as.character(e$basis), stringsAsFactors = FALSE)
        o <- order(out$merged_peak_rank, out$member_rank,
                   out$member_peak_rank, seq_len(nrow(out)))
        out <- out[o, , drop = FALSE]
        out$edge_rank <- seq_len(nrow(out))
        out
    })
    out <- do.call(.mzs_rbind_fill, rows)
    if (!is.null(out)) rownames(out) <- NULL
    out
}

#' Loss ledger rows for the averaged spectra.
#'
#' @noRd
.mzs_av_losses <- function(av, gk, merge_peak, pa) {
    if (is.null(pa) || !nrow(av))
        return(list())
    configured <- av$configured %in% TRUE
    used <- if (is.null(gk)) rep(FALSE, nrow(av)) else gk$use
    list(
        if (any(used))
            list(construct = paste("avSpectra: msPurityContributors,",
                                   "msPurityClusterContributors"),
                 disposition = "dropped_by_configuration",
                 affected_rows = sum(used), target_ref = "merge_member",
                 reason = paste("The inputs and peaks each average was made",
                                "from, recomputable from merge_member and",
                                "x_mspurity_merge_peak.")),
        if (any(configured))
            list(construct = "avSpectra: precursor purity",
                 disposition = "no_such_concept_in_source",
                 affected_rows = sum(configured),
                 reason = paste("Averages computed with AverageFragSpectraParam",
                                "have no spectrum-level purity; each peak",
                                "keeps x_mspurity_in_purity.")),
        if (!is.null(gk$graph) && any(!used & !configured))
            list(construct = "averageGraph",
                 disposition = "not_present_in_source",
                 affected_rows = sum(!used & !configured),
                 target_ref = "merge_member",
                 reason = paste("No complete record of these averages'",
                                "inputs, or their peaks changed after",
                                "averaging; their members are",
                                "reconstructed.")),
        if (!is.null(merge_peak) && any(merge_peak$basis == "unknown"))
            list(construct = "averageGraph: member peak positions",
                 disposition = "not_present_in_source",
                 affected_rows = sum(merge_peak$basis == "unknown"),
                 target_ref = "x_mspurity_merge_peak.member_peak_rank",
                 reason = paste("The provenance does not record which peak",
                                "of these inputs a peak came from.")))
}
