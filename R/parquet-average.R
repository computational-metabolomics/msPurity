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
.mzs_live_averages <- function(pa, files, uid, metadata, rows = NULL) {
    sp <- pa@avSpectra
    configured <- .pa_av_configured(sp)
    if (!any(configured))
        return(NULL)
    if (is.null(rows)) rows <- which(configured)
    view <- .pa_graph_view(pa, sp[rows])
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
        pid = -match(rows, which(configured)), fileid = fileid, spectrum_type = method,
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
#' @param counts the number of retained peaks in each exported spectrum.
#'
#' @noRd
.mzs_av_graph <- function(pa, av, counts, dir) {
    n <- nrow(av)
    none <- list(graph = NULL, node = rep(NA_integer_, n),
                 use = rep(FALSE, n))
    g <- .pa_graph_data(pa)
    sp <- pa@avSpectra
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
    use <- !is.na(node) & g$nodes$complete[node] %in% TRUE &
        g$nodes$nPeaks[node] == counts
    if (.pa_spectra_has_queue(sp))
        use[use] <- vapply(which(use), function(i)
            !.mzs_av_peaks_changed(sp[row[i]], counts[i]), logical(1))
    rank_columns <- intersect(c("mz", "cl"), Spectra::peaksVariables(sp))
    ranks <- function(i) {
        pk <- as.list(Spectra::peaksData(sp[row[i]], columns = rank_columns,
                                        BPPARAM = .msp_bpparam()))
        lapply(pk, function(m) {
            o <- if ("cl" %in% colnames(m)) order(m[, "mz"], m[, "cl"]) else order(m[, "mz"])
            order(o)
        })
    }
    # Inter-file edges may refer to any earlier intra-file average. Save
    # ranks in bounded blocks and read back just the parent block required.
    size <- .mzs_batch_size(TRUE)
    available <- which(use)
    position <- match(seq_len(n), available)
    .mzs_batches(length(available), size, function(i)
        saveRDS(ranks(available[i]), file.path(dir, paste0("rank-", i[1L])), compress = FALSE))
    last <- 0L
    cache <- NULL
    rank <- function(i) {
        if (!length(i) || is.na(i) || !use[i]) return(NULL)
        pos <- position[i]
        first <- as.integer((pos - 1L) %/% size * size + 1L)
        if (first != last) {
            cache <<- readRDS(file.path(dir, paste0("rank-", first)))
            last <<- first
        }
        cache[[pos - first + 1L]]
    }
    list(graph = g, node = node, use = use, rank = rank)
}

#' merge_member rows of the averages whose provenance is recorded: every
#' selected immediate member, whether it contributed a stored peak, and its
#' feature link. Measured members are referenced as scans, earlier averages
#' as spectra of the dataset.
#'
#' @noRd
.mzs_graph_merge_member <- function(av, gk, scans, scan_refs, rows = which(gk$use)) {
    g <- gk$graph
    rows <- lapply(rows, function(i) {
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
.mzs_merge_peak <- function(av, gk, rows = which(gk$use)) {
    g <- gk$graph
    rows <- lapply(rows, function(i) {
        e <- .pa_graph_step(g, g$nodes$key[gk$node[i]], peaks = TRUE)
        if (!nrow(e))
            return(NULL)
        d <- match(e$member, g$nodes$key[gk$node])
        derived <- is.na(e$linkId)
        member_peak <- as.integer(e$memberPeak)
        for (parent in unique(d[derived])) {
            j <- which(derived & d == parent)
            r <- gk$rank(parent)
            member_peak[j] <- if (is.null(r)) NA_integer_ else r[member_peak[j]]
        }
        out <- data.frame(
            merged_spectrum_id_ = av$spectrum_id_[i],
            merged_peak_rank = gk$rank(i)[e$mergedPeak],
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
    unknown <- if (inherits(merge_peak, "mzs_table_batches")) merge_peak$unknown else
        sum(merge_peak$basis == "unknown")
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
        if (unknown > 0)
            list(construct = "averageGraph: member peak positions",
                 disposition = "not_present_in_source",
                 affected_rows = unknown,
                 target_ref = "x_mspurity_merge_peak.member_peak_rank",
                 reason = paste("The provenance does not record which peak",
                                "of these inputs a peak came from.")))
}

# The scalar writer assigns positive pids to non-empty averages in the
# legacy group's iteration order. Keep their identities when adding empty
# spectra, which have no corresponding SQLite pid.
.mzs_scalar_index <- function(pa, groups, scaninfo, reserved) {
    sp <- pa@avSpectra
    if (!length(sp)) return(NULL)
    i <- which(!.pa_av_configured(sp) & !sp$is_null)
    if (!length(i)) return(NULL)
    group <- as.integer(as.character(sp$grpid[i]))
    at <- match(group, groups$grpid)
    if (anyNA(at)) .mzs_abort("semantic", "An averaged spectrum names an unknown feature.")
    method <- sub("^av_", "", sp$av_level[i])
    file <- as.integer(as.character(sp$sample[i]))
    purity <- total <- rep(NA_real_, length(i))
    size <- integer(length(i))
    pass <- 0
    present <- character()
    pv <- Spectra::peaksVariables(sp)
    .mzs_batches(length(i), .mzs_batch_size(TRUE), function(j) {
        pk <- as.list(Spectra::peaksData(sp[i[j]], columns = pv,
                                        BPPARAM = .msp_bpparam()))
        size[j] <<- vapply(pk, nrow, integer(1))
        if ("inPurity" %in% pv)
            purity[j] <<- vapply(pk, function(m) m[, "inPurity"][1L], 0)
        if ("total" %in% pv)
            total[j] <<- vapply(pk, function(m) {
                x <- m[, "total"]
                if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
            }, 0)
        if ("pass_flag" %in% pv)
            pass <<- pass + sum(vapply(pk, function(m) sum(!is.na(m[, "pass_flag"])), 0))
        present <<- union(present, pv[vapply(pv, function(c)
            any(vapply(pk, function(m) any(!is.na(m[, c])), logical(1))), logical(1))])
    })
    pid <- numeric(length(i))
    pid[size > 0] <- max(c(0, scaninfo$pid)) + seq_len(sum(size > 0))
    # Internal join keys for empty averages occupy a separate range from
    # configured averages (whose keys are -1 .. -length(configured)).
    pid[size == 0] <- -(reserved + seq_len(sum(size == 0)))
    list(meta = data.frame(pid = pid, fileid = file, spectrum_type = method,
        precursor_mz = groups$mz[at], retention_time = groups$rt[at],
        inPurity = purity, grpid = group), row = i, size = size,
        total = total, present = present, pass = pass, columns = pv)
}

.mzs_absent_ledger <- function(pa, base, dir) {
    parts <- .mzs_spool(dir, "loss_ledger")
    parts$add(base)
    count <- nrow(base)
    sp <- pa@avSpectra
    groups <- .pa_spectra_metadata(sp)$groups
    for (level in c("intra", "inter", "all")) {
        absent <- pa@params$absentAverages[[level]]
        if (is.null(absent)) {
            settings <- methods::slot(pa, paste0("av_", level, "_params"))
            if (!length(settings) || !is.null(settings$parameter_sets)) next
            # Older objects have no attempt log. Each recorded scalar call
            # visited all groups in its legacy list.
            at <- if (length(sp)) which(!.pa_av_configured(sp) &
                sp$av_level == paste0("av_", level)) else integer()
            missing <- groups
            absent <- data.frame(group = character(), sample = character())
            if (length(at)) {
                null <- at[sp$is_null[at] %in% TRUE]
                missing <- setdiff(groups, as.character(sp$grpid[at]))
                absent <- data.frame(group = as.character(sp$grpid[null]),
                                     sample = as.character(sp$sample[null]))
            }
            if (length(missing))
                absent <- rbind(absent, data.frame(group = missing, sample = NA_character_))
        }
        .mzs_batches(nrow(absent), .mzs_batch_size(), function(i) {
            ref <- paste0("av_", level, ": feature ", absent$group[i])
            sample <- absent$sample[i]
            known <- !is.na(sample)
            ref[known] <- paste0(ref[known], ", sample ", sample[known])
            parts$add(data.frame(loss_ledger_id_ = count + seq_along(i),
                source_construct = ref, disposition = "not_present_in_source",
                affected_rows = 1, target_ref = NA_character_,
                reason = "Averaging was attempted and produced no spectrum."))
            count <<- count + length(i)
        })
    }
    parts$finish()
}

# A run reads peak arrays only when its next batch is written.
.mzs_peak_reader <- function(sp, rows, columns) {
    force(rows)
    force(columns)
    function(i) {
        mats <- as.list(Spectra::peaksData(sp[rows[i]],
            columns = c("mz", "intensity", columns), BPPARAM = .msp_bpparam()))
        mats <- lapply(mats, function(m) {
            o <- if ("cl" %in% colnames(m)) order(m[, "mz"], m[, "cl"])
                 else order(m[, "mz"])
            m[o, , drop = FALSE]
        })
        peak <- function(c, type, optional) .mzs_peak_list(lapply(mats, function(m) {
            x <- m[, c]
            if (optional && all(is.na(x))) NULL else .mzs_as_r_type(unname(x), type)
        }), type)
        out <- list(mz = peak("mz", "f64", FALSE), intensity = peak("intensity", "f64", FALSE))
        for (c in columns) {
            spec <- .MZS_AV_PEAKS[[c]]
            out[[spec[1L]]] <- peak(c, spec[2L], TRUE)
        }
        out
    }
}

.mzs_live_derived <- function(pa, src, av, files, uid, metadata, dir) {
    sp <- pa@avSpectra
    configured <- which(.pa_av_configured(sp))
    present <- src$scalar$present
    size <- src$scalar$size
    rows <- which(!.pa_av_configured(sp))[src$scalar$row]
    pids <- src$scalar$meta$pid
    chunks <- list()
    .mzs_batches(length(configured), .mzs_batch_size(TRUE), function(i) {
        at <- configured[i]
        x <- .mzs_live_averages(pa, files, uid, metadata, at)
        chunks[[length(chunks) + 1L]] <<- x$av
        nonnull <- vapply(x$peaks, function(v) any(!is.na(v)), logical(1))
        present <<- union(present, names(x$peaks)[nonnull])
        size <<- c(size, tabulate(match(x$peaks$pid, x$av$pid), length(at)))
    })
    live <- if (length(chunks)) list(av = do.call(rbind, chunks))
    av <- .mzs_av_combine(av, live, files)
    k <- match(av$pid, c(pids, -seq_along(configured)))
    rows <- c(rows, configured)[k]
    if (anyNA(rows)) .mzs_abort("semantic", "An average has no stored spectrum.")
    derived <- .mzs_derived_runs(av, NULL)
    columns <- intersect(names(.MZS_AV_PEAKS), present)
    unmapped <- setdiff(Spectra::peaksVariables(sp), c("mz", "intensity",
        names(.MZS_AV_PEAKS), "frac", "pass_flag", "purity_pass_flag",
        "intensity_pass_flag"))
    if (length(rows) && length(unmapped))
        .mzs_abort("unsupported", "Unmapped averaged peak columns: ", paste(unmapped, collapse = ", "))
    for (j in seq_along(derived$runs)) {
        id <- derived$runs[[j]]$run_id
        derived$runs[[j]]$peaks <- NULL
        derived$runs[[j]]$read_peaks <- .mzs_peak_reader(sp, rows[av$run_id == id], columns)
        derived$peak_columns[[id]] <- unname(vapply(.MZS_AV_PEAKS[columns], `[[`, "", 1L))
    }
    av$empty <- size[k] == 0L
    empty <- which(av$empty & !av$configured)
    if (length(empty)) {
        if (is.null(av$spectrum_key)) av$spectrum_key <- NA_character_
        av$spectrum_key[empty] <- .pa_graph_sp_ids(sp)[rows[empty]]
    }
    list(av = av, derived = derived, graph = .mzs_av_graph(pa, av, size[k], dir))
}
