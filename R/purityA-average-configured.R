# Averaging with AverageFragSpectraParam (all-spectrum, intra- and
# inter-file). Inputs, the numerical worker and contributor records are those
# of the scalar route; this file selects feature links, groups compatible ones
# and describes the results.

# TRUE for averaged spectra computed with AverageFragSpectraParam.
.pa_av_configured <- function(sp) {
    vars <- Spectra::spectraVariables(sp)
    configured <- rep(FALSE, length(sp))
    for (field in intersect(c("msPurityParameterSet", "msPurityAverageGroup"), vars))
        configured <- configured | !is.na(sp[[field]])
    configured
}

.PA_AV_CONFIGURED_LEGACY <- paste(
    "Averages computed with AverageFragSpectraParam are available as Spectra",
    "only: use averagedSpectra(pa, legacy = FALSE). The legacy list holds one",
    "all-spectrum average per feature.")

# Configured-average variables, including older columns for scalar cleanup.
.PA_AV_CONFIGURED_COLUMNS <- c("msPurityParameterSet",
    "msPurityAverageGroup", "msPurityAverageStatus",
    "msPuritySelectedContributors", "msPurityPrecursorPath",
    "msPurityAncestorConditions", "msPurityAverageInputs")

# Averaging without param replaces configured averages of its own level (and
# inter-file averages, for intra), including settings recorded by a
# configured call that stored no averages; the legacy representation cannot
# keep configured averages of other levels.
.pa_average_scalar_ready <- function(pa, av_level) {
    sp <- pa@avSpectra
    configured <- .pa_av_configured(sp)
    replaced <- switch(av_level, all = "all", intra = c("intra", "inter"),
        inter = "inter")
    recorded <- vapply(replaced, function(level) !is.null(methods::slot(pa,
        paste0("av_", level, "_params"))$parameter_sets), logical(1))
    if (!any(configured) && !any(recorded))
        return(pa)
    remaining <- if (any(configured)) setdiff(unique(sp$av_level[configured]),
        paste0("av_", replaced)) else character()
    remaining <- paste(remaining, collapse = ", ")
    if (nzchar(remaining))
        stop("Averaging without param cannot keep the averages computed ",
            "with AverageFragSpectraParam at other levels (", remaining,
            "). Repeat those averages ",
            "without param first (possible for MS2-only feature links), or ",
            "use param for this level too.")
    pa <- .pa_processing_invalidate(pa, replaced)
    sp <- pa@avSpectra
    typed <- intersect(.PA_AV_CONFIGURED_COLUMNS,
        Spectra::spectraVariables(sp))
    if (length(typed))
        pa@avSpectra <- Spectra::selectSpectraVariables(sp,
            setdiff(Spectra::spectraVariables(sp), typed))
    pa
}

# Integer codes for exact comparison of acquisition conditions. Activation
# names and accessions are normalized as for select; NA matches only NA.
.pa_average_condition_codes <- function(metadata, axes) {
    lapply(stats::setNames(axes, axes), function(axis) {
        x <- metadata[[axis]]
        if (axis == "activationType") x <- .pa_activation_type(x)
        match(x, unique(x))
    })
}

# One row per feature link (grped_df row) with its position in the feature,
# source spectrum, acquisition metadata and precursor path.
.pa_average_links <- function(pa) {
    g <- pa@grped_df
    frag <- pa@fragSpectra
    source <- match(g$pid, frag$pid)
    if (anyNA(source))
        stop("Feature links refer to unavailable fragmentation spectra.")
    grpid <- as.character(g$grpid)
    links <- data.frame(featureLinkIndex = seq_len(nrow(g)), grpid = grpid,
        contributorIndex = stats::ave(seq_len(nrow(g)), grpid,
            FUN = seq_along),
        pid = g$pid, sample = if ("sample" %in% names(g))
            as.character(g$sample) else rep(NA_character_, nrow(g)),
        source = source,
        stringsAsFactors = FALSE)
    links <- cbind(links, .pa_filter_metadata(frag)[source, , drop = FALSE])
    rownames(links) <- NULL
    # Numeric order where the values are numbers, text order otherwise.
    variable <- function(name) {
        if (name %in% Spectra::spectraVariables(frag))
            utils::type.convert(as.character(frag[[name]][source]),
                as.is = TRUE)
        else rep(NA_real_, length(source))
    }
    links$fileid <- as.character(variable("fileid"))
    canonical <- order(variable("fileid"), variable("seqNum"),
        variable("acquisitionNum"), links$featureLinkIndex, na.last = TRUE,
        method = "radix")
    links$.rank <- order(canonical)
    paths <- .pa_average_paths(pa, links)
    links$.stages <- paths$stages
    links$.ancestors <- paths$ancestors
    attr(links, "ancestry") <- paths$scans
    links
}

# Precursor m/z of each stage after MS2, and a key for the conditions of the
# earlier stages, from the validated ancestry stored by purityA().
.pa_average_paths <- function(pa, links) {
    n <- nrow(links)
    stages <- rep(list(numeric()), n)
    ancestors <- rep("", n)
    higher <- which(links$msLevel > 2L)
    if (!length(higher))
        return(list(stages = stages, ancestors = ancestors, scans = NULL))
    .pa_msn_validate_preparation(pa)
    lineage <- pa@params$msnLineage
    if (is.null(lineage))
        stop("Averaging MS3/MS4 feature links requires the precursor paths ",
            "stored by purityA() and frag4feature() with msLevel including ",
            "3 or 4.")
    scans <- .pa_msn_stage_metadata(lineage$scans)
    rows <- match(links$pid[higher], scans$pid)
    if (anyNA(rows) || any(scans$pathStatus[rows] != "valid"))
        stop("MS3/MS4 feature links require validated precursor paths.")
    codes <- .pa_average_condition_codes(scans, names(.PA_FILTER_AXES))
    conditions <- do.call(paste, c(unname(codes), sep = "."))
    for (k in seq_along(higher)) {
        path <- match(scans$path[[rows[k]]], scans$pid)
        stages[[higher[k]]] <- scans$precursor_mz[path[-1L]]
        ancestors[higher[k]] <- paste(conditions[path[-length(path)]],
            collapse = "/")
    }
    list(stages = stages, ancestors = ancestors, scans = scans)
}

# Apply the fragment-clustering cutoff at every precursor stage after MS2.
.pa_average_compatible <- function(a, b, ppm, tolerance) {
    if (!length(a)) return(TRUE)
    limit <- .hc_distance(abs(a - b), (a + b) / 2, ppm, tolerance)
    all(limit$distance <= limit$height)
}

# Every link of a agrees with every link of b at each stage after MS2.
.pa_average_units_compatible <- function(links, a, b, spec) {
    for (i in a) for (j in b)
        if (!.pa_average_compatible(links$.stages[[i]], links$.stages[[j]],
                spec$ppm, spec$tolerance)) return(FALSE)
    TRUE
}

# Groups of compatible units for one parameter set; a unit is a vector of
# link rows (one link, or the contributing links of an intra-file average).
# Units are visited in canonical order of their first link; each joins the
# first group whose every link is compatible with all of its links, so group
# membership does not depend on input order. Units sharing a value of
# exclusive (the file, for inter-file averaging) never share a group.
# Returns unit indices per group.
.pa_average_unit_groups <- function(links, units, spec, exclusive = NULL) {
    axes <- setdiff(names(.PA_FILTER_AXES), spec$combine)
    codes <- .pa_average_condition_codes(links, axes)
    partition <- do.call(paste, c(list(links$grpid), unname(codes),
        list(links$.ancestors), sep = "\r"))
    visit <- order(vapply(units, function(unit) min(links$.rank[unit]),
        numeric(1)))
    groups <- list()
    for (u in visit) {
        unit <- units[[u]]
        fit <- 0L
        for (j in seq_along(groups)) {
            members <- unlist(units[groups[[j]]], use.names = FALSE)
            if (partition[members[1L]] != partition[unit[1L]]) next
            if (!is.null(exclusive) &&
                    exclusive[u] %in% exclusive[groups[[j]]]) next
            if (.pa_average_units_compatible(links, members, unit, spec)) {
                fit <- j
                break
            }
        }
        if (fit) groups[[fit]] <- c(groups[[fit]], u)
        else groups[[length(groups) + 1L]] <- u
    }
    groups
}

# Groups of compatible link rows, each link being its own unit.
.pa_average_groups <- function(links, rows, spec) {
    lapply(.pa_average_unit_groups(links, as.list(rows), spec),
        function(group) rows[group])
}

# Average one group with the scalar route's peak table and worker.
.pa_average_group <- function(pa, links, rows, spec, rmp, spectra, filtered,
    cores, blockSize) {
    grpid <- links$grpid[rows[1L]]
    feature <- which(as.character(pa@grped_df$grpid) == grpid)
    info <- pa@grped_df[feature, , drop = FALSE]
    info$index <- seq_len(nrow(info))
    index <- links$contributorIndex[rows]
    table <- .pa_average_peak_table(info[index, , drop = FALSE],
        spectra[[grpid]][index], filtered)
    table <- table[order(table$mz), ]
    result <- average_spectra(table, cores = cores, ppm = spec$ppm,
        minnum = spec$minnum, sumi = spec$sumi, minfrac = spec$minfrac,
        snthr = spec$snr, rathr = spec$ra, av_type = spec$av,
        record_inputs = TRUE, tolerance = spec$tolerance, blockSize = blockSize)
    peaks <- result$peaks
    if (rmp && !is.null(peaks))
        peaks <- peaks[peaks$pass_flag, ]
    list(peaks = peaks, effective_cores = result$effective_cores,
        clusters = .pa_average_used_members(result$members, peaks),
        graph = list(links = links$featureLinkIndex[rows],
            edges = .pa_average_peak_edges(result, peaks, ranks = index)))
}

# Median precursor m/z of each stage, and the earlier stages' conditions, of
# the feature links contributing to an average.
.pa_average_path_summary <- function(links, rows) {
    level <- links$msLevel[rows[1L]]
    if (!length(rows) || level < 3L)
        return(list(path = NULL, ancestors = NULL))
    mz <- do.call(rbind, links$.stages[rows])
    path <- data.frame(msLevel = seq.int(3L, level),
        precursorMz = apply(mz, 2L, stats::median))
    # Ancestor conditions are identical within a group by construction.
    scans <- attr(links, "ancestry")
    first <- scans$path[[match(links$pid[rows[1L]], scans$pid)]]
    first <- match(first[-length(first)], scans$pid)
    ancestors <- scans[first, names(.PA_FILTER_AXES), drop = FALSE]
    ancestors$activationType <- .pa_activation_type(ancestors$activationType)
    rownames(ancestors) <- NULL
    list(path = path, ancestors = ancestors)
}

# Scalar metadata of one average, from the links supporting its stored peaks.
.pa_average_row <- function(pa, links, rows, set, ordinal, result,
    level = "av_all", sample = NA_character_) {
    grpid <- links$grpid[rows[1L]]
    used_index <- unique(unlist(result$clusters, use.names = FALSE))
    summary <- .pa_average_summary(pa, grpid, result$clusters,
        normalise = TRUE)
    used_rows <- rows[links$contributorIndex[rows] %in% used_index]
    feature <- pa@grped_df[as.character(pa@grped_df$grpid) == grpid, ,
        drop = FALSE]
    used_feature <- feature[sort(used_index), , drop = FALSE]
    precursor <- if (links$msLevel[rows[1L]] == 2L)
        used_feature$precurMtchMZ else
        Spectra::precursorMz(pa@fragSpectra)[links$source[used_rows]]
    group <- paste0(set, "/", ordinal)
    list(grpid = grpid, av_level = level, fileid = sample, sample = sample,
        msLevel = links$msLevel[rows[1L]],
        precursorMz = if (length(used_index))
            stats::median(precursor, na.rm = TRUE) else NA_real_,
        rtime = if (length(used_index))
            stats::median(used_feature$retentionTime, na.rm = TRUE) else
            NA_real_,
        msPurityParameterSet = set, msPurityAverageGroup = group,
        msPurityAcquisitionValues = summary$values,
        msPurityAcquisitionConfiguration = summary$configuration)
}

.pa_average_configured <- function(pa, level, param, supplied, rmp, cores, blockSize) {
    specs <- .pa_average_resolve(param, supplied, rmp)
    started <- Sys.time()
    step <- switch(level, all = "averageAllFragSpectra",
        intra = "averageIntraFragSpectra", inter = "averageInterFragSpectra")
    if (level == "inter") pa <- .pa_graph_import(pa)
    input <- .pa_processing_input(pa, step)
    if (level == "inter")
        return(.pa_average_inter(pa, specs, rmp, cores, blockSize, started, input))
    pa <- .pa_processing_invalidate(pa,
        if (level == "intra") c("intra", "inter") else "all")
    links <- if (nrow(pa@grped_df)) .pa_average_links(pa) else
        data.frame(featureLinkIndex = integer(), grpid = character())
    metadata <- links[, intersect(names(.PA_FILTER_AXES), names(links)),
        drop = FALSE]
    matched <- .pa_parameter_matches(metadata,
        lapply(specs, function(spec) list(select = spec$select)))
    spectra <- .pa_grouped_legacy(pa, peak_identity = TRUE)
    filtered <- length(pa@filter_frag_params) > 0
    rows <- list()
    peaks <- list()
    graph <- list()
    effective_cores <- 0L
    for (set in names(specs)) {
        selected <- which(matched[[set]])
        for (grpid in unique(links$grpid[selected])) {
            feature <- selected[links$grpid[selected] == grpid]
            # Intra-file averages partition by file before grouping.
            parts <- if (level == "intra") {
                ordered <- feature[order(links$.rank[feature])]
                if (anyNA(links$fileid[ordered]))
                    stop("Intra-file averaging requires source file identities.")
                split(ordered, factor(links$fileid[ordered],
                    levels = unique(links$fileid[ordered])))
            } else list(feature)
            for (part in parts) {
                sample <- if (level == "intra") links$fileid[part[1L]] else
                    NA_character_
                groups <- .pa_average_groups(links, part, specs[[set]])
                for (ordinal in seq_along(groups)) {
                    result <- .pa_average_group(pa, links, groups[[ordinal]],
                        specs[[set]], rmp, spectra, filtered, cores, blockSize)
                    effective_cores <- max(effective_cores, result$effective_cores)
                    row <- .pa_average_row(pa, links,
                        groups[[ordinal]], set, ordinal, result,
                        level = paste0("av_", level), sample = sample)
                    rows[[length(rows) + 1L]] <- row
                    graph[[length(rows)]] <- result$graph
                    peaks[length(peaks) + 1L] <- list(result$peaks)
                }
            }
        }
    }
    pa <- .pa_average_store(pa, rows, peaks, graph)
    assignments <- links[, intersect(c("featureLinkIndex", "grpid", "pid",
            "sample", "fileid", "msLevel"), names(links)), drop = FALSE]
    assignments$parameter_set_names <- if (nrow(links)) lapply(
        seq_len(nrow(links)), function(i) names(specs)[vapply(matched,
            `[`, logical(1), i)]) else list()
    methods::slot(pa, paste0("av_", level, "_params")) <- list(rmp = rmp,
        cores = cores, blockSize = blockSize, parameter_sets = specs)
    parameters <- list(rmp = rmp, cores = cores, effective_cores = effective_cores,
        blockSize = blockSize, configuration = list(parameter_sets = specs,
            assignments = assignments,
            unused_parameter_sets = names(specs)[!vapply(matched, any,
                    logical(1))]))
    .pa_processing_record(pa, step, parameters, started, input)
}

# Selectors compared by meaning: sorted fields and values, numeric levels,
# energies and polarities, and activation names mapped as for select.
.pa_average_select_key <- function(select) {
    fields <- sort(names(select))
    lapply(stats::setNames(fields, fields), function(axis) {
        value <- select[[axis]]
        value <- if (axis %in% c("msLevel", "collisionEnergy", "polarity"))
            as.numeric(value) else if (axis == "activationType")
            .pa_activation_type(value) else as.character(value)
        sort(unique(value), na.last = TRUE)
    })
}

# Inter-file averaging of intra-file averages computed with param. Each
# inter-file set uses the select and combine of the intra-file set of the
# same name; its own thresholds and tolerances apply.
.pa_average_inter <- function(pa, specs, rmp, cores, blockSize, started, input) {
    sp <- pa@avSpectra
    configured <- .pa_av_configured(sp)
    if (!any(configured) || !any(configured & sp$av_level == "av_intra"))
        stop("Inter-file averaging with param requires intra-file averages ",
            "computed with param.")
    intra_sets <- pa@av_intra_params$parameter_sets
    for (set in names(specs)) {
        reference <- intra_sets[[set]]
        if (is.null(reference))
            stop("Inter-file parameter set ", set, " has no intra-file ",
                "parameter set of the same name.")
        if (!identical(.pa_average_select_key(reference$select),
                    .pa_average_select_key(specs[[set]]$select)) ||
                !setequal(reference$combine, specs[[set]]$combine))
            stop("Inter-file parameter set ", set, " must use the select ",
                "and combine of its intra-file parameter set.")
    }
    pa <- .pa_processing_invalidate(pa, "inter")
    sp <- pa@avSpectra
    intra <- which(.pa_av_configured(sp) & sp$av_level == "av_intra" &
        sp$msPurityParameterSet %in% names(specs))
    links <- .pa_average_links(pa)
    peak_list <- .msp_peak_list(sp)
    intra <- intra[vapply(peak_list[intra], function(x)
        !is.null(x) && nrow(x) > 0L, logical(1))]
    columns <- .pa_spectra_metadata(sp)$columns
    graph <- pa@params$averageGraph
    exact <- graph$nodes$exactLeaves[match(sp$spectrumId, graph$nodes$spectrumId)]
    prepared <- lapply(intra, function(row)
        .pa_average_inter_input(pa, sp, row, peak_list, columns, links, exact[row]))
    names(prepared) <- as.character(intra)
    intra <- intra[vapply(prepared, function(x) nrow(x$peaks) > 0L,
        logical(1))]
    rows <- list()
    peaks <- list()
    graph <- list()
    effective_cores <- 0L
    for (set in names(specs)) {
        spec <- specs[[set]]
        candidates <- intra[sp$msPurityParameterSet[intra] == set]
        for (grpid in unique(sp$grpid[candidates])) {
            members <- candidates[sp$grpid[candidates] == grpid]
            units <- lapply(members, function(row)
                prepared[[as.character(row)]]$link_rows)
            for (k in seq_along(units))
                if (!.pa_average_units_compatible(links, units[[k]],
                        units[[k]], spec))
                    stop("The precursors of intra-file average ",
                        sp$spectrumId[members[k]], " differ by more than ",
                        "the inter-file tolerance of parameter set ", set,
                        ". Use at least the intra-file ppm and tolerance.")
            # Each file contributes at most one input to an average.
            files <- vapply(units, function(unit) {
                file <- unique(links$fileid[unit])
                if (length(file) != 1L || is.na(file))
                    stop("Each intra-file input requires one source file.")
                file
            }, character(1))
            groups <- .pa_average_unit_groups(links, units, spec,
                exclusive = files)
            for (ordinal in seq_along(groups)) {
                group <- members[groups[[ordinal]]]
                result <- .pa_average_inter_group(sp, group, spec, rmp,
                    prepared, cores, blockSize)
                effective_cores <- max(effective_cores, result$effective_cores)
                link_rows <- sort(unlist(units[groups[[ordinal]]],
                    use.names = FALSE))
                row <- .pa_average_row(pa, links, link_rows, set, ordinal,
                    result, level = "av_inter")
                rows[[length(rows) + 1L]] <- row
                graph[[length(rows)]] <- result$graph
                peaks[length(peaks) + 1L] <- list(result$peaks)
            }
        }
    }
    pa <- .pa_average_store(pa, rows, peaks, graph)
    pa@av_inter_params <- list(rmp = rmp, cores = cores, blockSize = blockSize,
        parameter_sets = specs)
    parameters <- list(rmp = rmp, cores = cores, effective_cores = effective_cores,
        blockSize = blockSize, configuration = list(parameter_sets = specs))
    .pa_processing_record(pa, "averageInterFragSpectra", parameters,
        started, input)
}

# Preserve positions in the producing spectrum before selecting passing peaks.
# Those peaks also determine the measured links eligible for inter-file grouping.
.pa_average_inter_input <- function(pa, sp, row, peak_list, columns, links, exact) {
    peaks <- .pa_av_frame(peak_list[[row]], columns, sp$col_types[row],
        sp$row_names[row])
    passing <- which(peaks$pass_flag)
    peaks <- peaks[passing, , drop = FALSE]
    peaks$memberPeak <- passing
    peaks$basis <- rep("derived", length(passing))
    id <- sp$spectrumId[row]
    peaks$sample <- rep(id, length(passing))
    if (!length(passing))
        return(list(peaks = peaks, clusters = list(), link_rows = integer()))
    if (!isTRUE(exact))
        stop("Historical intra-file average ", id, " lacks exact peak membership. ",
            "Repeat intra-file averaging before averaging across files.")
    origins <- averagePeakMembers(pa, id, peak = passing, depth = "original")
    # Occurrence IDs survive reordering of the current feature-link table.
    live <- match(origins$linkId, pa@grped_df$linkId)
    if (anyNA(live)) stop("Inter-file inputs refer to unavailable feature links.")
    clusters <- lapply(split(links$contributorIndex[live], origins$mergedPeak),
        unique)
    if (!all(as.character(passing) %in% names(clusters)))
        stop("Inter-file input has incomplete peak membership. ",
            "Repeat intra-file averaging before averaging across files.")
    list(peaks = peaks, clusters = clusters, link_rows = unique(live))
}

# Each selected file average is one numerical input. The worker records the
# immediate peak connections; older summary fields are projections of them.
.pa_average_inter_group <- function(sp, members, spec, rmp, prepared, cores, blockSize) {
    inputs <- unname(prepared[as.character(members)])
    tables <- lapply(seq_along(inputs), function(k) {
        p <- inputs[[k]]$peaks
        p$memberRank <- rep.int(k, nrow(p))
        p
    })
    table <- do.call(rbind, tables)
    table <- table[order(table$mz), , drop = FALSE]
    result <- average_spectra(table, indx = "sample", cores = cores,
        ppm = spec$ppm, minnum = spec$minnum, sumi = spec$sumi,
        minfrac = spec$minfrac, snthr = spec$snr, rathr = spec$ra,
        av_type = spec$av, record_inputs = TRUE, tolerance = spec$tolerance,
        blockSize = blockSize)
    peaks <- result$peaks
    if (rmp && !is.null(peaks))
        peaks <- peaks[peaks$pass_flag, , drop = FALSE]
    edges <- .pa_average_peak_edges(result, peaks, ranks = seq_along(inputs))
    used <- split(seq_len(nrow(edges)), factor(edges$mergedPeak,
        levels = seq_len(NROW(peaks))))
    clusters <- lapply(used, function(idx) {
        sources <- Map(function(rank, peak)
            inputs[[rank]]$clusters[[as.character(peak)]],
            edges$rank[idx], edges$memberPeak[idx])
        unique(unlist(sources, use.names = FALSE))
    })
    if (length(clusters))
        names(clusters) <- as.character(peaks$cl)
    else
        clusters <- list()
    list(peaks = peaks, clusters = clusters, effective_cores = result$effective_cores,
        graph = list(inputs = sp$spectrumId[members], edges = edges))
}

# Store configured averages after the averages of other levels.
.pa_average_store <- function(pa, rows, peaks, graph) {
    carried <- pa@avSpectra
    if ("msPurityAverageStatus" %in% Spectra::spectraVariables(carried)) {
        carried$msPurityAverageStatus <- NULL
        pa@avSpectra <- carried
    }
    if (!length(rows)) return(pa)
    scalar <- setdiff(names(rows[[1L]]), "msPurityAcquisitionValues")
    sv <- as.data.frame(lapply(stats::setNames(scalar, scalar), function(nm)
        unlist(lapply(rows, `[[`, nm), use.names = FALSE)),
        stringsAsFactors = FALSE)
    sv$is_null <- vapply(peaks, is.null, logical(1))
    sv$col_types <- vapply(peaks, function(df) if (is.null(df))
        NA_character_ else paste(vapply(df, function(x) class(x)[1], ""),
            collapse = ","), "")
    sv$row_names <- vapply(peaks, function(df) if (is.null(df))
        NA_character_ else .pa_encode_rownames(df), "")
    cols <- .pa_spectra_metadata(carried)$columns
    for (df in peaks) if (is.null(cols) && !is.null(df)) cols <- colnames(df)
    sp <- .pa_av_memory(sv, peaks, cols)
    sp <- .pa_av_scalar_columns(sp, seq_along(sp),
        lapply(rows, `[[`, "msPurityAcquisitionValues"))
    meta <- .pa_spectra_metadata(carried)
    meta$groups <- union(meta$groups, unique(sv$grpid))
    meta$columns <- cols
    out <- if (length(carried)) c(carried, sp) else sp
    .pa_spectra_metadata(out) <- meta
    pa@avSpectra <- out
    .pa_graph_store(pa, length(carried) + seq_along(rows), graph)
}
