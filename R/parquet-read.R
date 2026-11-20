# A read-only view of a single results generation. Cache tables only within a
# call, so a later read sees newly committed results or a relocated study.
.mzs_reader <- function(path, sourcePaths = NULL, m = MsBackendParquet::readManifest(path)) {
    if (!identical(m$role, "results"))
        .mzs_abort("unsupported", "'", path, "' is not a results dataset.")
    if (!identical(.mzs_semver_major(m$results$version), .MZS_RESULTS_MAJOR))
        .mzs_abort("unsupported", "Unsupported results version ", m$results$version)
    cache <- new.env(parent = emptyenv())
    table <- function(name, column = NULL, values = NULL) {
        entry <- m$results$tables[[name]]
        if (is.null(entry) || isTRUE(entry$omitted)) return(.mzs_empty_table(name))
        if (!is.null(column)) return(.mzs_read_table(path, m, name, column, values))
        if (is.null(cache[[name]])) cache[[name]] <- .mzs_read_table(path, m, name)
        cache[[name]]
    }
    list(path = path, manifest = m, table = table, sourcePaths = sourcePaths)
}

.mzs_original_ids <- function(reader, table, column) {
    si <- reader$table("source_identifier", "source_column", column)
    si <- si[si$target_table == table, ]
    stats::setNames(si$source_value, si$target_key)
}

.mzs_select <- function(value, selector) {
    if (is.null(selector)) rep(TRUE, length(value)) else value %in% selector
}

# Compute list lengths in Arrow. Older native runs have no scalar peak count.
.mzs_peak_counts <- function(path, manifest, info) {
    out <- integer(nrow(info))
    for (r in manifest$runs) {
        i <- which(info$run_id == r$run_id)
        if (!length(i)) next
        ds <- arrow::open_dataset(.mzs_run_dir(path, r))
        ids <- info$spectrum_id_[i]
        ds <- dplyr::filter(ds, spectrum_id_ %in% ids)
        counts <- dplyr::collect(dplyr::transmute(ds, spectrum_id_,
            n = arrow_list_value_length(mz)))
        k <- match(ids, counts$spectrum_id_)
        if (anyNA(k)) .mzs_abort("reference", "Missing spectrum in run ", r$run_id)
        out[i] <- counts$n[k]
    }
    out
}
utils::globalVariables(c("arrow_list_value_length", "mz"))

.mzs_absent_attempts <- function(reader) {
    ledger <- reader$table("loss_ledger")
    pattern <- "^av_(intra|inter|all): feature ([^,]+)(, sample (.+))?$"
    i <- which(ledger$reason == "Averaging was attempted and produced no spectrum." &
               grepl(pattern, ledger$source_construct))
    labels <- ledger$source_construct[i]
    groups <- .mzs_original_ids(reader, "feature", "grpid")
    sample <- sub(pattern, "\\4", labels)
    sample[sample == ""] <- NA_character_
    group <- sub(pattern, "\\2", labels)
    data.frame(loss_ledger_id_ = ledger$loss_ledger_id_[i],
        method = sub(pattern, "\\1", labels),
        feature_id_ = as.numeric(names(groups)[match(group, groups)]),
        grpid = as.numeric(group), sample = sample)
}

.mzs_stored_info <- function(reader, featureId, msLevel, method, runId) {
    m <- reader$manifest
    rr <- m$results$runs
    methods <- vapply(rr, `[[`, "", "method")
    runs <- names(rr)[.mzs_select(methods, method) & .mzs_select(names(rr), runId)]
    columns <- c("id", "ms_level", "selected_ion_mz", "scan_polarity",
        "x_mspurity_precursor_purity", names(.MZS_AV_CONFIGURED_TYPES))
    d <- .mzs_read_spectra_meta(reader$path, m, columns, runs %||% character())
    sf <- reader$table("spectrum_feature")
    sf <- sf[sf$ms2_source == "self", ]
    d$feature_id_ <- sf$feature_id_[match(d$spectrum_id_, sf$ms2_spectrum_id_)]
    d <- d[.mzs_select(d$feature_id_, featureId) &
           .mzs_select(d$ms_level, msLevel), , drop = FALSE]
    groups <- .mzs_original_ids(reader, "feature", "grpid")
    ft <- reader$table("feature")
    d$grpid <- as.numeric(groups[as.character(d$feature_id_)])
    d$retention_time <- ft$retention_time_in_seconds[match(d$feature_id_, ft$feature_id_)]
    d$method <- unname(methods[d$run_id])
    d$native_id <- as.character(d$id)
    d$id <- NULL
    d$ms_level <- as.integer(d$ms_level)
    d$spectrum_id_ <- as.numeric(d$spectrum_id_)
    d$precursor_mz <- as.numeric(d$selected_ion_mz)
    d$selected_ion_mz <- NULL
    for (name in names(.MZS_AV_CONFIGURED_TYPES))
        d[[name]] <- .mzs_as_r_type(d[[name]], .MZS_AV_CONFIGURED_TYPES[[name]])
    d$scan_polarity <- as.integer(d$scan_polarity)
    d$x_mspurity_precursor_purity <- as.numeric(d$x_mspurity_precursor_purity)
    d
}

.mzs_inventory <- function(reader, featureId, msLevel, method, runId) {
    m <- reader$manifest
    d <- .mzs_stored_info(reader, featureId, msLevel, method, runId)
    d$peak_count <- .mzs_peak_counts(reader$path, m, d)
    d$availability <- ifelse(d$peak_count == 0L, "empty", "available")
    d$availability <- as.character(d$availability)
    d$loss_ledger_id_ <- rep(NA_real_, nrow(d))
    d$sample <- rep(NA_character_, nrow(d))
    a <- .mzs_absent_attempts(reader)
    a <- a[.mzs_select(a$method, method) & .mzs_select(a$feature_id_, featureId), ]
    if (nrow(a)) {
        absent <- d[rep(NA_integer_, nrow(a)), , drop = FALSE]
        absent[names(a)] <- a
        # These ledger entries are emitted only for scalar MS2 averaging.
        absent$ms_level <- 2L
        absent$availability <- "unavailable"
        absent$run_id <- paste0("av_", a$method)
        intra <- a$method == "intra"
        fileids <- .mzs_original_ids(reader, "assay", "fileid")
        assay <- reader$table("assay")
        source <- assay$run_id[match(names(fileids)[match(a$sample, fileids)],
                                     assay$assay_id_)]
        absent$run_id[intra] <- ifelse(is.na(source[intra]), NA_character_,
                                       paste0("av_intra_", source[intra]))
        absent$run_id[!absent$run_id %in% names(m$results$runs)] <- NA_character_
        absent <- absent[.mzs_select(absent$run_id, runId) &
                         .mzs_select(absent$ms_level, msLevel), , drop = FALSE]
        d <- rbind(d, absent)
    }
    rownames(d) <- NULL
    d
}

# Keep the rank of each immediate selection. In an inter-file average, a
# selected intra-file average can contribute just a subset of its peaks.
.mzs_graph_reader <- function(reader, roots) {
    mm <- reader$table("merge_member")
    by_id <- split(seq_len(nrow(mm)), mm$merged_spectrum_id_)
    members <- function(id) mm[by_id[[as.character(id)]], , drop = FALSE]
    # Peak edges are needed only for exact original-measurement counts.
    # Resolving raw scans needs just the much smaller membership table.
    lineage <- function(ids) {
        seen <- numeric()
        while (length(ids)) {
            seen <- c(seen, ids)
            ix <- unlist(by_id[as.character(ids)], use.names = FALSE)
            parents <- mm$member_spectrum_id_[ix[mm$member_source[ix] == "self"]]
            ids <- setdiff(unique(parents), seen)
        }
        seen
    }
    family <- lineage(roots[!is.na(roots)])
    edges <- by_peak <- NULL
    # Original feature-link occurrences contributing to the specified stored
    # peaks. NULL peak positions request every stored peak of the root.
    used <- function(id, positions = NULL, visiting = numeric()) {
        if (id %in% visiting) .mzs_abort("semantic", "Cycle in averaging membership.")
        rows <- members(id)
        if (!nrow(rows) || !all(rows$members_complete %in% TRUE) ||
            is.null(rows$x_mspurity_contributes) || anyNA(rows$x_mspurity_contributes))
            return(NA_character_)
        if (is.null(by_peak)) {
            edges <<- reader$table("x_mspurity_merge_peak", "merged_spectrum_id_", family)
            by_peak <<- split(seq_len(nrow(edges)), edges$merged_spectrum_id_)
        }
        e <- edges[by_peak[[as.character(id)]], , drop = FALSE]
        if (!is.null(positions) && !all(positions %in% e$merged_peak_rank))
            return(NA_character_)
        if (!is.null(positions)) e <- e[e$merged_peak_rank %in% positions, ]
        if (!nrow(e)) {
            if (is.null(positions) && any(rows$x_mspurity_contributes))
                return(NA_character_)
            return(character())
        }
        ranks <- unique(e$member_rank)
        rows <- rows[match(ranks, rows$member_rank), ]
        if (anyNA(rows$member_source))
            .mzs_abort("semantic", "Unresolved member rank for averaged spectrum ", id)
        measured <- rows$member_source != "self"
        out <- character(sum(measured))
        if (any(measured)) {
            # Multiple peaks from one feature-link occurrence count once.
            out <- paste(id, ranks[measured], sep = ":")
            link <- rows$x_mspurity_feature_link_index[measured]
            if (length(link)) out[!is.na(link)] <- paste0("link:", link[!is.na(link)])
        }
        for (i in which(!measured)) {
            p <- e$member_peak_rank[e$member_rank == ranks[i]]
            if (anyNA(p)) return(NA_character_)
            out <- c(out, used(rows$member_spectrum_id_[i], p, c(visiting, id)))
        }
        unique(out)
    }
    list(lineage = lineage, members = members, used = used, indices = by_id, family = family)
}

.mzs_source_state <- function(reader, refs) {
    refs <- unique(refs[c("source", "run_id")])
    refs$status <- rep("resolved", nrow(refs))
    refs$condition <- vector("list", nrow(refs))
    for (i in seq_len(nrow(refs))) {
        key <- refs$source[i]
        err <- tryCatch({
            .mzs_open_source(reader$path, reader$manifest, key,
                             reader$sourcePaths, refs$run_id[i])
            NULL
        }, mzstack_capability = identity, mzstack_reference = identity,
           mzstack_stale = identity)
        if (!is.null(err)) {
            refs$status[i] <- if (inherits(err, "mzstack_stale")) "stale"
                else if (inherits(err, "mzstack_capability")) "external" else "missing"
            refs$condition[[i]] <- err
        }
    }
    rownames(refs) <- NULL
    refs
}

.mzs_member_refs <- function(mm) {
    raw <- mm[mm$member_source != "self", ]
    data.frame(source = raw$member_source, run_id = raw$member_run_id,
        spectrum_id_ = raw$member_spectrum_id_, native_id = raw$member_native_id,
        scan_number = raw$x_mspurity_member_scan_number %||% rep(NA_integer_, nrow(raw)))
}

.mzs_with_members <- function(reader, info, graph) {
    mm <- reader$table("merge_member")
    n <- nrow(info)
    info$contributing_inputs <- info$used_measurements <- rep(NA_integer_, n)
    refs <- vector("list", n)
    for (i in seq_len(n)) {
        id <- info$spectrum_id_[i]
        ids <- if (is.na(id)) numeric() else graph$lineage(id)
        ix <- unlist(graph$indices[as.character(ids)], use.names = FALSE)
        ref <- unique(.mzs_member_refs(mm[ix, ])[c("source", "run_id")])
        run <- reader$manifest$results$runs[[info$run_id[i]]]
        if (!is.null(run$source))
            ref <- unique(rbind(ref, data.frame(source = run$source, run_id = run$source_run_id)))
        refs[[i]] <- ref
        if (is.na(id)) next
        members <- graph$members(id)
        if (nrow(members) && all(members$members_complete %in% TRUE) &&
            !is.null(members$x_mspurity_contributes))
            info$contributing_inputs[i] <- sum(members$x_mspurity_contributes)
        if (identical(info$availability[i], "empty")) {
            info$contributing_inputs[i] <- info$used_measurements[i] <- 0L
        } else {
            used <- graph$used(id)
            if (!anyNA(used)) info$used_measurements[i] <- length(used)
        }
    }
    all_refs <- .mzs_member_refs(mm[FALSE, ])[c("source", "run_id")]
    if (length(refs)) all_refs <- do.call(rbind, refs)
    states <- .mzs_source_state(reader, all_refs)
    key <- function(d) paste(d$source, d$run_id, sep = "\r")
    info$source_runs <- lapply(refs, function(r)
        states[match(key(r), key(states)), , drop = FALSE])
    info
}

#' Read processed spectra and their provenance from a Parquet results dataset
#'
#' @param path Directory written by `createDatabase(format = "parquet")` or
#'   `convertSqliteToParquet()`.
#' @param featureId Exported numeric feature identifiers (`feature_id_`).
#' @param msLevel Integer MS levels.
#' @param method Averaging methods: `"intra"`, `"inter"` or `"all"`.
#' @param runId Character results run identifiers.
#' @param sourcePaths Named paths to relocated sources, keyed by manifest source.
#' @param includeRaw Also return the referenced study scans.
#' @details Values within a selector are alternatives. Selectors intersect,
#'   with NULL selecting all values. Stored spectra follow manifest run and
#'   spectrum order. An empty average has a stored ID and zero peaks. An
#'   unavailable attempt has no spectrum and appears only in metadata.
#' @return `processedSpectrumInfo()` returns a data frame including exported
#'   and xcms feature identifiers, acquisition metadata, availability, immediate
#'   contributing input counts and original measured occurrence counts for
#'   retained peaks (`contributing_inputs` and `used_measurements`). Unknown counts are NA. Its `source_runs` list column
#'   holds source/run, resolution status and the original typed condition
#'   (`mzstack_capability` for an external source, `mzstack_reference` for
#'   a missing or wrong dataset, and `mzstack_stale` for a changed run).
#'   `readProcessedSpectra()` returns a Spectra with an MsBackendParquet backend.
#'   With `includeRaw = TRUE`, it returns a list of `processed` Spectra and
#'   `raw` Spectra grouped by source. Unresolvable raw references signal their
#'   mzStack condition. Scan-only datasets return their selected measured scans.
#'   `readProcessedProvenance()` returns dataset activities, tool settings and
#'   history, software, selected feature links, measured inputs, membership,
#'   peak edges, source status and loss records. Dataset-wide activity and loss
#'   records retain their original scope.
#' @seealso [createDatabase()], [averageMembers()], [averagePeakMembers()]
#' @examples
#' db <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite",
#'                   package = "msPurity")
#' path <- tempfile("results-")
#' convertSqliteToParquet(db, path)
#' info <- processedSpectrumInfo(path)
#' feature <- info$feature_id_[1L]
#' sp <- readProcessedSpectra(path, featureId = feature)
#' cols <- intersect(c("mz", "intensity", "sn"), Spectra::peaksVariables(sp))
#' head(Spectra::peaksData(sp[1L], columns = cols)[[1L]])
#' head(readProcessedProvenance(path, featureId = feature)$mergeMembers)
#' @export
processedSpectrumInfo <- function(path, featureId = NULL, msLevel = NULL,
                                  method = NULL, runId = NULL, sourcePaths = NULL) {
    r <- .mzs_reader(path, sourcePaths)
    info <- .mzs_inventory(r, featureId, msLevel, method, runId)
    .mzs_with_members(r, info, .mzs_graph_reader(r, info$spectrum_id_))
}

.mzs_spectra_by_id <- function(path, ids) {
    be <- MsBackendParquet::MsBackendParquet()
    if (length(ids)) {
        be <- Spectra::backendInitialize(be, path = path)
        keys <- as.character(ids)
        if (anyNA(keys) || !all(keys %in% Spectra::spectraNames(be)))
            .mzs_abort("reference", "Dataset '", path,
                       "' does not contain all requested spectrum identifiers.")
        be <- be[keys]
    }
    Spectra::Spectra(be)
}

.mzs_processed_selection <- function(reader, featureId, msLevel, method, runId) {
    info <- .mzs_stored_info(reader, featureId, msLevel, method, runId)
    list(info = info, spectra = .mzs_spectra_by_id(reader$path, info$spectrum_id_))
}

.mzs_raw_spectra <- function(r, refs, nativeOnly = FALSE) {
    out <- list()
    for (key in unique(refs$source)) {
        rows <- refs[refs$source == key, ]
        src <- .mzs_open_source(r$path, r$manifest, key, r$sourcePaths, unique(rows$run_id))
        if (nativeOnly && any(vapply(src$manifest$runs, `[[`, "", "kind") != "native"))
            .mzs_abort("unsupported", "Matching scan peaks requires native runs in source '",
                       key, "'.")
        runs <- .mzs_runs_frame(src$manifest)
        k <- match(rows$run_id, runs$run_id)
        if (anyNA(k) || anyNA(rows$spectrum_id_) ||
            any(rows$spectrum_id_ < runs$uid_base[k] |
                rows$spectrum_id_ >= runs$uid_base[k] + runs$n_spectra[k]))
            .mzs_abort("reference", "Invalid spectrum reference in source '", key,
                       "'.", data = list(source = key))
        out[[key]] <- .mzs_spectra_by_id(src$path, unique(rows$spectrum_id_))
    }
    out
}

.mzs_read_scan_refs <- function(sc) {
    data.frame(source = sc$ms2_source, run_id = sc$ms2_run_id,
        spectrum_id_ = sc$ms2_spectrum_id_, native_id = sc$ms2_native_id,
        scan_number = sc$scan_number %||% sc$x_mspurity_scan_number %||%
            rep(NA_integer_, nrow(sc)))
}

# External sources do not have a dataset spectrum_id_. Include the run and
# fall back to the recorded native identifier to keep their links distinct.
.mzs_ref_key <- function(refs) {
    id <- ifelse(!is.na(refs$spectrum_id_), as.character(refs$spectrum_id_),
                 ifelse(!is.na(refs$native_id), paste0("native:", refs$native_id),
                        ifelse(!is.na(refs$scan_number), paste0("scan:", refs$scan_number),
                               NA_character_)))
    ifelse(is.na(id), NA_character_, paste(refs$source, refs$run_id, id, sep = "\r"))
}

.mzs_selected_inputs <- function(r, info, featureId, msLevel, method, runId) {
    mm <- r$table("merge_member")
    ids <- stats::na.omit(info$spectrum_id_)
    graph <- .mzs_graph_reader(r, ids)
    family <- graph$family
    members <- mm[mm$merged_spectrum_id_ %in% family, , drop = FALSE]
    members <- members[order(members$merged_spectrum_id_, members$member_rank), ]
    refs <- unique(.mzs_member_refs(members))
    scans <- r$table("x_mspurity_scan")
    links <- r$table("spectrum_feature")
    msn <- r$table("x_mspurity_msn_feature")
    scan_key <- .mzs_ref_key(.mzs_read_scan_refs(scans))
    link_key <- .mzs_ref_key(.mzs_read_scan_refs(links))
    if (!length(r$manifest$results$runs) && is.null(method) && is.null(runId)) {
        level <- scans$ms_level %||% rep(2L, nrow(scans))
        keep <- .mzs_select(level, msLevel)
        if (!is.null(featureId))
            keep <- keep & (scan_key %in% stats::na.omit(link_key[links$feature_id_ %in% featureId]) |
                scans$scan_annotation_id_ %in%
                msn$scan_annotation_id_[msn$feature_id_ %in% featureId])
        scans <- scans[keep, , drop = FALSE]
        refs <- unique(.mzs_read_scan_refs(scans))
    } else {
        scans <- scans[scan_key %in% stats::na.omit(.mzs_ref_key(refs)), , drop = FALSE]
    }
    links <- links[(links$ms2_source == "self" & links$ms2_spectrum_id_ %in% family) |
                   link_key %in% stats::na.omit(.mzs_ref_key(refs)), , drop = FALSE]
    msn <- msn[msn$scan_annotation_id_ %in% scans$scan_annotation_id_, , drop = FALSE]
    list(refs = refs, family = family, mergeMembers = members, measured = scans,
         featureLinks = links, msnFeatureLinks = msn)
}

#' @rdname processedSpectrumInfo
#' @export
readProcessedSpectra <- function(path, featureId = NULL, msLevel = NULL,
                                 method = NULL, runId = NULL, sourcePaths = NULL,
                                 includeRaw = FALSE) {
    r <- .mzs_reader(path, sourcePaths)
    selected <- .mzs_processed_selection(r, featureId, msLevel, method, runId)
    if (!includeRaw) return(selected$spectra)
    inputs <- .mzs_selected_inputs(r, selected$info, featureId, msLevel, method, runId)
    list(processed = selected$spectra, raw = .mzs_raw_spectra(r, inputs$refs))
}

#' @rdname processedSpectrumInfo
#' @export
readProcessedProvenance <- function(path, featureId = NULL, msLevel = NULL,
                                    method = NULL, runId = NULL, sourcePaths = NULL) {
    r <- .mzs_reader(path, sourcePaths)
    info <- .mzs_inventory(r, featureId, msLevel, method, runId)
    inputs <- .mzs_selected_inputs(r, info, featureId, msLevel, method, runId)
    edges <- r$table("x_mspurity_merge_peak", "merged_spectrum_id_", inputs$family)
    edges <- edges[order(edges$merged_spectrum_id_, edges$edge_rank), , drop = FALSE]
    tools <- r$table("tool_provenance")
    parse_json <- function(x) lapply(x, function(s)
        if (!is.na(s)) jsonlite::fromJSON(s, simplifyVector = FALSE))
    tools$settings <- parse_json(tools$parameters)
    tools$history <- parse_json(tools$x_mspurity_history %||% rep(NA_character_, nrow(tools)))
    loss <- r$table("loss_ledger")
    loss <- loss[!loss$reason %in% "Averaging was attempted and produced no spectrum." |
                 loss$loss_ledger_id_ %in% info$loss_ledger_id_, , drop = FALSE]
    activities <- r$manifest$provenance
    software <- r$table("software")
    filtered <- !all(vapply(list(featureId, msLevel, method, runId), is.null, logical(1)))
    if (filtered && !nrow(info) && !nrow(inputs$refs)) {
        activities <- activities[FALSE]
        tools <- tools[FALSE, , drop = FALSE]
        software <- software[FALSE, , drop = FALSE]
        loss <- loss[FALSE, , drop = FALSE]
    }
    list(activities = activities, tools = tools, software = software,
         outputs = info, measured = inputs$measured, featureLinks = inputs$featureLinks,
         msnFeatureLinks = inputs$msnFeatureLinks, mergeMembers = inputs$mergeMembers,
         mergePeaks = edges, sources = .mzs_source_state(r, inputs$refs), losses = loss)
}
