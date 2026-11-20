# An unknown leaf occupies one row. This lets a known contributor establish
# overlap even if another branch of an averaging graph has missing records.
.sm_measurements <- function(n = 1L) {
    x <- rep(NA_character_, n)
    data.frame(dataOrigin = x, datasetUid = x, run_id = x, spectrum_id_ = x,
        native_id = x, scanNumber = x)
}

.sm_origin <- function(x) {
    # A bare file name or an in-memory backend label has no shared namespace.
    file_uri <- grepl("^file://(localhost)?/", x, ignore.case = TRUE)
    x[file_uri] <- utils::URLdecode(sub("^file://(localhost)?", "", x[file_uri],
        ignore.case = TRUE))
    drive_uri <- file_uri & grepl("^/[A-Za-z]:/", x)
    x[drive_uri] <- substring(x[drive_uri], 2L)
    windows <- grepl("^[A-Za-z]:[/\\\\]|^\\\\\\\\", x)
    x[windows] <- chartr("\\", "/", x[windows])
    local <- startsWith(x, "/")
    absolute <- local | windows | grepl("^[A-Za-z][A-Za-z0-9+.-]*://", x)
    x[!absolute | is.na(absolute)] <- NA_character_
    rows <- which(local)
    paths <- unique(x[rows])
    if (length(paths))
        x[rows] <- normalizePath(paths, winslash = "/", mustWork = FALSE)[match(x[rows], paths)]
    x
}

.sm_live_leaves <- function(graph, key, peaks = NULL) {
    node <- match(key, graph$nodes$key)
    if (!isTRUE(graph$nodes$complete[node])) return(NA_integer_)
    edges <- .pa_graph_step(graph, key, peaks = TRUE)
    if (!is.null(peaks)) edges <- edges[edges$mergedPeak %in% peaks, , drop = FALSE]
    measured <- !is.na(edges$linkId)
    links <- edges$linkId[measured]
    for (parent in unique(edges$member[!measured])) {
        pos <- unique(edges$memberPeak[!measured & edges$member == parent])
        if (anyNA(pos)) links <- c(links, NA_integer_)
        pos <- pos[!is.na(pos)]
        if (length(pos)) links <- c(links, .sm_live_leaves(graph, parent, pos))
    }
    unique(links)
}

.sm_link_measurements <- function(refs) {
    value <- function(...) .sm_metadata_value(refs, ...)
    data.frame(dataOrigin = .sm_origin(value("dataOrigin")), datasetUid = value("datasetUid"),
        run_id = value("run_id"), spectrum_id_ = value("spectrum_id_"),
        native_id = value("spectrumId", "nativeId"), scanNumber = value("acquisitionNum"))
}

.sm_live_measurements <- function(sp) {
    # The public averagedSpectra view has already traced the retained peaks.
    # Per-row records also survive combining spectra from independent graphs.
    records <- if ("msPurityContributors" %in% Spectra::spectraVariables(sp))
        sp$msPurityContributors else vector("list", length(sp))
    out <- lapply(records, function(refs) {
        if (is.data.frame(refs) && nrow(refs)) .sm_link_measurements(refs)
    })
    missing <- which(vapply(out, is.null, logical(1)))
    if (length(missing)) {
        subset <- sp[missing]
        graph <- .pa_graph_data(subset)
        keys <- .pa_graph_sp_keys(subset, graph)
        out[missing] <- lapply(keys, function(key) {
            if (is.na(key)) return(.sm_measurements())
            ids <- .sm_live_leaves(graph, key)
            if (is.null(graph$links)) return(.sm_measurements(length(ids)))
            .sm_link_measurements(graph$links[match(ids, graph$links$linkId), , drop = FALSE])
        })
    }
    out
}

.sm_disk_measurements <- function(reader, ids) {
    graph <- .mzs_graph_reader(reader, ids)
    mm <- reader$table("merge_member")
    assay <- reader$table("assay")
    paths <- .mzs_original_ids(reader, "assay", "filepth")
    assay$origin <- unname(paths[as.character(assay$assay_id_)])
    assay_key <- paste(assay$source, assay$run_id, sep = "\r")
    sources <- reader$manifest$sources %||% list()
    names(sources) <- vapply(sources, `[[`, "", "key")
    lapply(ids, function(id) {
        ix <- graph$usedRows(id)
        known <- ix[!is.na(ix)]
        ref <- .mzs_member_refs(mm[known, , drop = FALSE])
        entries <- sources[as.character(ref$source)]
        origin <- assay$origin[match(paste(ref$source, ref$run_id, sep = "\r"), assay_key)]
        missing <- is.na(origin)
        origin[missing] <- vapply(entries[missing], function(s) s$location %||% NA_character_, "")
        out <- data.frame(dataOrigin = .sm_origin(origin),
            datasetUid = unname(vapply(entries, function(s) s$uid %||% NA_character_, "")),
            run_id = as.character(ref$run_id), spectrum_id_ = as.character(ref$spectrum_id_),
            native_id = as.character(ref$native_id), scanNumber = as.character(ref$scan_number))
        if (anyNA(ix)) out <- rbind(out, .sm_measurements())
        unique(out)
    })
}

.sm_source_reader <- function(sp, meta, role = "experimental") {
    n <- nrow(meta)
    type <- rep("measured", n)
    type[meta$consensus] <- "consensus"
    if (role == "library") {
        type[] <- "library"
        return(list(type = type, read = function(rows) rep(list(.sm_measurements()), length(rows))))
    }
    manifests <- list()
    storage <- uid <- rep(NA_character_, n)
    if (.pa_is_parquet(sp) && n) {
        storage <- as.character(Spectra::dataStorage(sp))
        for (path in unique(storage)) {
            manifest <- MsBackendParquet::readManifest(path)
            manifests[[path]] <- manifest
            rows <- which(storage == path)
            uid[rows] <- manifest$uid %||% NA_character_
            if (identical(manifest$role, "results"))
                type[rows[meta$runId[rows] %in% names(manifest$results$runs)]] <- "consensus"
        }
    }
    read <- function(rows) {
        out <- vector("list", length(rows))
        m <- which(type[rows] == "measured")
        if (length(m)) {
            i <- rows[m]
            scans <- data.frame(dataOrigin = .sm_origin(meta$origin[i]), datasetUid = uid[i],
                run_id = meta$runId[i], spectrum_id_ = meta$id[i], native_id = meta$nativeId[i],
                scanNumber = meta$scan[i])
            out[m] <- split(scans, seq_along(i))
        }
        consensus <- which(type[rows] == "consensus")
        live <- consensus[is.na(storage[rows[consensus]])]
        if (length(live)) out[live] <- .sm_live_measurements(sp[rows[live]])
        disk <- setdiff(consensus, live)
        for (group in split(disk, storage[rows[disk]])) {
            path <- storage[rows[group[1L]]]
            manifest <- manifests[[path]]
            if (identical(manifest$role, "results")) {
                reader <- .mzs_reader(path, m = manifest)
                out[group] <- .sm_disk_measurements(reader, as.numeric(meta$id[rows[group]]))
            } else out[group] <- rep(list(.sm_measurements()), length(group))
        }
        lapply(out, unique)
    }
    list(type = type, read = read)
}

# Each comparison is three-valued. A known different file rules out a
# shared scan even when one or both scan identifiers are missing.
.sm_shared <- function(query, reference) {
    shared_file <- shared_scan <- FALSE
    for (i in seq_len(nrow(query))) {
        q <- query[i, ]
        same_file <- q$dataOrigin == reference$dataOrigin
        same_dataset <- !is.na(q$datasetUid) & q$datasetUid == reference$datasetUid
        common_run_scope <- which(same_dataset & !is.na(q$run_id) & !is.na(reference$run_id))
        same_run <- q$run_id == reference$run_id[common_run_scope]
        within_run <- common_run_scope[same_run]
        other_run <- common_run_scope[!same_run]
        same_file[within_run] <- TRUE
        same_file[other_run[is.na(same_file[other_run])]] <- FALSE

        same_id <- q$scanNumber == reference$scanNumber
        native <- q$native_id == reference$native_id
        same_id[!is.na(native)] <- native[!is.na(native)]
        id <- q$spectrum_id_ == reference$spectrum_id_[within_run]
        same_id[within_run[!is.na(id)]] <- id[!is.na(id)]
        shared_file <- shared_file | any(same_file)
        # A matching scan identifier only has meaning inside a known file.
        same_scan <- same_file & same_id
        same_scan[is.na(same_file)] <- NA
        shared_scan <- shared_scan | any(same_scan)
        if (isTRUE(shared_scan)) break
    }
    c(sharedMeasurement = shared_scan, sharedSourceFile = shared_file)
}

.sm_hit_sources <- function(hits, query, reference) {
    if (is.null(hits)) {
        empty <- data.frame(side = character(), inputIndex = integer(), .sm_measurements(0L),
            fileComplete = logical(), measurementComplete = logical())
        return(list(hits = NULL, table = empty))
    }
    participants <- function(side, input, rows) {
        scans <- input$read(rows)
        tables <- Map(function(i, d) {
            file <- !is.na(d$dataOrigin) | (!is.na(d$datasetUid) & !is.na(d$run_id))
            scan <- !is.na(d$native_id) | !is.na(d$scanNumber) |
                (!is.na(d$datasetUid) & !is.na(d$run_id) & !is.na(d$spectrum_id_))
            data.frame(side = rep(side, nrow(d)), inputIndex = rep(i, nrow(d)),
                d, fileComplete = rep(all(file), nrow(d)),
                measurementComplete = rep(all(file & scan), nrow(d)))
        }, rows, scans)
        list(scans = scans, table = do.call(rbind, tables))
    }
    q <- unique(hits$queryIndex)
    r <- unique(hits$referenceIndex)
    left <- participants("query", query, q)
    right <- participants("reference", reference, r)
    table <- rbind(left$table, right$table)
    rownames(table) <- NULL
    hits$queryType <- query$type[hits$queryIndex]
    hits$referenceType <- reference$type[hits$referenceIndex]
    if (all(hits$referenceType == "library")) {
        hits$sharedMeasurement <- hits$sharedSourceFile <- NA
    } else {
        qi <- match(hits$queryIndex, q)
        ri <- match(hits$referenceIndex, r)
        flags <- vapply(seq_len(nrow(hits)), function(i)
            .sm_shared(left$scans[[qi[i]]], right$scans[[ri[i]]]), logical(2L))
        hits$sharedMeasurement <- flags[1L, ]
        hits$sharedSourceFile <- flags[2L, ]
    }
    list(hits = hits, table = table)
}
