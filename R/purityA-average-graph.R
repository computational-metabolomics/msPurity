# Object-local spectrum keys and immediate input records. Acquisition PIDs
# remain available even when they cannot themselves be used as integer keys.
.pa_graph_new <- function(pa) {
    pid <- unique(c(if (length(pa@spectra)) pa@spectra$pid,
        if (length(pa@fragSpectra)) pa@fragSpectra$pid) %||% integer())
    if (any(!is.finite(pid))) stop("A measured spectrum has no finite pid.")
    direct <- all(pid >= 1 & pid == floor(pid) & pid < .Machine$integer.max)
    keys <- if (direct) as.integer(pid) else seq_along(pid)
    list(version = 2L, lastKey = max(c(0L, keys)), lastLink = 0L,
        scans = data.frame(key = keys, pid = pid),
        nodes = data.frame(key = integer(), spectrumId = character(), nPeaks = integer(),
            complete = logical(), exactLeaves = logical(), activity = character()),
        mergeMembers = data.frame(merged = integer(), rank = integer(),
            member = integer(), linkId = integer()),
        mergePeaks = data.frame(merged = integer(), mergedPeak = integer(),
            rank = integer(), memberPeak = integer(), basis = factor(character(), levels = c(
                "source_reference", "filter_input", "derived", "unknown"))),
        links = NULL)
}

.pa_graph_keys <- function(last, n) {
    if (as.double(last) + n > .Machine$integer.max)
        stop("The object's integer provenance key space is exhausted.")
    last + seq_len(n)
}

# Earlier graph snapshots kept the integer key on each averaged spectrum.
# Expose its decimal representation when the row has no existing spectrumId.
.pa_graph_sp_ids <- function(sp) {
    frame <- is.data.frame(sp)
    vars <- if (frame) names(sp) else Spectra::spectraVariables(sp)
    ids <- if ("spectrumId" %in% vars) as.character(sp$spectrumId)
        else rep(NA_character_, if (frame) nrow(sp) else length(sp))
    ids[!is.na(ids) & !nzchar(ids)] <- NA_character_
    if ("msPurityKey" %in% vars) {
        missing <- is.na(ids)
        ids[missing] <- as.character(sp$msPurityKey[missing])
    }
    ids
}

.pa_graph_upgrade <- function(g, sp) {
    if (is.null(g) || identical(g$version, 2L)) return(g)
    if (!identical(g$version, 1L)) stop("Unknown averaging graph version.")
    ids <- as.character(g$nodes$key)
    if ("msPurityKey" %in% Spectra::spectraVariables(sp)) {
        pos <- match(sp$msPurityKey, g$nodes$key)
        known <- !is.na(pos)
        ids[pos[known]] <- .pa_graph_sp_ids(sp)[known]
    }
    if (anyDuplicated(ids)) stop("Historical average spectrumId values are ambiguous.")
    g$nodes$spectrumId <- ids
    g$version <- 2L
    g
}

.pa_graph_sp_keys <- function(sp, g) {
    g <- .pa_graph_upgrade(g, sp)
    if (is.null(g)) return(rep(NA_integer_, length(sp)))
    key <- g$nodes$key[match(.pa_graph_sp_ids(sp), g$nodes$spectrumId)]
    # The same native ID can occur in a measured row of a mixed Spectra.
    vars <- Spectra::spectraVariables(sp)
    if ("av_level" %in% vars) {
        key[is.na(sp$av_level)] <- NA_integer_
    } else if ("msPurityProcessingId" %in% vars) {
        activity <- g$nodes$activity[match(key, g$nodes$key)]
        key[is.na(activity) | is.na(sp$msPurityProcessingId) |
            activity != sp$msPurityProcessingId] <- NA_integer_
    } else {
        key[] <- NA_integer_
    }
    key
}

.pa_graph_set_ids <- function(sp, ids) {
    sp$spectrumId <- ids
    if ("msPurityKey" %in% Spectra::spectraVariables(sp)) sp$msPurityKey <- NULL
    sp
}

# A preserved historical label can already be a decimal string. Leave it
# intact and allocate a fresh contiguous block for this execution.
.pa_graph_average_keys <- function(last, n, reserved) {
    repeat {
        keys <- .pa_graph_keys(last, n)
        collisions <- keys[as.character(keys) %in% reserved]
        if (!length(collisions)) return(keys)
        last <- max(collisions)
    }
}

# Allocate identifiers for older inputs without assigning unrecorded edges.
.pa_graph_import <- function(pa) {
    sp <- pa@avSpectra
    graph <- .pa_graph_upgrade(pa@params$averageGraph, sp) %||% .pa_graph_new(pa)
    keys <- .pa_graph_sp_keys(sp, graph)
    ids <- .pa_graph_sp_ids(sp)
    if (anyDuplicated(ids[!is.na(ids)]))
        stop("Historical average spectrumId values are ambiguous.")
    missing <- which(is.na(keys))
    if (length(missing)) {
        fresh <- .pa_graph_average_keys(graph$lastKey, length(missing), ids)
        graph$lastKey <- max(fresh)
        unnamed <- is.na(ids[missing])
        ids[missing[unnamed]] <- as.character(fresh[unnamed])
        activity <- if ("msPurityProcessingId" %in% Spectra::spectraVariables(sp))
            sp$msPurityProcessingId else rep(NA_character_, length(sp))
        cols <- c("key", "spectrumId", "nPeaks", "complete", "exactLeaves", "activity")
        graph$nodes <- rbind(graph$nodes[, cols, drop = FALSE], data.frame(key = fresh,
            spectrumId = ids[missing], nPeaks = .msp_lengths(sp[missing]),
            complete = FALSE, exactLeaves = FALSE, activity = activity[missing]))
        graph <- .pa_graph_index(graph)
    }
    if (length(sp)) pa@avSpectra <- .pa_graph_set_ids(sp, ids)
    pa@params$averageGraph <- graph
    pa
}

# Mint each occurrence at association, before filtering or sorting can remove
# its original position. Older saved objects receive keys on first use.
.pa_graph_mint_links <- function(pa) {
    graph <- pa@params$averageGraph %||% .pa_graph_new(pa)
    id <- .pa_graph_keys(graph$lastLink, nrow(pa@grped_df))
    pa@grped_df$linkId <- id
    if (length(id)) graph$lastLink <- max(id)
    pa@params$averageGraph <- graph
    pa
}

.pa_graph_link_rows <- function(g, refs) {
    idx <- match(g$linkId, refs$linkId)
    fields <- intersect(c("pid", "grpid", "cid", "sample"), names(g))
    same <- vapply(fields, function(nm) identical(as.character(refs[[nm]][idx]),
        as.character(g[[nm]])), logical(1))
    if (anyNA(idx) || anyDuplicated(g$linkId) || !all(same))
        stop("Feature-link identities no longer match their association. Repeat frag4feature().")
    idx
}

.pa_graph_register_links <- function(pa) {
    if (!nrow(pa@grped_df)) return(pa)
    if (is.null(pa@grped_df$linkId)) pa <- .pa_graph_mint_links(pa)
    g <- pa@grped_df
    graph <- pa@params$averageGraph
    if (!is.null(graph$links)) {
        idx <- .pa_graph_link_rows(g, graph$links)
        if ("msPurityProcessingId" %in% Spectra::spectraVariables(pa@fragSpectra)) {
            tip <- pa@fragSpectra$msPurityProcessingId
            graph$links$msPurityProcessingId[idx] <- tip[match(g$pid, pa@fragSpectra$pid)]
            pa@params$averageGraph <- graph
        }
        return(pa)
    }
    assoc <- pa@params$processing$active$frag4feature %||% NA_character_
    original <- if (!is.na(assoc)) pa@params$processing$records[[assoc]]$output$feature_links
    if (is.null(original$linkId)) original <- g
    link_index <- .pa_graph_link_rows(g, original)
    within_feature <- stats::ave(seq_len(nrow(original)), as.character(original$grpid),
        FUN = seq_along)
    rows <- match(g$pid, pa@fragSpectra$pid)
    if (anyNA(rows)) stop("A feature link has no measured input spectrum.")
    refs <- data.frame(linkId = g$linkId, associationId = assoc,
        featureLinkIndex = link_index, contributorIndex = as.integer(within_feature[link_index]),
        grpid = as.character(g$grpid), sample = if (is.null(g$sample))
            NA_character_ else as.character(g$sample), pid = g$pid)
    refs$cid <- if (is.null(g$cid)) NA_integer_ else g$cid
    for (nm in intersect(c("fileid", "seqNum", "spectrumId", "nativeId",
            "acquisitionNum", "msPurityProcessingId"),
            Spectra::spectraVariables(pa@fragSpectra)))
        refs[[nm]] <- pa@fragSpectra[[nm]][rows]
    graph$links <- .pa_graph_qualify_links(refs, pa, rows)
    pa@params$averageGraph <- graph
    pa
}

# Keep source identities on the graph so its measured leaves can be compared
# with scans outside this object, even after subsetting the averages.
.pa_graph_qualify_links <- function(refs, pa, rows) {
    for (nm in intersect(c("dataOrigin", "run_id", "spectrum_id_"),
            Spectra::spectraVariables(pa@fragSpectra)))
        refs[[nm]] <- pa@fragSpectra[[nm]][rows]
    refs$datasetUid <- rep(.pa_study(pa)$uid %||% NA_character_, nrow(refs))
    refs
}

.pa_graph_index <- function(g) {
    # Rows are appended in key order. Offsets keep small queries local even
    # when other features have millions of recorded peak edges.
    for (field in c("mergeMembers", "mergePeaks")) {
        counts <- tabulate(match(g[[field]]$merged, g$nodes$key), nrow(g$nodes))
        g$nodes[[paste0(field, "End")]] <- cumsum(counts)
    }
    g
}

.pa_graph_prune <- function(pa) {
    g <- .pa_graph_upgrade(pa@params$averageGraph, pa@avSpectra)
    if (is.null(g)) return(pa)
    keep <- .pa_graph_sp_keys(pa@avSpectra, g)
    if (all(g$nodes$key %in% keep)) return(pa)
    g$nodes <- g$nodes[g$nodes$key %in% keep, , drop = FALSE]
    for (nm in c("mergeMembers", "mergePeaks")) {
        tab <- g[[nm]]
        g[[nm]] <- tab[tab$merged %in% keep, , drop = FALSE]
    }
    pa@params$averageGraph <- .pa_graph_index(g)
    pa
}

.pa_graph_positions <- function(x, upper) {
    is.numeric(x) && !anyNA(x) && all(x >= 1 & x <= upper & x == floor(x))
}

# Each group supplies its selected measured link rows or earlier spectrum
# keys plus the exact edges recorded by the numerical worker in output order.
.pa_graph_store <- function(pa, rows, groups) {
    if (length(rows) != length(groups)) stop("Missing averaging input records.")
    if (!length(rows)) return(.pa_graph_prune(pa))
    g <- .pa_graph_upgrade(pa@params$averageGraph, pa@avSpectra) %||% .pa_graph_new(pa)
    pa@params$averageGraph <- g
    if (any(vapply(groups, function(x) !is.null(x$links), logical(1))))
        pa <- .pa_graph_register_links(pa)
    pa <- .pa_graph_prune(pa)
    g <- pa@params$averageGraph
    sp <- pa@avSpectra
    ids_out <- .pa_graph_sp_ids(sp)
    keys <- .pa_graph_average_keys(g$lastKey, length(rows), ids_out)
    sizes <- .msp_lengths(sp[rows])
    complete <- exact <- logical(length(rows))
    members <- edges <- vector("list", length(rows))
    for (i in seq_along(groups)) {
        x <- groups[[i]]
        complete[i] <- isTRUE(x$complete %||% TRUE)
        if (!is.null(x$links)) {
            if (!.pa_graph_positions(x$links, nrow(pa@grped_df)))
                stop("Invalid measured feature-link reference.")
            ids <- pa@grped_df$linkId[x$links]
            member <- g$scans$key[match(pa@grped_df$pid[x$links], g$scans$pid)]
            exact[i] <- complete[i]
        } else {
            input <- x$inputs %||% integer()
            parent <- if (is.character(input)) match(input, g$nodes$spectrumId)
                else match(input, g$nodes$key)
            member <- g$nodes$key[parent]
            if (anyNA(parent)) stop("An averaging input key is unresolved.")
            ids <- rep(NA_integer_, length(member))
            exact[i] <- complete[i] && all(g$nodes$exactLeaves[parent])
        }
        if (anyNA(member)) stop("An averaging input has no spectrum key.")
        if (length(member)) members[[i]] <- data.frame(merged = keys[i],
            rank = seq_along(member), member = member, linkId = ids)
        e <- x$edges
        if (!is.null(e) && nrow(e)) {
            if (anyNA(e$basis) || any(!e$basis %in% levels(g$mergePeaks$basis)))
                stop("Unknown peak coordinate basis.")
            if (!.pa_graph_positions(e$rank, length(member)) ||
                !.pa_graph_positions(e$mergedPeak, sizes[i]) ||
                !.pa_graph_positions(stats::na.omit(e$memberPeak), .Machine$integer.max) ||
                any(is.na(e$memberPeak) & e$basis != "unknown"))
                stop("Invalid position in recorded averaging peak edges.")
            if (is.null(x$links) && any(e$memberPeak >
                    g$nodes$nPeaks[parent[e$rank]], na.rm = TRUE))
                stop("A referenced peak lies outside its input average.")
            e$basis <- factor(e$basis, levels = levels(g$mergePeaks$basis))
            edges[[i]] <- data.frame(merged = keys[i],
                e[, c("mergedPeak", "rank", "memberPeak", "basis"), drop = FALSE])
        }
    }
    g$mergeMembers <- rbind(g$mergeMembers, do.call(rbind, members))
    g$mergePeaks <- rbind(g$mergePeaks, do.call(rbind, edges))
    # Offsets are rebuilt after appending to the two tables.
    cols <- c("key", "spectrumId", "nPeaks", "complete", "exactLeaves", "activity")
    g$nodes <- rbind(g$nodes[, cols, drop = FALSE], data.frame(key = keys,
        spectrumId = as.character(keys), nPeaks = sizes, complete = complete, exactLeaves = exact,
        activity = NA_character_))
    g$lastKey <- max(keys)
    ids_out[rows] <- as.character(keys)
    sp <- .pa_graph_set_ids(sp, ids_out)
    # Historical rows retain their older contributor records. For new exact
    # results the public Spectra view obtains these lists from the graph.
    for (nm in intersect(c("msPurityContributors", "msPurityClusterContributors",
            "msPurityAcquisitionValues", "msPurityAcquisitionStatus",
            "msPurityPrecursorPath", "msPurityAncestorConditions"),
            Spectra::spectraVariables(sp))) {
        val <- sp[[nm]]
        val[rows[exact]] <- rep(list(NULL), sum(exact))
        if (all(vapply(val, is.null, logical(1)))) sp[[nm]] <- NULL
        else sp[[nm]] <- val
    }
    pa@avSpectra <- sp
    pa@params$averageGraph <- .pa_graph_index(g)
    pa
}

.pa_graph_data <- function(object) {
    if (methods::is(object, "purityA"))
        return(.pa_graph_upgrade(object@params$averageGraph, object@avSpectra))
    if (!methods::is(object, "Spectra")) stop("Expected a purityA or Spectra object.")
    meta <- .pa_spectra_metadata(object)
    graphs <- meta[names(meta) %in% "msPurityAverageGraph"]
    if (!length(graphs)) return(NULL)
    g <- graphs[[1L]]
    if (length(graphs) > 1L && !all(vapply(graphs, identical, logical(1), g)))
        stop("The spectra refer to different provenance graphs.")
    # Some backends keep only the first object's metadata on concatenation.
    g <- .pa_graph_upgrade(g, object)
    keys <- .pa_graph_sp_keys(object, g)
    rows <- which(!is.na(keys))
    known <- g$nodes$activity[match(keys[rows], g$nodes$key)]
    activity <- if ("msPurityProcessingId" %in% Spectra::spectraVariables(object))
        object$msPurityProcessingId[rows] else rep(NA_character_, length(rows))
    if (any(is.na(known) != is.na(activity)) || any(known != activity, na.rm = TRUE))
        stop("A spectrum's processing activity is absent from or differs from this graph.")
    g
}

.pa_graph_rows <- function(g, table, node) {
    end <- g$nodes[[paste0(table, "End")]]
    start <- if (node == 1L) 0L else end[node - 1L]
    g[[table]][start + seq_len(end[node] - start), , drop = FALSE]
}

.pa_graph_step <- function(g, key, peaks = FALSE) {
    node <- match(key, g$nodes$key)
    if (is.na(node)) stop("Unresolved derived spectrum key: ", key)
    if (!g$nodes$complete[node])
        stop("Immediate averaging provenance is unknown for historical spectrum ",
            g$nodes$spectrumId[node], ".")
    m <- .pa_graph_rows(g, "mergeMembers", node)
    if (any(is.na(m$linkId) & m$member >= key))
        stop("A derived input must refer to an earlier spectrum key.")
    e <- .pa_graph_rows(g, "mergePeaks", node)
    if (!peaks) {
        m$contributes <- m$rank %in% e$rank
        return(m)
    }
    e$basis <- as.character(e$basis)
    e$member <- m$member[e$rank]
    e$linkId <- m$linkId[e$rank]
    e
}

# Expand each distinct intermediate spectrum once at this level. Matching on
# its stored peak position excludes ancestors whose only peak was dropped.
.pa_graph_walk <- function(g, key, peaks, original, peak = NULL) {
    out <- .pa_graph_step(g, key, peaks)
    if (peaks && !is.null(peak)) {
        if (!.pa_graph_positions(peak, g$nodes$nPeaks[match(key, g$nodes$key)]))
            stop("Peak position lies outside the requested average.")
        out <- out[out$mergedPeak %in% peak, , drop = FALSE]
    }
    if (original) while (anyNA(out$linkId)) {
        derived <- which(is.na(out$linkId))
        parts <- split(derived, factor(out$member[derived],
            levels = unique(out$member[derived])))
        expanded <- lapply(parts, function(ix) {
            parent <- out$member[ix[1L]]
            child <- .pa_graph_step(g, parent, peaks)
            if (peaks) {
                if (anyNA(out$memberPeak[ix]))
                    stop("The intermediate peak position is unknown.")
                by_pos <- split(seq_len(nrow(child)), child$mergedPeak)
                hits <- by_pos[as.character(out$memberPeak[ix])]
                j <- unlist(hits, use.names = FALSE)
                ix <- rep(ix, lengths(hits))
                child <- child[j, , drop = FALSE]
                child$mergedPeak <- out$mergedPeak[ix]
            } else {
                ix <- rep(ix, each = nrow(child))
                child <- child[rep(seq_len(nrow(child)), length.out = length(ix)), , drop = FALSE]
            }
            child$merged <- rep.int(key, length(ix))
            child$.order <- ix
            child
        })
        keep <- which(!is.na(out$linkId))
        measured <- out[keep, , drop = FALSE]
        measured$.order <- keep
        out <- do.call(rbind, c(list(measured), unname(expanded)))
        out <- out[order(out$.order), , drop = FALSE]
        out$.order <- NULL
    }
    if (original && !peaks) out$contributes <- NULL
    rownames(out) <- NULL
    out
}

# Resolve public identifiers against the supplied subset before traversing its graph.
.pa_graph_resolve <- function(sp, g, spectrumId, peaks = FALSE) {
    if (!is.character(spectrumId) || !length(spectrumId) || anyNA(spectrumId) ||
            any(!nzchar(spectrumId)))
        stop("Supply a non-missing character spectrumId for each requested average.")
    if (is.null(g)) stop("Averaging provenance is unknown for this object.")
    keys <- .pa_graph_sp_keys(sp, g)
    ids <- .pa_graph_sp_ids(sp)
    absent <- setdiff(spectrumId, ids[!is.na(keys)])
    if (length(absent)) {
        if ("av_level" %in% Spectra::spectraVariables(sp) &&
                absent[1L] %in% ids[!is.na(sp$av_level)])
            stop("Averaging provenance is unknown for spectrumId ", absent[1L], ".")
        stop("Select an averaged spectrumId present in the supplied object.")
    }
    key <- g$nodes$key[match(spectrumId, g$nodes$spectrumId)]
    if (peaks) {
        selected <- which(keys %in% key)
        spectra <- sp[selected]
        changed <- any(.msp_lengths(spectra) !=
            g$nodes$nPeaks[match(keys[selected], g$nodes$key)])
        if (!changed && .pa_spectra_has_queue(spectra))
            changed <- !identical(.msp_peaks(spectra),
                .msp_peaks(Spectra::reset(spectra)))
        if (changed)
            stop("Peak processing of averaged spectra is not tracked. ",
                "Query the unmodified averaging result.")
    }
    key
}

# Expose the public ID at each end of an edge. Measured IDs are file-local
# and may be missing, so keep their scan, file and feature-link references.
.pa_graph_member_ids <- function(g, out) {
    link <- match(out$linkId, g$links$linkId)
    measured <- !is.na(out$linkId)
    out$merged <- g$nodes$spectrumId[match(out$merged, g$nodes$key)]
    member <- g$nodes$spectrumId[match(out$member, g$nodes$key)]
    member[measured] <- if (is.null(g$links$spectrumId)) NA_character_
        else as.character(g$links$spectrumId[link[measured]])
    out$member <- member
    for (nm in c("pid", "associationId", "featureLinkIndex", "contributorIndex",
            "cid", "grpid", "fileid", "sample")) {
        value <- g$links[[nm]]
        out[[nm]] <- if (is.null(value)) rep(NA, nrow(out)) else value[link]
    }
    out
}

.pa_graph_query <- function(object, spectrumId, peaks, depth, peak = NULL) {
    if (length(spectrumId) != 1L) stop("Select one averaged spectrumId.")
    sp <- if (methods::is(object, "purityA")) object@avSpectra else object
    g <- .pa_graph_data(object)
    key <- .pa_graph_resolve(sp, g, spectrumId, peaks = peaks)
    .pa_graph_member_ids(g, .pa_graph_walk(g, key, peaks, depth == "original", peak))
}

#' Trace the immediate or original inputs to an averaged spectrum
#'
#' A new average receives a decimal string in `spectrumId`, unique among the
#' averages of the same object. Independent objects can assign the same ID.
#' Measured spectra keep their file-native identifiers. Derived IDs persist
#' after serialization and repeated averaging on the same object assigns fresh
#' IDs. `msPurityProcessingId` identifies the execution.
#'
#' `averageMembers()` follows selected spectrum inputs, including those
#' supplying no stored peak. At immediate depth, `contributes` reports
#' whether that input supplied a peak to the stored output.
#' `averagePeakMembers()` follows only the input
#' peaks of the requested stored output peaks. When failing peaks were kept,
#' use their `pass_flag` to select passing positions. An original selected
#' input can lack a path to any retained final peak.
#'
#' @param object A purityA object or averaged Spectra returned by
#'   `averagedSpectra(pa, legacy = FALSE)`, including a subset.
#' @param spectrumId Character `spectrumId` of an average in the supplied object.
#' @param depth `"immediate"` for one step or `"original"` to follow measured
#'   inputs through earlier averages.
#' @param peak Optional vector of 1-based stored output peak positions.
#'   The default includes all stored peaks.
#' @return A data frame of memberships or peak connections. `merged` is the
#'   requested `spectrumId` and `member` is the input's `spectrumId`. For
#'   measured inputs the native ID can repeat across files or be missing.
#'   `pid` identifies the measured scan, `fileid` its file and
#'   `linkId` the selected feature-link occurrence. `rank` is the input's
#'   rank at its immediately consuming step. Peak connections additionally
#'   give `mergedPeak` in the requested spectrum, `memberPeak` at the input,
#'   and its coordinate `basis`. A filter-input position need not equal a
#'   raw acquisition position. An unknown input peak position is reported as
#'   missing. Unrecorded historical edges cannot be traversed.
#' @details Peak queries require the unmodified averaged peak arrays. Subsequent
#'   peak processing with Spectra does not update these records. To inspect
#'   particular peaks, pass their original averaged positions to \code{peak}.
#'   \code{averageMembers()} queries the recorded selection even when peak
#'   arrays are omitted. Spectrum subsetting and reordering preserve the records.
#' @examples
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'     "2_frag4feature_pa.rds", package = "msPurity"))
#' pa <- averageAllFragSpectra(pa, cores = 1)
#' s <- averagedSpectra(pa, legacy = FALSE)
#' averageMembers(s, s$spectrumId[1L])
#' averagePeakMembers(s, s$spectrumId[1L], depth = "original")
#' @export
averageMembers <- function(object, spectrumId, depth = c("immediate", "original")) {
    .pa_graph_query(object, spectrumId, FALSE, match.arg(depth))
}

#' @rdname averageMembers
#' @export
averagePeakMembers <- function(object, spectrumId, peak = NULL,
        depth = c("immediate", "original")) {
    .pa_graph_query(object, spectrumId, TRUE, match.arg(depth), peak)
}

# Reconstruct the Spectra list columns from the same retained graph leaves.
# Historical rows without exact edges keep their existing records.
.pa_graph_view <- function(pa, sp = pa@avSpectra) {
    g <- .pa_graph_data(pa)
    if (is.null(g)) return(sp)
    historical <- NROW(g$links) > 0L && is.null(g$links$datasetUid)
    if (historical)
        g$links <- .pa_graph_qualify_links(g$links, pa,
            match(g$links$pid, pa@fragSpectra$pid))
    old_keys <- "msPurityKey" %in% Spectra::spectraVariables(sp)
    if (old_keys) sp <- .pa_graph_set_ids(sp, .pa_graph_sp_ids(sp))
    if (historical || old_keys) {
        meta <- .pa_spectra_metadata(sp)
        meta$msPurityAverageGraph <- g
        .pa_spectra_metadata(sp) <- meta
    }
    if (!length(sp)) return(sp)
    keys <- .pa_graph_sp_keys(sp, g)
    exact <- which(g$nodes$exactLeaves[match(keys, g$nodes$key)] %in% TRUE)
    if (!length(exact)) return(sp)
    vars <- Spectra::spectraVariables(sp)
    field <- function(nm) if (nm %in% vars) sp[[nm]] else rep(list(NULL), length(sp))
    contributors <- field("msPurityContributors")
    clusters <- field("msPurityClusterContributors")
    values <- field("msPurityAcquisitionValues")
    status <- field("msPurityAcquisitionStatus")
    paths <- field("msPurityPrecursorPath")
    ancestors <- field("msPurityAncestorConditions")
    configured <- .pa_av_configured(sp)
    # The scalar legacy adapter represents a feature with no intra-file
    # samples by one NULL entry. An actual empty average has an empty list.
    no_sample <- if (all(c("av_level", "sample", "is_null") %in% vars))
        !configured & sp$av_level == "av_intra" & is.na(sp$sample) & sp$is_null
    else rep(FALSE, length(sp))
    higher <- configured & Spectra::msLevel(sp) %in% c(3L, 4L)
    refs <- g$links
    refcols <- setdiff(names(refs), c("cid", "associationId", "linkId"))
    source <- match(refs$pid, pa@fragSpectra$pid)
    if (any(higher[exact])) {
        links <- .pa_average_links(pa)
        live <- match(refs$linkId, pa@grped_df$linkId)
    }
    labels <- if ("cl" %in% Spectra::peaksVariables(sp)) Spectra::peaksData(
        sp[exact], columns = "cl", BPPARAM = .msp_bpparam()) else NULL
    for (j in seq_along(exact)) {
        i <- exact[j]
        out <- .pa_graph_walk(g, keys[i], TRUE, TRUE)
        ref <- match(out$linkId, refs$linkId)
        orig <- refs$contributorIndex[ref]
        pos <- unique(ref[order(orig)])
        tab <- refs[pos, refcols, drop = FALSE]
        rownames(tab) <- NULL
        contributors[i] <- list(tab)
        if (anyNA(source[pos])) stop("An average refers to unavailable measured spectra.")
        acquisition <- .pa_acquisition_values(pa@fragSpectra, source[pos], configured[i])
        values[i] <- list(acquisition)
        status[i] <- list(.pa_acquisition_status(acquisition))
        if (higher[i]) {
            if (anyNA(live[pos])) stop("An average refers to unavailable feature links.")
            summary <- .pa_average_path_summary(links, live[pos])
            paths[i] <- list(summary$path)
            ancestors[i] <- list(summary$ancestors)
        }
        ids <- if (length(labels)) as.character(labels[[j]][, 1L]) else character()
        map <- if (length(ids)) stats::setNames(lapply(split(orig,
            factor(out$mergedPeak, levels = seq_along(ids))), unique), ids)
        else if (isTRUE(no_sample[i])) NULL else list()
        clusters[i] <- list(map)
    }
    sp$msPurityContributors <- contributors
    sp$msPurityClusterContributors <- clusters
    sp$msPurityAcquisitionValues <- values
    sp$msPurityAcquisitionStatus <- status
    if (any(configured[exact])) {
        sp$msPurityPrecursorPath <- paths
        sp$msPurityAncestorConditions <- ancestors
    }
    sp
}

#' Extract the spectra and recorded inputs behind selected averages
#'
#' Includes the requested averages and earlier averages they selected, the
#' associated measured spectra and filtered inputs, their feature links and
#' immediate membership tables. For disk-backed measured spectra the backend's
#' usual file dependencies apply after serialization.
#' The averaged peak arrays must be the unmodified averaging outputs.
#'
#' @param object A purityA object containing recorded average provenance.
#' @param spectrumId Character vector of `spectrumId` values of the averages to extract.
#' @return A list with `averages` (including their derived ancestors),
#'   `measured`, `fragments`, `featureLinks`, `mergeMembers`, `mergePeaks`,
#'   `lineage` and `processing`. Measured spectra include available MSn
#'   precursor scans along recorded paths. In `lineage`, `parentLineageScanId`
#'   joins to the measured parent's `lineageScanId`, and `parentPid` joins to
#'   its `pid`. The averaged Spectra carries the restricted graph, so
#'   [averageMembers()] and [averagePeakMembers()] also work on its subsets.
#'   The `merged` and derived `member` values match `averages$spectrumId`.
#'   Measured inputs have a `pid`, `fileid` and `linkId`.
#'   Resolve `rank` in `mergePeaks` against `(merged, rank)` in `mergeMembers`.
#' @examples
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'     "2_frag4feature_pa.rds", package = "msPurity"))
#' pa <- averageAllFragSpectra(pa, cores = 1)
#' s <- averagedSpectra(pa, legacy = FALSE)
#' bundle <- averageProvenance(pa, s$spectrumId[1L])
#' head(bundle$mergeMembers)
#' @export
averageProvenance <- function(object, spectrumId) {
    if (!methods::is(object, "purityA")) stop("Expected a purityA object.")
    g <- .pa_graph_data(object)
    sp <- object@avSpectra
    keep <- todo <- unique(.pa_graph_resolve(sp, g, spectrumId, peaks = TRUE))
    while (length(todo)) {
        rows <- match(todo, g$nodes$key)
        if (anyNA(rows) || !all(g$nodes$complete[rows]))
            stop("The selected spectra have unknown historical provenance.")
        m <- g$mergeMembers
        parents <- m$member[m$merged %in% todo & is.na(m$linkId)]
        todo <- setdiff(parents, keep)
        keep <- c(keep, todo)
    }
    g$nodes <- g$nodes[g$nodes$key %in% keep, , drop = FALSE]
    for (nm in c("mergeMembers", "mergePeaks")) {
        tab <- g[[nm]]
        g[[nm]] <- tab[tab$merged %in% keep, , drop = FALSE]
    }
    m <- g$mergeMembers
    g$scans <- g$scans[g$scans$key %in% m$member[!is.na(m$linkId)], , drop = FALSE]
    g$links <- g$links[g$links$linkId %in% m$linkId, , drop = FALSE]
    g <- .pa_graph_index(g)
    sp <- .pa_graph_view(object, sp[.pa_graph_sp_keys(sp, g) %in% keep])
    meta <- .pa_spectra_metadata(sp)
    meta$msPurityAverageGraph <- g
    .pa_spectra_metadata(sp) <- meta
    pid <- g$scans$pid
    lineage <- object@params$msnLineage$scans
    if (!is.null(lineage)) {
        paths <- lineage$path[lineage$pid %in% pid]
        pid <- union(pid, unlist(paths, use.names = FALSE))
        lineage <- lineage[lineage$pid %in% pid, , drop = FALSE]
    }
    links <- g$links
    live <- object@grped_df
    extra <- setdiff(names(live), names(links))
    if (length(extra) && NROW(links) > 0L)
        links[extra] <- live[match(links$linkId, live$linkId), extra, drop = FALSE]
    peaks <- g$mergePeaks
    peaks$merged <- g$nodes$spectrumId[match(peaks$merged, g$nodes$key)]
    list(averages = sp, measured = object@spectra[object@spectra$pid %in% pid],
        fragments = object@fragSpectra[object@fragSpectra$pid %in% g$scans$pid],
        featureLinks = links, mergeMembers = .pa_graph_member_ids(g, m), mergePeaks = peaks,
        lineage = lineage, processing = processingHistory(sp))
}
