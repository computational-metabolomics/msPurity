#' Select a component or branch of a scan or peak tree
#'
#' Follow the recorded edges of a scan or precursor-fragment tree. Select a
#' scan by acquisition number or a peak occurrence by acquisition number and
#' stored peak position. An exact node ID works for either type of tree.
#'
#' @inheritParams as_igraph
#' @param scan One nonnegative acquisition number.
#' @param run Exact run identifier. Required for scan lookup in a multi-run
#'   view. For ordinary mzML input the run identifier is usually the file path.
#' @param peak Positive, one-based `peak_index` in the prepared peak tree's
#'   stored position domain. Required with `scan` on a peak tree.
#' @param node One exact node ID (`lineageScanId` or `node_id`) as an
#'   alternative to `run`, `scan` and `peak`. Use this for reference-only nodes
#'   or views lacking scan numbers or peak positions.
#' @param mode `"all"` follows edges in either direction, `"out"` follows
#'   children and `"in"` follows parents.
#' @param depth Maximum number of edges from the start, or `Inf` for all
#'   reachable nodes. Zero selects the starting node.
#' @param includeParent In `"out"` mode, also keep the starting node's
#'   immediate parents that are present in the view. Their other children
#'   are outside this traversal.
#' @param max_nodes,max_edges Maximum numbers of nodes and edges in the
#'   selection. Exceeding either limit raises an error.
#' @return A `MassSpectralTrees` result for a measured scan tree, or the
#'   selected `MST` or `PFT` view for a prepared result. Node and edge rows
#'   retain their original order and annotations. Edges have both endpoints
#'   in the selection. Parent references on nodes still describe the original
#'   acquisition. Use [as_igraph()], `plot()` or [writeGraphML()] on the result.
#' @details The default selects the whole weakly connected component.
#'   Selection uses the stored tables and works with the raw files offline.
#'   For prepared results, saved spectrum metadata can supply acquisition
#'   numbers. These are carried onto the returned nodes for subsequent lookup.
#' @seealso [massSpectralTrees()], [Spectra::filterPrecursorScan()]
#' @examples
#' file <- system.file("extdata", "msn-standards",
#'     "msPurity_L_kynurenine_161166_MSn.mzML.gz", package = "msPurity")
#' spectra <- Spectra::Spectra(file, source = Spectra::MsBackendMzR(),
#'     BPPARAM = BiocParallel::SerialParam())
#' trees <- massSpectralTrees(spectra)
#' branch <- subsetTree(trees, scan = 977L, mode = "out",
#'     includeParent = TRUE)
#' branch$nodes$acquisitionNum
#' example <- system.file("extdata", "msn-standards", "kynurenine-pft.rds",
#'     package = "msPurity")
#' products <- subsetTree(readRDS(example), scan = 985L, peak = 44L,
#'     mode = "out", depth = 1L)
#' products$nodes[, c("acquisitionNum", "peak_index", "mz")]
#' @export
subsetTree <- function(x, scan = NULL, run = NULL, peak = NULL,
        node = NULL, mode = c("all", "out", "in"), depth = Inf,
        includeParent = FALSE, structure = NULL,
        max_nodes = 2000L, max_edges = 5000L) {
    mode <- match.arg(mode)
    if (!identical(depth, Inf) && !.tree_whole_number(depth))
        stop("depth must be a nonnegative whole number or Inf.", call. = FALSE)
    if (!is.logical(includeParent) || length(includeParent) != 1L ||
            is.na(includeParent) || (includeParent && mode != "out"))
        stop("includeParent must be TRUE or FALSE and requires mode = 'out' when TRUE.",
            call. = FALSE)
    .tree_check_limits(max_nodes, max_edges)
    view <- .tree_display_view(x, structure)
    index <- .tree_table_index(view)
    view$nodes <- .tree_add_scan_numbers(view$nodes, x[["spectrum_metadata"]])
    start <- .tree_start_row(view, index$ids, node, run, scan, peak)
    keep <- .tree_reachable(index, start, mode, depth, max_nodes)
    if (includeParent)
        keep[index$from[index$to == start]] <- TRUE
    edges <- keep[index$from] & keep[index$to]
    if (sum(keep) > max_nodes || sum(edges) > max_edges)
        stop("The selected tree exceeds max_nodes or max_edges.", call. = FALSE)
    view$nodes <- view$nodes[keep, , drop = FALSE]
    view$edges <- view$edges[edges, , drop = FALSE]
    view
}

.tree_whole_number <- function(value) {
    is.numeric(value) && is.null(dim(value)) && length(value) == 1L &&
        is.finite(value) && value >= 0 && value == floor(value)
}

.tree_check_limits <- function(max_nodes, max_edges) {
    if (!.tree_whole_number(max_nodes) || !.tree_whole_number(max_edges))
        stop("max_nodes and max_edges must be nonnegative whole numbers.",
            call. = FALSE)
}

# Prepared views may store acquisition numbers only in the enclosing result.
.tree_add_scan_numbers <- function(nodes, spectra) {
    if (!all(c("spectrum_id", "acquisitionNum") %in% names(spectra)) ||
            !"spectrum_id" %in% names(nodes))
        return(nodes)
    ids <- nodes[["spectrum_id"]]
    keys <- spectra[["spectrum_id"]]
    if (!is.null(spectra[["run_id"]]) && !is.null(nodes[["run_id"]])) {
        group <- nodes[["run_id"]]
        groups <- spectra[["run_id"]]
    } else {
        group <- nodes[["file"]]
        if (is.null(group)) group <- nodes[["dataOrigin"]]
        groups <- spectra[["dataOrigin"]]
    }
    known <- !is.na(ids)
    valid <- !is.na(keys)
    if (!is.null(group) && !is.null(groups)) {
        # An integer group prefix keeps arbitrary native IDs distinct.
        levels <- unique(groups)
        known <- known & !is.na(group) & group %in% levels
        valid <- valid & !is.na(groups)
        ids <- paste(match(group, levels), ids, sep = ":")
        keys <- paste(match(groups, levels), keys, sep = ":")
    }
    ambiguous <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
    valid <- valid & !ambiguous
    pos <- match(ids, keys[valid])
    pos[!known] <- NA_integer_
    numbers <- spectra$acquisitionNum[valid][pos]
    existing <- nodes[["acquisitionNum"]]
    if (!is.null(existing)) numbers[!is.na(existing)] <- existing[!is.na(existing)]
    nodes$acquisitionNum <- numbers
    if (!any(c("run_id", "file", "dataOrigin") %in% names(nodes))) {
        runs <- spectra[["run_id"]]
        if (is.null(runs)) runs <- spectra[["dataOrigin"]]
        if (!is.null(runs)) nodes$run_id <- runs[valid][pos]
    }
    nodes
}

.tree_start_row <- function(view, ids, node, run, scan, peak) {
    if (!is.null(node)) {
        if (!is.character(node) || length(node) != 1L || is.na(node) ||
                !is.null(run) || !is.null(scan) || !is.null(peak))
            stop("Supply one exact node ID alone, or use run, scan and peak.",
                call. = FALSE)
        start <- match(node, ids)
        if (is.na(start)) stop("The node ID is absent from this view.", call. = FALSE)
        return(start)
    }
    if (!.tree_whole_number(scan))
        stop("Supply one nonnegative scan acquisition number or an exact node ID.",
            call. = FALSE)
    nodes <- view$nodes
    if (is.null(nodes[["acquisitionNum"]]))
        stop("Scan numbers are unavailable in this view. Use an exact node ID.",
            call. = FALSE)
    runs <- nodes[["run_id"]]
    if (is.null(runs)) runs <- nodes[["dataOrigin"]]
    if (is.null(runs)) runs <- nodes[["file"]]
    if (is.null(runs) || !any(!is.na(runs) & nzchar(runs)))
        stop("Run identities are incomplete. Use an exact node ID.", call. = FALSE)
    if (is.null(run)) {
        if (anyNA(runs) || any(!nzchar(runs)) || length(unique(runs)) > 1L)
            stop("Supply run for scan lookup in a multi-run or partially identified view.",
                call. = FALSE)
    } else if (!is.character(run) || length(run) != 1L || is.na(run) || !nzchar(run))
        stop("run must be one exact run identifier.", call. = FALSE)
    hit <- nodes$acquisitionNum == scan
    if (!is.null(run)) hit <- hit & runs == run
    if (inherits(view, "PFT")) {
        if (!.tree_whole_number(peak) || peak < 1 || is.null(nodes[["peak_index"]]))
            stop("Supply a positive stored peak position or an exact peak node ID.",
                call. = FALSE)
        hit <- hit & nodes$peak_index == peak
    } else if (!is.null(peak))
        stop("peak is only used on a precursor-fragment tree.", call. = FALSE)
    start <- which(hit)
    if (length(start) != 1L)
        stop("The lookup must match exactly one stored node. Check the run, ",
            "scan and peak, or use an exact node ID.", call. = FALSE)
    start
}

.tree_reachable <- function(index, start, mode, depth, max_nodes) {
    seen <- rep(FALSE, length(index$ids))
    seen[start] <- TRUE
    if (depth == 0 || !length(index$from)) return(seen)
    from <- switch(mode, out = index$from, `in` = index$to, all = c(index$from, index$to))
    to <- switch(mode, out = index$to, `in` = index$from, all = c(index$to, index$from))
    neighbours <- split(to, factor(from, levels = seq_along(seen)))
    frontier <- start
    count <- 1L
    distance <- 0L
    while (length(frontier) && distance < depth) {
        next_rows <- unique(unlist(neighbours[frontier], use.names = FALSE))
        frontier <- next_rows[!seen[next_rows]]
        count <- count + length(frontier)
        if (count > max_nodes)
            stop("The selected tree exceeds max_nodes.", call. = FALSE)
        seen[frontier] <- TRUE
        distance <- distance + 1L
    }
    seen
}
