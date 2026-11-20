#' List the direct precursor-product links of a peak tree
#'
#' @param x A [precursorFragmentTrees()] result, a saved peak tree, or a
#'   peak view returned by [subsetTree()].
#' @return A data frame in stored fragmentation-edge order. Edge columns,
#'   including exact `source` and `target` IDs, are retained. An edge annotation
#'   with a reserved output name receives an `edge_` prefix. Added columns
#'   give the product run, and each endpoint's acquisition number, original
#'   peak position and measured m/z. `delta_mz` is the signed difference
#'   `precursor_mz - fragment_mz`. It is missing when either endpoint is
#'   unmeasured. For ions retaining the same known charge and adduct, multiply
#'   `delta_mz` by the absolute charge to screen for a neutral loss in Da.
#'   Older views without scan numbers retain their IDs.
#'   Repeated masses and parallel observations each retain their own rows.
#' @examples
#' f <- system.file("extdata", "msn-standards", "kynurenine-pft.rds", package = "msPurity")
#' links <- fragmentLinks(readRDS(f))
#' subset(links, abs(fragment_mz - 118.10124) <= 118.10124 * 5 / 1e6)
#' # Candidate loss for a singly charged series retaining the same adduct.
#' subset(links, abs(delta_mz - 18.010565) <= 0.002)
#' @export
fragmentLinks <- function(x) {
    view <- .tree_display_view(x, "PFT")
    index <- .tree_table_index(view)
    nodes <- .tree_add_scan_numbers(view$nodes, x[["spectrum_metadata"]])
    if (!all(c("mz", "measured") %in% names(nodes)))
        stop("The peak tree needs mz and measured node columns.", call. = FALSE)
    keep <- view$edges$relationship == "fragmentation"
    links <- view$edges[keep, , drop = FALSE]
    derived <- c("run_id", "precursor_scan", "precursor_peak", "precursor_mz",
        "fragment_scan", "fragment_peak", "fragment_mz", "delta_mz")
    for (field in intersect(derived, names(links))) {
        annotation <- paste0("edge_", field)
        while (annotation %in% names(links)) annotation <- paste0("edge_", annotation)
        names(links)[names(links) == field] <- annotation
    }
    from <- index$from[keep]
    to <- index$to[keep]
    run <- nodes[["run_id"]]
    if (is.null(run)) run <- nodes[["file"]]
    if (is.null(run)) run <- nodes[["dataOrigin"]]
    links$run_id <- if (is.null(run)) rep(NA_character_, nrow(links)) else run[to]
    from[!nodes$measured[from] %in% TRUE] <- NA_integer_
    to[!nodes$measured[to] %in% TRUE] <- NA_integer_
    scans <- nodes[["acquisitionNum"]]
    if (is.null(scans)) scans <- rep(NA_real_, nrow(nodes))
    positions <- nodes[["peak_index"]]
    if (is.null(positions)) positions <- rep(NA_integer_, nrow(nodes))
    links$precursor_scan <- scans[from]
    links$precursor_peak <- positions[from]
    links$precursor_mz <- nodes$mz[from]
    links$fragment_scan <- scans[to]
    links$fragment_peak <- positions[to]
    links$fragment_mz <- nodes$mz[to]
    links$delta_mz <- links$precursor_mz - links$fragment_mz
    links
}
