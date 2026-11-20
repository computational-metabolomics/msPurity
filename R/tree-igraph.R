#' Convert, plot and export spectral trees
#'
#' Use the node and edge tables of a scan or precursor-fragment tree to
#' create a directed graph, draw it or export it as GraphML.
#'
#' @param x A `MassSpectralTrees` result from [massSpectralTrees()], a
#'   `PrecursorFragmentTrees` result from [precursorFragmentTrees()], or its
#'   scan or peak view.
#' @param structure Optional `"MST"` or `"PFT"`. A scan tree defaults to
#'   `"MST"` and a precursor-fragment result defaults to `"PFT"`. Select
#'   `"MST"` on a precursor-fragment result to use its stored scan tree.
#'
#' @details Vertex names are the scan lineage IDs or the stored node IDs of
#'   the selected view. Each peak occurrence has its own vertex. All node
#'   and edge columns are retained as attributes, including the endpoints and
#'   relationship types. Edges follow the direction recorded in the tables.
#'   Isolated nodes and empty results are retained.
#'
#'   A node annotation named `name` is stored as `annotation_name` because
#'   igraph reserves `name` for vertex identities. Having both annotation
#'   columns is an error.
#'
#' @return `as_igraph()` returns a directed igraph with graph attribute
#'   `structure` set to `"MST"` or `"PFT"`. The plot methods return the
#'   displayed graph invisibly. `writeGraphML()` returns the output path
#'   invisibly.
#' @order -1
#' @export
#' @examples
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'     data <- S4Vectors::DataFrame(
#'         spectrumId = c("scan=2", "scan=3"),
#'         dataOrigin = "example.mzML", msLevel = c(2L, 3L),
#'         acquisitionNum = c(2L, 3L), scanIndex = c(1L, 2L),
#'         precScanNum = c(NA_integer_, 2L), rtime = c(1, 2),
#'         precursorMz = c(500, 150))
#'     data$mz <- list(150, 75)
#'     data$intensity <- list(10, 5)
#'     spectra <- Spectra::Spectra(data, source = Spectra::MsBackendMemory())
#'     trees <- massSpectralTrees(spectra, parentReference = "acquisitionNum")
#'     graph <- as_igraph(trees)
#'     igraph::as_data_frame(graph, what = "edges")
#'     plot(trees)
#'     path <- tempfile(fileext = ".graphml")
#'     writeGraphML(trees, path)
#'     restored <- igraph::read_graph(path, format = "graphml")
#'     igraph::V(restored)$lineageScanId
#'     unlink(path)
#' }
as_igraph <- function(x, structure = NULL) {
    view <- .tree_display_view(x, structure)
    if (!requireNamespace("igraph", quietly = TRUE))
        stop("Install the optional package 'igraph' to convert tree tables.",
            call. = FALSE)

    .tree_table_graph(view)
}

.tree_display_view <- function(x, structure) {
    peak_tree <- inherits(x, c("PrecursorFragmentTrees", "PFT"))
    if (!peak_tree && !inherits(x, c("MassSpectralTrees", "MST")))
        stop("Supply a MassSpectralTrees or PrecursorFragmentTrees result ",
            "or its MST or PFT view.", call. = FALSE)
    primary <- if (peak_tree) "PFT" else "MST"
    type <- if (is.null(structure)) primary else
        match.arg(structure, c("MST", "PFT"))
    stored <- x[[tolower(type), exact = TRUE]]
    if (type == "MST" && inherits(stored, "MassSpectralTrees")) return(stored)
    if (!is.null(stored) && inherits(stored, type)) return(stored)
    if (type != primary || inherits(x, "PrecursorFragmentTrees"))
        stop("The requested ", type, " view is absent from this tree.",
            call. = FALSE)
    x
}

.tree_table_index <- function(x) {
    nodes <- x$nodes
    edges <- x$edges
    lineage <- inherits(x, "MassSpectralTrees")
    key <- if (lineage) "lineageScanId" else "node_id"
    endpoints <- if (lineage) c("parentLineageScanId", "lineageScanId") else
        c("source", "target")
    if (!is.data.frame(nodes) || !is.data.frame(edges) ||
            !key %in% names(nodes) ||
            !all(c(endpoints, "relationship") %in% names(edges)) ||
            anyDuplicated(names(nodes)) || anyDuplicated(names(edges)))
        stop("Tree tables must contain unique columns and node identifiers.",
            call. = FALSE)
    ids <- nodes[[key]]
    if (!is.character(ids) || anyNA(ids) || any(!nzchar(ids)) ||
            anyDuplicated(ids))
        stop("Node ", key, " values must be unique, nonmissing text.",
            call. = FALSE)
    rows <- lapply(edges[endpoints], function(values) {
        if (!is.character(values))
            stop("Edge identifiers must reference present nodes.",
                call. = FALSE)
        index <- match(values, ids)
        if (anyNA(index))
            stop("Edge identifiers must reference present nodes.", call. = FALSE)
        index
    })
    if (!is.character(edges$relationship) || anyNA(edges$relationship) ||
            any(!nzchar(edges$relationship)))
        stop("Every edge must have a nonmissing relationship type.",
            call. = FALSE)
    list(ids = ids, endpoints = endpoints, from = rows[[1L]], to = rows[[2L]])
}

.tree_table_graph <- function(x) {
    index <- .tree_table_index(x)
    nodes <- x$nodes
    edges <- x$edges
    if ("name" %in% names(nodes)) {
        if ("annotation_name" %in% names(nodes))
            stop("Rename the node name or annotation_name column before ",
                "igraph conversion.", call. = FALSE)
        names(nodes)[names(nodes) == "name"] <- "annotation_name"
    }
    vertices <- data.frame(.vertex_id = index$ids, nodes, check.names = FALSE)
    relationships <- data.frame(.from = edges[[index$endpoints[1L]]],
        .to = edges[[index$endpoints[2L]]], edges, check.names = FALSE)
    graph <- igraph::graph_from_data_frame(relationships,
        directed = TRUE, vertices = vertices)
    igraph::set_graph_attr(graph, "structure", if (inherits(x, "PFT")) "PFT" else "MST")
}
