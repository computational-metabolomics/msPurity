#' @rdname as_igraph
#' @param file GraphML output path in an existing directory.
#' @param overwrite Whether to replace an existing file after export succeeds.
#' @param max_nodes,max_edges Maximum numbers of vertices and edges to export.
#'   Increase these limits for larger trees.
#' @details GraphML stores scalar text, logical and numeric attributes.
#'   Factors are exported as text. An attribute with missing values has a
#'   companion logical flag. The graph attribute
#'   \code{SpectralTrees_missing_json} maps attributes to their flags for
#'   vertices and edges. \code{SpectralTrees_metadata_json} contains the
#'   selected tree's metadata. Text and attribute names are written as UTF-8.
#' @export
writeGraphML <- function(x, file, structure = NULL, overwrite = FALSE,
        max_nodes = 2000L, max_edges = 5000L) {
    view <- .tree_display_view(x, structure)
    .tree_check_limits(max_nodes, max_edges)
    if (is.data.frame(view$nodes) && nrow(view$nodes) > max_nodes ||
            is.data.frame(view$edges) && nrow(view$edges) > max_edges)
        stop("The tree exceeds max_nodes or max_edges.", call. = FALSE)
    if (!is.character(file) || length(file) != 1L || is.na(file) ||
            !nzchar(file) || !is.logical(overwrite) || length(overwrite) != 1L ||
            is.na(overwrite))
        stop("Supply one output path and a logical overwrite value.", call. = FALSE)
    file <- path.expand(file)
    if (!dir.exists(dirname(file)) || dir.exists(file))
        stop("Supply a file in an existing directory.", call. = FALSE)
    file <- file.path(normalizePath(dirname(file)), basename(file))
    if (!overwrite && file.exists(file))
        stop("Output exists. Use overwrite = TRUE to replace it.", call. = FALSE)

    graph <- as_igraph(view)
    vertices <- .tree_graphml_attributes(igraph::vertex_attr(graph))
    edges <- .tree_graphml_attributes(igraph::edge_attr(graph))
    igraph::vertex_attr(graph) <- vertices$values
    igraph::edge_attr(graph) <- edges$values
    missing <- list(vertex = vertices$missing, edge = edges$missing)
    encode <- function(x) as.character(jsonlite::toJSON(x, auto_unbox = TRUE,
        null = "null", na = "null", digits = NA))
    graph <- igraph::set_graph_attr(graph, "SpectralTrees_metadata_json",
        encode(view$metadata))
    graph <- igraph::set_graph_attr(graph, "SpectralTrees_missing_json",
        encode(missing))
    graph <- igraph::set_graph_attr(graph, "SpectralTrees_export_version", "1")

    tmp <- tempfile(".tree-graphml-", dirname(file))
    on.exit(unlink(tmp), add = TRUE)
    igraph::write_graph(graph, tmp, format = "graphml")
    # A hard link makes the no-overwrite check atomic even for concurrent writers.
    published <- if (overwrite) file.rename(tmp, file) else file.link(tmp, file)
    if (!published) stop("Could not publish GraphML at ", file, call. = FALSE)
    invisible(file)
}

.tree_graphml_text <- function(value) {
    text <- iconv(value, to = "UTF-8", sub = NA_character_)
    controls <- "[\\x01-\\x08\\x0B\\x0C\\x0E-\\x1F\uFFFE\uFFFF]"
    if (anyNA(text) || any(grepl(controls, text, perl = TRUE)))
        stop("A GraphML attribute contains invalid XML 1.0 text.", call. = FALSE)
    text
}

.tree_graphml_attributes <- function(attributes) {
    names(attributes) <- .tree_graphml_text(names(attributes))
    missing <- list()
    for (name in names(attributes)) {
        value <- attributes[[name]]
        if (is.factor(value)) value <- as.character(value)
        if (is.object(value) || !is.null(dim(value)) ||
                !typeof(value) %in% c("character", "integer", "double", "logical"))
            stop("GraphML attribute '", name,
                "' must contain scalar text, logical or numeric values.", call. = FALSE)
        if (is.numeric(value) && any(is.infinite(value)))
            stop("GraphML attribute '", name, "' contains infinite values.", call. = FALSE)
        absent <- is.na(value)
        if (any(absent)) {
            flag <- paste0(name, "_missing")
            while (flag %in% names(attributes)) flag <- paste0(flag, "_")
            attributes[[flag]] <- absent
            missing[[name]] <- flag
            value[absent] <- switch(typeof(value), character = "", logical = FALSE, 0)
        }
        if (is.character(value)) value <- .tree_graphml_text(value)
        attributes[[name]] <- value
    }
    list(values = attributes, missing = missing)
}
