#' @rdname as_igraph
#' @param labels Optional character vector in node order. A single `NA`
#'   hides every label. Scan labels include the precursor m/z and precursor
#'   purity when present. Peak labels give m/z and peak position. Labels
#'   are hidden above 30 vertices. All vertices are drawn.
#' @param ... Additional igraph plotting arguments, such as `layout`,
#'   `vertex.size`, `vertex.label.cex` or `main`.
#' @details Square vertices identify nodes marked as context, reference only
#'   or unmeasured in the input. A dashed arrow marks a recorded reference
#'   or a context edge. The plot uses a neutral grey palette.
#' @export
plot.MassSpectralTrees <- function(x, ..., structure = NULL, labels = NULL) {
    graph <- as_igraph(x, structure)
    type <- igraph::graph_attr(graph, "structure")
    title <- if (type == "MST") "Mass spectral tree" else "Precursor-fragment tree"
    options <- list(...)
    count <- igraph::vcount(graph)
    if (length(labels) == 1L && is.atomic(labels) && is.na(labels))
        labels <- rep(NA_character_, count)
    if (!is.null(labels) && (!is.character(labels) || length(labels) != count))
        stop("Supply one label per node or NA to hide labels.", call. = FALSE)
    if (!count) {
        graphics::plot.new()
        if (is.null(options$main)) options$main <- paste(title, "(empty)")
        graphics::title(main = options$main)
        return(invisible(graph))
    }

    nodes <- igraph::vertex_attr(graph)
    edges <- igraph::edge_attr(graph)
    context <- rep(FALSE, count)
    if (!is.null(nodes$context)) context <- nodes$context %in% TRUE
    if (!is.null(nodes$referenceOnly))
        context <- context | nodes$referenceOnly %in% TRUE
    if (!is.null(nodes$measured)) context <- context | nodes$measured %in% FALSE
    dashed <- edges$relationship == "recorded_reference"
    if (!is.null(edges$context)) dashed <- dashed | edges$context %in% TRUE
    if (is.null(labels)) labels <- .tree_plot_labels(nodes, type)
    if (is.null(options$layout)) {
        if (igraph::ecount(graph)) {
            layout <- igraph::layout_as_tree(graph, mode = "out")
            if (type == "PFT") layout <- cbind(-layout[, 2L], layout[, 1L])
            options$layout <- layout
        } else {
            options$layout <- igraph::layout_in_circle(graph)
        }
    }
    # Leave space above precursor nodes for labels and to the right of products.
    above <- type == "PFT" & igraph::degree(graph, mode = "out") > 0
    defaults <- list(main = title, vertex.label = labels,
        vertex.shape = ifelse(context, "square", "circle"),
        vertex.color = ifelse(context, "grey80", "grey97"),
        vertex.frame.color = "grey35", vertex.size = 8,
        vertex.label.color = "grey15", vertex.label.family = "sans",
        vertex.label.cex = 0.8, vertex.label.adj = c(0, 0.5),
        vertex.label.degree = ifelse(above, -0.8 * pi, 0),
        vertex.label.dist = ifelse(above, 5, 1.6),
        edge.color = "grey50", edge.arrow.size = 0.35,
        edge.lty = ifelse(dashed, 2L, 1L),
        margin = c(0.17, 0.24, 0.17, 0.64), asp = 0)
    args <- utils::modifyList(defaults, options, keep.null = TRUE)
    do.call(graphics::plot, c(list(x = graph), args))
    invisible(graph)
}

#' @rdname as_igraph
#' @export
plot.PrecursorFragmentTrees <- plot.MassSpectralTrees

#' @rdname as_igraph
#' @export
plot.SpectralTree <- plot.MassSpectralTrees

.tree_plot_labels <- function(nodes, type) {
    if (length(nodes$name) > 30L) return(NA_character_)
    labels <- nodes$name
    if (!is.null(nodes$label)) {
        known <- !is.na(nodes$label)
        labels[known] <- nodes$label[known]
    }
    if (type == "PFT") {
        if (!is.numeric(nodes$mz)) return(labels)
        known <- is.finite(nodes$mz)
        labels[known] <- sprintf("m/z %.5f", nodes$mz[known])
        if (!is.null(nodes$peak_index)) {
            peaks <- known & !is.na(nodes$peak_index)
            labels[peaks] <- paste0(labels[peaks], "\npeak ", nodes$peak_index[peaks])
        }
        return(labels)
    }
    if (!is.null(nodes$spectrumId)) {
        known <- !is.na(nodes$spectrumId)
        labels[known] <- nodes$spectrumId[known]
    }
    if (!is.null(nodes$acquisitionNum)) {
        known <- !is.na(nodes$acquisitionNum)
        labels[known] <- paste("Scan", nodes$acquisitionNum[known])
        if (!is.null(nodes$msLevel)) {
            known <- known & !is.na(nodes$msLevel)
            labels[known] <- paste0("MS", nodes$msLevel[known], "  ", labels[known])
        }
    } else {
        labels <- sub(".*\\[acquisition ([0-9]+)\\].*", "Scan \\1", labels)
    }
    if (is.numeric(nodes$precursorMz)) {
        known <- is.finite(nodes$precursorMz)
        labels[known] <- paste0(labels[known],
            sprintf("\nm/z %.5f", nodes$precursorMz[known]))
    }
    if (is.numeric(nodes$inPurity)) {
        known <- is.finite(nodes$inPurity)
        labels[known] <- paste0(labels[known],
            sprintf("\npurity %.1f%%", 100 * nodes$inPurity[known]))
    }
    labels
}
