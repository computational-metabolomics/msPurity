.selection_scan_tree <- function() {
    nodes <- data.frame(lineageScanId = c("b1", "a3", "a1", "a4", "a2", "a5"),
        run_id = c("b", rep("a", 5)), acquisitionNum = c(1L, 3L, 1L, 4L, 2L, 5L),
        parentLineageScanId = c(NA, "a2", NA, "a1", "a1", "a3"),
        inPurity = c(NA, 0, NA, 1, 0.83, NA), passed = c(NA, FALSE, NA, TRUE, TRUE, NA))
    edges <- data.frame(parentLineageScanId = c("a3", "a1", "a2", "a1"),
        lineageScanId = c("a5", "a4", "a3", "a2"), relationship = "acquisition")
    structure(list(nodes = nodes, edges = edges, metadata = list(sourceType = "Spectra")),
        class = "MassSpectralTrees")
}

test_that("run-qualified selection retains row order, annotations and isolates", {
    forest <- .selection_scan_tree()
    expect_error(subsetTree(forest, scan = 3L), "multi-run")
    expect_error(subsetTree(forest, run = "other", scan = 3L), "exactly one")
    tree <- subsetTree(forest, run = "a", scan = 3L)
    expect_s3_class(tree, "MassSpectralTrees")
    expect_identical(tree$nodes, forest$nodes[2:6, ])
    expect_identical(tree$edges, forest$edges)
    expect_identical(tree$metadata, forest$metadata)
    isolate <- subsetTree(forest, run = "b", scan = 1L, max_edges = 0L)
    expect_identical(isolate$nodes, forest$nodes[1L, ])
    expect_identical(isolate$edges, forest$edges[FALSE, ])
})

test_that("branch direction, depth and parent context follow stored edges", {
    tree <- .selection_scan_tree()
    out <- subsetTree(tree, node = "a3", mode = "out", depth = 1L,
        includeParent = TRUE)
    expect_identical(out$nodes, tree$nodes[c(2L, 5L, 6L), ])
    expect_identical(out$edges, tree$edges[c(1L, 3L), ])
    up <- subsetTree(tree, node = "a3", mode = "in", depth = 1L)
    expect_identical(up$nodes$lineageScanId, c("a3", "a2"))
    both <- subsetTree(tree, node = "a2", depth = 1L)
    expect_identical(both$nodes$lineageScanId, c("a3", "a1", "a2"))
    only <- subsetTree(out, node = "a2", depth = 0)
    expect_identical(only$nodes$parentLineageScanId, "a1")
    expect_equal(nrow(only$edges), 0L)
    expect_identical(subsetTree(only, node = "a2", mode = "in"), only)
})

test_that("selection budgets apply to the returned nodes and induced edges", {
    tree <- .selection_scan_tree()
    tree$edges <- rbind(tree$edges, tree$edges[3L, ])
    exact <- subsetTree(tree, node = "a3", mode = "out", depth = 0L,
        includeParent = TRUE, max_nodes = 2L, max_edges = 2L)
    expect_equal(nrow(exact$edges), 2L)
    expect_error(subsetTree(tree, node = "a3", mode = "out", depth = 0L,
        includeParent = TRUE, max_nodes = 1L), "max_nodes")
    expect_error(subsetTree(tree, node = "a3", mode = "in", depth = 1L,
        max_edges = 1L), "max_edges")
    expect_error(subsetTree(tree, node = "a3", max_nodes = 3L), "max_nodes")
    expect_error(subsetTree(tree, node = "b1", max_nodes = 0L), "max_nodes")
})

test_that("peak lookup uses the saved run, spectrum and stored position", {
    nodes <- data.frame(node_id = c("p", "q", "r", "s"), run_id = c("a", "a", "b", "a"),
        file = c("f1", "f1", "f2", "f1"), spectrum_id = c("x", "x", "x", "y"),
        peak_index = c(4L, 5L, 4L, 1L), mz = c(150, 150, 150, 98))
    edges <- data.frame(source = "p", target = "s", relationship = "fragmentation",
        observation_id = "event")
    view <- structure(list(nodes = nodes, edges = edges, metadata = list()),
        class = c("PFT", "SpectralTree"))
    spectra <- data.frame(spectrum_id = c("y", "x", "x"), dataOrigin = c("f1", "f2", "f1"),
        acquisitionNum = c(14L, 9L, 12L))
    saved <- structure(list(pft = view, spectrum_metadata = spectra),
        class = "PrecursorFragmentTrees")
    branch <- subsetTree(saved, run = "a", scan = 12L, peak = 4L, mode = "out")
    expect_identical(class(branch), class(view))
    expect_identical(branch$nodes$node_id, c("p", "s"))
    expect_identical(branch$nodes$acquisitionNum, c(12L, 14L))
    expect_identical(branch$edges, edges)
    expect_identical(subsetTree(branch, scan = 14L, peak = 1L), branch)
    expect_identical(subsetTree(saved, run = "b", scan = 9L, peak = 4L)$nodes$node_id, "r")
    expect_error(subsetTree(saved, scan = 12L, peak = 4L), "multi-run")
    expect_error(subsetTree(saved, run = "a", scan = 12L), "peak position")
    expect_error(subsetTree(view, scan = 12L, peak = 4L), "Scan numbers")
    expect_identical(subsetTree(view, node = "q")$nodes$node_id, "q")
    saved$spectrum_metadata <- rbind(spectra, spectra[3L, ])
    expect_error(subsetTree(saved, run = "a", scan = 12L, peak = 4L), "exactly one")
})

test_that("invalid starts and unsupported traversal options fail explicitly", {
    tree <- .selection_scan_tree()
    expect_error(subsetTree(tree, node = "absent"), "absent")
    expect_error(subsetTree(tree, node = "a1", scan = 1L), "alone")
    expect_error(subsetTree(tree, run = "a", scan = 1L, peak = 1L), "only used")
    expect_error(subsetTree(tree, node = "a1", depth = -1), "depth")
    expect_error(subsetTree(tree, node = "a1", includeParent = TRUE), "includeParent")
})

test_that("saved scan metadata resolves identified peaks beside an external reference", {
    peaks <- data.frame(node_id = c("p", "q", "ref"), spectrum_id = c("s1", "s2", NA),
        acquisitionNum = c(8L, NA_integer_, NA_integer_), peak_index = c(3L, 1L, NA))
    edges <- data.frame(source = "p", target = "q", relationship = "fragmentation")
    view <- structure(list(nodes = peaks, edges = edges), class = c("PFT", "SpectralTree"))
    metadata <- data.frame(spectrum_id = c("s2", "s1"), dataOrigin = "file1",
        acquisitionNum = c(9L, 8L))
    saved <- structure(list(pft = view, spectrum_metadata = metadata),
        class = "PrecursorFragmentTrees")
    selected <- subsetTree(saved, run = "file1", scan = 9L, peak = 1L, mode = "in")
    expect_identical(selected$nodes$node_id, c("p", "q"))
    expect_identical(selected$nodes$run_id, rep("file1", 2L))
    expect_identical(selected$nodes$acquisitionNum, c(8L, 9L))
    expect_identical(subsetTree(saved, node = "ref")$nodes$node_id, "ref")
})

test_that("the seven measured product occurrences retain their directed GraphML IDs", {
    skip_if_not_installed("igraph")
    file <- system.file("extdata", "msn-standards", "kynurenine-pft.rds",
        package = "msPurity", mustWork = TRUE)
    tree <- subsetTree(readRDS(file), scan = 985L, peak = 44L,
        mode = "out", depth = 1L)
    parent <- "peak:kynurenine/985:44"
    products <- paste0("peak:kynurenine/987:", 1:7)
    expect_identical(tree$nodes$node_id, c(parent, products))
    expect_identical(tree$nodes$peak_index, c(44L, 1:7))
    path <- tempfile(fileext = ".graphml")
    on.exit(unlink(path))
    writeGraphML(tree, path)
    graph <- igraph::read_graph(path, format = "graphml")
    expect_true(igraph::is_directed(graph))
    expect_identical(igraph::V(graph)$name, c(parent, products))
    expect_identical(igraph::as_edgelist(graph), cbind(rep(parent, 7L), products,
        deparse.level = 0))
})
