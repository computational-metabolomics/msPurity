.d09_peak_result <- function() {
    peaks <- data.frame(node_id = c("target:a", "a:1", "a:2", "b:1"),
        mz = c(200, 150, 150, 80), peak_index = c(NA, 1L, 2L, 1L),
        measured = c(FALSE, TRUE, TRUE, TRUE), context = c(TRUE, FALSE, FALSE, FALSE))
    products <- data.frame(source = c("target:a", "target:a"),
        target = c("a:1", "a:2"), relationship = "fragmentation",
        observation_id = c("o1", "o2"))
    pft <- structure(list(nodes = peaks, edges = products, metadata = list(ppm = 12)),
        class = c("PFT", "SpectralTree"))
    scans <- data.frame(node_id = c("scan:a", "scan:b"), msLevel = 1:2)
    links <- data.frame(source = "scan:a", target = "scan:b", relationship = "acquisition")
    mst <- structure(list(nodes = scans, edges = links, metadata = list()),
        class = c("MST", "SpectralTree"))
    structure(list(pft = pft, mst = mst), class = "PrecursorFragmentTrees")
}

test_that("plot dispatch draws scan and peak views and returns the displayed graph", {
    skip_if_not_installed("igraph")
    trees <- .d09_peak_result()
    path <- tempfile(fileext = ".pdf")
    on.exit(unlink(path))
    grDevices::pdf(path)
    tryCatch({
        peaks <- withVisible(plot(trees, labels = c("target", "first", "second", "isolated")))
        expect_false(peaks$visible)
        expect_identical(igraph::V(peaks$value)$name, c("target:a", "a:1", "a:2", "b:1"))
        expect_identical(igraph::V(peaks$value)$mz, c(200, 150, 150, 80))
        expect_identical(igraph::E(peaks$value)$observation_id, c("o1", "o2"))
        expect_identical(igraph::graph_attr(peaks$value, "structure"), "PFT")
        scans <- plot(trees, structure = "MST")
        expect_identical(igraph::graph_attr(scans, "structure"), "MST")
        expect_identical(igraph::as_edgelist(scans), matrix(c("scan:a", "scan:b"), nrow = 1L))
        expect_error(as_igraph(trees$mst, "PFT"), "absent")
        expect_identical(igraph::E(plot(trees$pft, labels = NA_character_))$relationship,
            rep("fragmentation", 2L))
        expect_error(plot(trees, labels = "single"), "one label per node")
    }, finally = grDevices::dev.off())
    expect_gt(file.size(path), 0)
})

test_that("default scan labels display available child scores, including zero purity", {
    nodes <- list(name = c("s1", "s2", "s3"), acquisitionNum = c(966L, 977L, 981L),
        msLevel = 1:3, precursorMz = c(NA, 209.09215, 192.06577), inPurity = c(NA, 1, 0))
    expect_identical(.tree_plot_labels(nodes, "MST"), c("MS1  Scan 966",
        "MS2  Scan 977\nm/z 209.09215\npurity 100.0%",
        "MS3  Scan 981\nm/z 192.06577\npurity 0.0%"))
})

test_that("isolated scans survive conversion and empty trees retain the GraphML schema", {
    skip_if_not_installed("igraph")
    data <- S4Vectors::DataFrame(spectrumId = "s1", dataOrigin = "run",
        acquisitionNum = 1L, msLevel = 1L, mz = list(150), intensity = list(20))
    spectra <- Spectra::Spectra(data, source = Spectra::MsBackendMemory())
    survey <- as_igraph(massSpectralTrees(spectra))
    expect_identical(igraph::V(survey)$name, "run|s1")
    expect_equal(igraph::ecount(survey), 0L)

    tree <- massSpectralTrees(spectra[integer()])
    pdf <- tempfile(fileext = ".pdf")
    path <- tempfile(fileext = ".graphml")
    on.exit(unlink(c(pdf, path)))
    grDevices::pdf(pdf)
    tryCatch(expect_equal(igraph::vcount(plot(tree)), 0L), finally = grDevices::dev.off())
    writeGraphML(tree, path, max_nodes = 0L, max_edges = 0L)
    restored <- igraph::read_graph(path, format = "graphml")
    expect_equal(igraph::vcount(restored), 0L)
    expect_equal(igraph::ecount(restored), 0L)
    expect_true(all(names(tree$nodes) %in% igraph::vertex_attr_names(restored)))
    expect_true(all(names(tree$edges) %in% igraph::edge_attr_names(restored)))
})

test_that("GraphML round trips UTF-8 keys, values and typed missing annotations", {
    skip_if_not_installed("igraph")
    trees <- .d09_peak_result()
    latin1 <- iconv("caf\u00e9", from = "UTF-8", to = "latin1")
    trees$pft$nodes[[latin1]] <- c(latin1, "", NA_character_, "<peak & mass>")
    trees$pft$nodes$measured <- c(FALSE, TRUE, NA, TRUE)
    trees$pft$nodes$measured_missing <- "existing annotation"
    trees$pft$edges[[latin1]] <- c(latin1, NA_character_)
    path <- tempfile(fileext = ".graphml")
    on.exit(unlink(path))
    expect_identical(writeGraphML(trees, path), normalizePath(path))
    graph <- igraph::read_graph(path, format = "graphml")
    expect_true(igraph::is_directed(graph))
    expect_identical(igraph::V(graph)$name, c("target:a", "a:1", "a:2", "b:1"))
    expect_identical(igraph::as_edgelist(graph),
        cbind(c("target:a", "target:a"), c("a:1", "a:2")))
    expect_identical(igraph::E(graph)$relationship, rep("fragmentation", 2L))
    expect_identical(igraph::vertex_attr(graph, "caf\u00e9"),
        c("caf\u00e9", "", "", "<peak & mass>"))
    expect_identical(igraph::edge_attr(graph, "caf\u00e9"), c("caf\u00e9", ""))
    flags <- jsonlite::fromJSON(igraph::graph_attr(graph, "SpectralTrees_missing_json"))
    expect_identical(igraph::vertex_attr(graph, flags$vertex[["caf\u00e9"]]),
        c(FALSE, FALSE, TRUE, FALSE))
    expect_identical(igraph::edge_attr(graph, flags$edge[["caf\u00e9"]]), c(FALSE, TRUE))
    expect_identical(igraph::V(graph)$measured_missing, rep("existing annotation", 4L))
    expect_identical(igraph::vertex_attr(graph, flags$vertex$measured), c(FALSE, FALSE, TRUE, FALSE))
    expect_identical(igraph::vertex_attr(graph, flags$vertex$peak_index), c(TRUE, FALSE, FALSE, FALSE))
    expect_equal(jsonlite::fromJSON(igraph::graph_attr(graph, "SpectralTrees_metadata_json"))$ppm, 12)
})

test_that("GraphML limits and failed replacements preserve the destination", {
    skip_if_not_installed("igraph")
    trees <- .d09_peak_result()
    path <- tempfile(fileext = ".graphml")
    on.exit(unlink(path))
    expect_error(writeGraphML(trees, path, max_nodes = 3L), "exceeds")
    expect_error(writeGraphML(trees, path, max_edges = 1L), "exceeds")
    expect_false(file.exists(path))
    writeLines("original", path)
    expect_error(writeGraphML(trees, path), "overwrite")
    trees$pft$nodes$label <- "invalid\u0001"
    expect_error(writeGraphML(trees, path, overwrite = TRUE), "invalid XML")
    expect_identical(readLines(path), "original")
    trees$pft$edges$nested <- I(list(1:2, 1:3))
    trees$pft$nodes$label <- NULL
    expect_error(writeGraphML(trees, path, overwrite = TRUE), "scalar")
    expect_identical(readLines(path), "original")
    writeGraphML(trees, path, structure = "MST", overwrite = TRUE)
    expect_identical(igraph::V(igraph::read_graph(path, format = "graphml"))$name,
        c("scan:a", "scan:b"))
})
