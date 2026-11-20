.d09_graph_fixture <- function() {
    nodes <- data.frame(
        lineageScanId = c("run-a:2", "run-a:3", "run-a:4", "run-b:2"),
        spectrumId = c("scan=2", "scan=3", "scan=4", "scan=2"),
        sourceSpectrumId = c("source-2", "source-3", "source-4", NA_character_),
        msLevel = c(2L, 3L, 4L, 2L), rtime = c(10, 11, 12, 20),
        run_id = c("run-a", "run-a", "run-a", "run-b"),
        parentStatus = c("root", "valid", "valid", "unresolved_reference"),
        stringsAsFactors = FALSE)
    edges <- data.frame(
        parentLineageScanId = c("run-a:2", "run-a:3"),
        lineageScanId = c("run-a:3", "run-a:4"),
        relationship = "acquisition", parentReferenceType = "acquisitionNum",
        rawParentReference = c("2", "3"), parentStatus = "valid",
        stringsAsFactors = FALSE)
    structure(list(nodes = nodes, edges = edges, metadata = list()),
        class = "MassSpectralTrees")
}

test_that("igraph conversion reports its unavailable optional dependency", {
    convert <- as_igraph
    environment(convert) <- list2env(list(
        requireNamespace = function(package, quietly) FALSE),
        parent = environment(convert))

    expect_error(convert(.d09_graph_fixture()), "optional package 'igraph'")
})

test_that("igraph retains measured scan identities, directed edges and attributes", {
    skip_if_not_installed("igraph")
    data <- S4Vectors::DataFrame(
        pid = c(7L, 40L, 41L, 60L), lineageScanId = c("a2", "a3", "a4", "b2"),
        spectrumId = c("scan=2", "scan=3", "scan=4", "scan=2"),
        dataOrigin = c(rep("a.mzML", 3), "b.mzML"), msLevel = c(2:4, 2L),
        acquisitionNum = c(2:4, 2L), scanIndex = c(1:3, 1L),
        precScanNum = c(NA_integer_, 2L, 3L, NA_integer_),
        rtime = c(1:3, 1), precursorMz = c(500, 150, 75, 500))
    spectra <- Spectra::Spectra(data, source = Spectra::MsBackendMemory())
    tree <- massSpectralTrees(spectra, parentReference = "acquisitionNum")
    graph <- as_igraph(tree)

    expect_true(igraph::is_directed(graph))
    expect_identical(igraph::graph_attr(graph, "structure"), "MST")
    expect_identical(igraph::V(graph)$name, c("a2", "a3", "a4", "b2"))
    expect_identical(igraph::as_edgelist(graph), cbind(c("a2", "a3"), c("a3", "a4")))
    expect_identical(igraph::V(graph)$pid, c(7, 40, 41, 60))
    expect_identical(igraph::vertex_attr(graph)[names(tree$nodes)], as.list(tree$nodes))
    expect_identical(igraph::edge_attr(graph)[names(tree$edges)], as.list(tree$edges))
})

test_that("igraph conversion rejects missing and duplicated lineage keys", {
    skip_if_not_installed("igraph")
    trees <- .d09_graph_fixture()
    duplicate <- trees
    duplicate$nodes$lineageScanId[2L] <- duplicate$nodes$lineageScanId[1L]
    expect_error(as_igraph(duplicate), "must be unique")
    missing <- trees
    missing$nodes$lineageScanId[1L] <- NA_character_
    expect_error(as_igraph(missing), "must be unique, nonmissing")
    dangling <- trees
    dangling$edges$parentLineageScanId[1L] <- "absent"
    expect_error(as_igraph(dangling), "reference present nodes")
    untyped <- trees
    untyped$edges$relationship[1L] <- NA_character_
    expect_error(as_igraph(untyped), "nonmissing relationship type")
    expect_error(as_igraph(list()), "MassSpectralTrees")
})

test_that("igraph vertex names cannot replace lineage identity", {
    skip_if_not_installed("igraph")
    trees <- .d09_graph_fixture()
    trees$nodes$name <- c("parent", "intermediate", "child", "isolated")

    graph <- as_igraph(trees)

    expect_identical(igraph::V(graph)$name,
        c("run-a:2", "run-a:3", "run-a:4", "run-b:2"))
    expect_identical(igraph::V(graph)$annotation_name,
        c("parent", "intermediate", "child", "isolated"))
    trees$nodes$annotation_name <- "existing"
    expect_error(as_igraph(trees), "Rename the node name or annotation_name")
})
