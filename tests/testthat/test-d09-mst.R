.mst_two_runs <- function() {
    data <- S4Vectors::DataFrame(
        lineageScanId = c("a1", "a2", "a3", "a4", "b1", "b2"),
        dataOrigin = c(rep("first.mzML", 4L), rep("second.mzML", 2L)),
        msLevel = c(1:4, 1:2),
        spectrumId = c("s1", "s2", "s3", "s4", "s1", "s2"),
        acquisitionNum = c(11:14, 11:12),
        scanIndex = c(1:4, 1:2),
        precScanNum = c(NA_integer_, 11:13, NA_integer_, 11L),
        rtime = c(1:4, 1:2),
        precursorMz = c(NA, 300, 200, 100, NA, 350))
    data$mz <- rep(list(c(50, 100)), nrow(data))
    data$intensity <- rep(list(c(10, 20)), nrow(data))
    Spectra::Spectra(data, source = Spectra::MsBackendMemory())
}

test_that("MSTs resolve MS1 to MS4 within a run using metadata alone", {
    spectra <- Spectra::addProcessing(.mst_two_runs(), function(x, ...) {
        stop("Queued peak processing was executed.")
    })
    tree <- massSpectralTrees(spectra)

    expect_s3_class(tree, "MassSpectralTrees")
    expect_identical(tree$nodes$lineageScanId,
        c("a1", "a2", "a3", "a4", "b1", "b2"))
    expect_identical(tree$nodes$parentStatus,
        c("root", "valid", "valid", "valid", "root", "valid"))
    expect_identical(tree$edges$parentLineageScanId, c("a1", "a2", "a3", "b1"))
    expect_identical(tree$edges$lineageScanId, c("a2", "a3", "a4", "b2"))
})

test_that("native spectrum references resolve after reordering or subsetting", {
    spectra <- .mst_two_runs()
    spectra$lineageScanId <- NULL
    spectra$precScanNum <- rep(NA_integer_, 6L)
    spectra$spectrumRef <- c(NA, "s1", "s2", "s3", NA, "s1")
    ids <- massSpectralTrees(spectra)$nodes$lineageScanId

    order <- c(6L, 2L, 5L, 4L, 1L, 3L)
    reordered <- massSpectralTrees(spectra[order])
    expect_identical(anyDuplicated(ids), 0L)
    expect_identical(reordered$nodes$lineageScanId, ids[order])
    expect_identical(reordered$edges$parentLineageScanId, ids[c(5, 1, 3, 2)])
    expect_identical(reordered$edges$lineageScanId, ids[c(6, 2, 4, 3)])

    # The first file's MS1 scan cannot supply the second file's missing parent.
    partial <- massSpectralTrees(spectra[c(1, 6)])
    expect_identical(partial$nodes$lineageScanId, ids[c(1, 6)])
    expect_identical(partial$nodes$parentStatus, c("root", "unresolved_reference"))
    expect_identical(nrow(partial$edges), 0L)
})

test_that("acquisition numbers supply stable identities when native keys are absent", {
    spectra <- .mst_two_runs()
    spectra$lineageScanId <- NULL
    spectra$spectrumId <- rep(NA_character_, 6L)
    spectra$scanIndex <- rep(NA_integer_, 6L)

    complete <- massSpectralTrees(spectra)
    ids <- complete$nodes$lineageScanId
    second_run <- massSpectralTrees(spectra[5:6])

    expect_identical(anyDuplicated(ids), 0L)
    expect_identical(second_run$nodes$lineageScanId, ids[5:6])
    expect_identical(second_run$edges$parentLineageScanId, ids[5L])
    expect_identical(second_run$edges$lineageScanId, ids[6L])
})

.mst_ambiguous_reference <- function() {
    spectra <- .mst_two_runs()[1:3]
    spectra$msLevel <- c(2L, 2L, 3L)
    spectra$acquisitionNum <- c(2L, 1L, 3L)
    spectra$precScanNum <- c(NA_integer_, NA_integer_, 1L)
    spectra
}

test_that("retained pids preserve D06 lineage when native scan keys are absent", {
    spectra <- .mst_two_runs()[1:2]
    spectra$lineageScanId <- NULL
    spectra$spectrumId <- rep(NA_character_, 2L)
    spectra$scanIndex <- rep(NA_integer_, 2L)
    spectra$pid <- c(30L, 50L)
    pa <- methods::new("purityA", spectra = spectra)

    direct <- massSpectralTrees(spectra)
    processed <- massSpectralTrees(pa)
    expect_identical(direct$nodes$lineageScanId,
        c("first.mzML|pid%3A30", "first.mzML|pid%3A50"))
    expect_identical(direct$edges, processed$edges)
    expect_identical(direct$nodes, processed$nodes)
})

test_that("explicit parent choices survive in purityA while the source is current", {
    spectra <- .mst_ambiguous_reference()
    expect_error(massSpectralTrees(spectra), "Ambiguous parent")
    by_index <- massSpectralTrees(spectra, parentReference = "scanIndex")
    expect_identical(by_index$edges$parentLineageScanId, "a1")
    expect_identical(by_index$edges$lineageScanId, "a3")

    spectra$pid <- c(14L, 18L, 25L)
    pa <- methods::new("purityA", spectra = spectra)
    pa@params$msnLineage <- msPurity:::.pa_msn_preparation(pa, "acquisitionNum")

    tree <- massSpectralTrees(pa)
    expect_identical(tree$edges$parentLineageScanId, "a2")
    expect_identical(tree$edges$lineageScanId, "a3")
    expect_identical(tree$edges$parentReferenceType, "acquisitionNum")

    # Editing the child reference invalidates the retained parent choice.
    pa@spectra$precScanNum <- c(NA_integer_, NA_integer_, 2L)
    expect_error(massSpectralTrees(pa), "no longer matches acquisition metadata")
})

test_that("averages with inherited acquisition metadata are refused", {
    spectra <- .mst_two_runs()[2:3]
    spectra$av_level <- rep("all", 2L)
    expect_error(massSpectralTrees(spectra), "[Aa]verag|[Cc]onsensus|[Mm]easured")
})
