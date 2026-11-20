.d09_peak_scans <- function() {
    d <- S4Vectors::DataFrame(dataOrigin = "runA", lineageScanId = c("a", "b", "c"),
        acquisitionNum = c(10L, 20L, 30L), scanIndex = 1:3, msLevel = 1:3,
        precursorMz = c(NA, 200, 80), precScanNum = c(0L, 10L, 20L),
        inPurity = c(NA, 0.91, 0.85), centroided = TRUE, rtime = 1:3)
    d$mz <- list(c(100, 200), c(60, 80), c(35, 35))
    d$intensity <- list(c(1, 20), c(6, 8), c(2, 4))
    Spectra::Spectra(d, source = Spectra::MsBackendMemory())
}

test_that("purity annotations come from the current source in stored scan order", {
    s <- .d09_peak_scans()
    s$pid <- c(100L, 104L, 120L)
    pa <- methods::new("purityA", spectra = s)
    pa@params$msnLineage <- msPurity:::.pa_msn_preparation(pa)
    pa@spectra <- s[c(3, 1, 2)]
    pa@spectra$inPurity <- c(0.47, NA, 0.6)
    expect_identical(massSpectralTrees(pa)$nodes$inPurity, c(NA_real_, 0.6, 0.47))
    expect_identical(precursorFragmentTrees(pa)$selections$source, c("peak:a:2", "peak:b:2"))
})

test_that("ambiguous, unmatched and external parents stay explicitly unmeasured", {
    s <- .d09_peak_scans()
    s$precursorMz <- c(NA, 99, 70)
    tree <- precursorFragmentTrees(s, ppm = 0, tolerance = 10)
    expect_identical(tree$selections$match_status, c("matched", "ambiguous"))
    links <- fragmentLinks(tree)
    expect_identical(links$precursor_mz, c(100, 100, NA_real_, NA_real_))
    expect_identical(links$delta_mz, c(40, 20, NA_real_, NA_real_))
    expect_identical(fragmentLinks(tree)$precursor_scan, c(10, 10, NA_real_, NA_real_))
    tight <- precursorFragmentTrees(s[2:3], ppm = 0)
    expect_identical(tight$selections$match_status, c("unresolved_reference", "unmatched"))
    expect_identical(tight$selections$rawParentReference, c("10", "20"))
    expect_true(all(is.na(fragmentLinks(tight)$precursor_mz)))
})

test_that("empty peak arrays and survey-only inputs have typed empty links", {
    s <- .d09_peak_scans()
    only_survey <- precursorFragmentTrees(s[1])
    expect_identical(only_survey$mst$nodes$acquisitionNum, 10)
    expect_identical(fragmentLinks(only_survey)$fragment_mz, numeric())
    expect_identical(fragmentLinks(only_survey)$delta_mz, numeric())
    empty_input <- precursorFragmentTrees(s[integer()])
    expect_identical(fragmentLinks(empty_input)$target, character())
    expect_s3_class(empty_input$mst, "MassSpectralTrees")
    expect_identical(empty_input$mst$nodes$lineageScanId, character())
    expect_identical(empty_input$mst$nodes$parentStatus, character())
    expect_identical(empty_input$mst$edges$parentLineageScanId, character())
    expect_identical(empty_input$mst$edges$lineageScanId, character())
    empty <- Spectra::filterMzRange(s, mz = c(999, 1000))
    x <- precursorFragmentTrees(empty)
    expect_identical(x$selections$acquisitionNum, c(20, 30))
    expect_identical(fragmentLinks(x)$precursor_peak, integer())
    expect_identical(nrow(x$mst$nodes), 3L)
})

test_that("filtered products retain original peak positions and precursor context", {
    s <- .d09_peak_scans()
    matrices <- lapply(Spectra::peaksData(s), function(a)
        cbind(a, original_peak_index = c(3L, 9L), pass_flag = c(1, 0)))
    sv <- as.data.frame(Spectra::spectraData(s, columns =
        setdiff(Spectra::spectraVariables(s), Spectra::peaksVariables(s))))
    flagged <- msPurity:::.msp_memory(sv, matrices, colnames(matrices[[1]]))
    tree <- precursorFragmentTrees(flagged, passingOnly = TRUE)
    links <- fragmentLinks(tree)
    expect_identical(links$source, c("peak:a:9", "peak:b:9"))
    expect_identical(links$target, c("peak:b:3", "peak:c:3"))
    expect_identical(tree$pft$nodes$pass_flag[tree$pft$nodes$context], c(0, 0))
    queued <- Spectra::filterMzRange(s, mz = c(75, 250))
    expect_identical(fragmentLinks(precursorFragmentTrees(queued))$target, "peak:b:2")
    expect_error(precursorFragmentTrees(s, passingOnly = TRUE), "pass_flag")
})

test_that("the precursor window is configurable and requires a supported measurement", {
    s <- .d09_peak_scans()
    # The first target is 10 ppm above its measured parent peak.
    s$precursorMz <- c(NA, 200.002, 80)
    expect_identical(precursorFragmentTrees(s)$selections$match_status, c("unmatched", "matched"))
    expect_identical(precursorFragmentTrees(s, ppm = 20)$selections$match_status, c("matched", "matched"))
    expect_error(precursorFragmentTrees(s, tolerance = -1), "nonnegative")
    s$nSelectedIons <- c(0L, 1L, 2L)
    expect_error(precursorFragmentTrees(s), "Multiple isolated")
    s$nSelectedIons <- NULL
    s$centroided <- c(FALSE, TRUE, TRUE)
    expect_error(precursorFragmentTrees(s), "centroided")
})

test_that("native four-scan peak assignments agree with the independent saved tree", {
    f <- system.file("extdata", "msn-standards", package = "msPurity")
    s <- Spectra::Spectra(file.path(f, "msPurity_L_kynurenine_161166_MSn.mzML.gz"),
        source = Spectra::MsBackendMzR(), BPPARAM = BiocParallel::SerialParam())
    x <- precursorFragmentTrees(s[match(c(966, 977, 985, 987), Spectra::acquisitionNum(s))])
    expect_identical(x$selections$parent_peak_index, c(6L, 32L, 44L))
    key <- function(v) {
        v <- v[order(v$fragment_scan, v$fragment_peak),
            c("precursor_scan", "precursor_peak", "precursor_mz", "fragment_scan", "fragment_peak", "fragment_mz")]
        rownames(v) <- NULL
        v
    }
    saved <- readRDS(file.path(f, "kynurenine-pft.rds"))
    expect_equal(key(fragmentLinks(x)), key(fragmentLinks(saved)), tolerance = 0)
    links <- fragmentLinks(x)
    hit <- subset(links, abs(delta_mz - 18.010565) <= 0.002)
    expect_identical(hit$precursor_scan, 966)
    expect_identical(hit$precursor_peak, 6L)
    expect_identical(hit$fragment_scan, 977)
    expect_identical(hit$fragment_peak, 33L)
    expect_equal(hit$delta_mz, 18.01051331, tolerance = 1e-8)
    expect_true(any(links$delta_mz < 0, na.rm = TRUE))
    expect_identical(nrow(x$pft$nodes), 122L)
    products <- subsetTree(x, scan = 985L, peak = 44L, mode = "out", depth = 1L)
    expect_identical(fragmentLinks(products)$fragment_peak, 1:7)
})

test_that("duplicate scan numbers in different runs have distinct peak occurrences", {
    first <- .d09_peak_scans()
    first$lineageScanId <- NULL
    first$spectrumId <- paste0("scan=", c(10, 20, 30))
    second <- first
    second$dataOrigin <- rep("runB", 3L)
    trees <- precursorFragmentTrees(c(first, second))
    expect_error(subsetTree(trees, scan = 20L, peak = 2L), "multi-run")
    branch <- subsetTree(trees, run = "runB", scan = 20L, peak = 2L, mode = "out")
    links <- fragmentLinks(branch)
    expect_identical(links$precursor_mz, c(80, 80))
    expect_identical(links$fragment_mz, c(35, 35))
    expect_identical(links$fragment_scan, c(30, 30))
    expect_identical(links$observation_id,
        c("observation:runB|scan%3D30:1", "observation:runB|scan%3D30:2"))
    expect_identical(links$run_id, c("runB", "runB"))
    expect_identical(links$target, c("peak:runB|scan%3D30:1", "peak:runB|scan%3D30:2"))
    branch$edges$precursor_mz <- c(79.99, 79.99)
    branch$edges$delta_mz <- c(45.01, 44.99)
    annotated <- fragmentLinks(branch)
    expect_identical(annotated$edge_delta_mz, c(45.01, 44.99))
    expect_identical(annotated$delta_mz, c(45, 45))
    expect_identical(annotated$edge_precursor_mz, c(79.99, 79.99))
    expect_identical(annotated$precursor_mz, c(80, 80))
})

test_that("a stored scan tree with native references reproduces the peak view", {
    s <- .d09_peak_scans()
    s$spectrumId <- paste0("scan=", c(10, 20, 30))
    s$spectrumRef <- c(NA, "scan=10", "scan=20")
    s$precScanNum <- rep(NA_integer_, 3L)
    mst <- massSpectralTrees(s)
    direct <- precursorFragmentTrees(s)
    reused <- precursorFragmentTrees(mst, spectra = s[c(3, 1, 2)])
    expect_identical(reused, direct)
    expect_identical(reused$selections$parent_peak_index, c(2L, 2L))
    scan <- subsetTree(reused, scan = 20L, depth = 0, structure = "MST")
    expect_identical(scan$nodes$inPurity, 0.91)
    expect_error(fragmentLinks(mst), "PFT")
    s$precursorMz <- c(NA, 210, 80)
    expect_error(precursorFragmentTrees(mst, spectra = s), "metadata disagree")
})
