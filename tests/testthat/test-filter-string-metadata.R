test_that("serial filters describe the final activation event", {
    strings <- c(
        "FTMS + p ESI Full ms2 205.0971@cid30.00 [50.00-216.00]",
        "ITMS + p ESI Full ms3 205.0971@hcd43.75 163.5276@cid15.00 [50-180]",
        "FTMS + p ESI Full ms3 205.0971@cid30 163.5276@hcd15 [50-180]",
        "FTMS + p ESI Full ms3 205.0971@etd30 150.1234@cid20 [50-180]",
        paste("ITMS + p ESI Full ms4 205.0971@hcd43.75",
            "163.5276@hcd15 206.4793@cid30 [50-220]"),
        paste("ITMS + p ESI Full ms4 205.0971@cid10",
            "150.1234@cid20 100.1234@hcd30 [30-120]"))
    result <- msPurity:::.pa_filter_string_metadata(strings,
        c(2L, 3L, 3L, 3L, 4L, 4L))
    expect_identical(result$activationType,
        c("CID", "CID", "HCD", "CID", "CID", "HCD"))
    expect_identical(result$scanFilterAnalyzer,
        c("FTMS", "ITMS", "FTMS", "FTMS", "ITMS", "ITMS"))
})

test_that("unsupported filters leave current activation unknown", {
    strings <- c(
        "FTMS Full ms3 205.0971@cid30 150.1234@etd20 [50-180]",
        "FTMS Full ms3 205.0971@cid30 [50-180]",
        "FTMS Full ms2 205.0971@hcd30@cid10 [50-216]",
        "FTMS Full ms2 205.0971@hcd30 supplemental [50-216]",
        "FTMS Full ms2 205.0971@hcdXX [50-216]",
        "FTMS Full ms2 205.0971@hcd30 [50-216]",
        "FTMS Full ms2 205.0971@hcd30 [50-216]",
        "FTMS Full ms [50-216]", NA_character_)
    levels <- c(3L, 3L, 2L, 2L, 2L, 3L, NA_integer_, 1L, 2L)
    result <- msPurity:::.pa_filter_string_metadata(strings, levels)
    expect_true(all(is.na(result$activationType)))
    expect_identical(result$scanFilterAnalyzer,
        c(rep("FTMS", 8L), NA_character_))
})

test_that("scan-filter fallback preserves known acquisition metadata", {
    source <- data.frame(msLevel = c(2L, 3L),
        filterString = c("FTMS Full ms2 205.0971@hcd30 [50-216]",
            "FTMS Full ms3 205.0971@cid10 150.1234@hcd20 [50-180]"),
        activationType = c("CID", NA_character_),
        activationMethod = c("CID", NA_character_),
        scanFilterAnalyzer = c("ITMS", NA_character_),
        massAnalyzer = c("declared analyzer", NA_character_),
        collisionEnergy = c(NA_real_, 25),
        collisionEnergyUnit = c("NCE", NA_character_))
    original <- source
    result <- msPurity:::.pa_filter_metadata(source)
    expect_identical(result$activationType, c("CID", "HCD"))
    expect_identical(result$scanFilterAnalyzer, c("ITMS", "FTMS"))
    expect_identical(result$massAnalyzer, source$massAnalyzer)
    expect_identical(result$collisionEnergy, source$collisionEnergy)
    expect_identical(result$collisionEnergyUnit, source$collisionEnergyUnit)
    expect_identical(source, original)
    source$filterString <- NULL
    expect_identical(msPurity:::.pa_filter_metadata(source)$activationType,
        c("CID", NA_character_))
})

test_that("D02 MS2 metadata are stored and reused without XML preparation", {
    path <- system.file("extdata", "msn-standards",
        "msPurity_L_tryptophan_6305_MSn.mzML.gz", package = "msPurity",
        mustWork = TRUE)
    raw <- Spectra::Spectra(path, source = Spectra::MsBackendMzR())
    pa <- purityA(raw, cores = 1L, plotP = FALSE)
    stored <- purityTable(pa, legacy = FALSE)
    expect_length(stored, 78L)
    expect_identical(tabulate(match(stored$activationType, c("CID", "HCD")),
            nbins = 2L), c(29L, 49L))
    expect_identical(stored$scanFilterAnalyzer, rep("FTMS", 78L))
    source_rows <- match(stored$spectrumId, raw$spectrumId)
    expect_identical(stored$collisionEnergy, raw$collisionEnergy[source_rows])

})
