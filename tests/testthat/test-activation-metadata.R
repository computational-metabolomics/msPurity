.activation_metadata_spectra <- function(labels, accessions) {
    metadata <- data.frame(pid = seq_along(labels), msLevel = 2L,
        activationType = labels, activationTypeAccession = accessions)
    peaks <- rep(list(cbind(mz = 146.0598, intensity = 100)), length(labels))
    msPurity:::.msp_memory(metadata, peaks)
}

test_that("activation selectors share the recognized label and CV mapping", {
    labels <- c("cid", "beam-type collision-induced dissociation",
        "electron transfer dissociation", NA_character_, "CID/HCD")
    accessions <- c("MS:1000133", "MS:1000422", "MS:1000598",
        "MS:1000422", NA_character_)
    source <- .activation_metadata_spectra(labels, accessions)
    param <- list(
        cid = FilterFragSpectraParam(select = list(
            activationType = "collision-induced dissociation")),
        hcd = FilterFragSpectraParam(select = list(
            activationType = "MS:1000422")),
        etd = FilterFragSpectraParam(select = list(activationType = "ETD")))
    fallback <- list(ilim = 0, plim = 0.8, ra = 0, snr = 3,
        snmeth = "median", rmp = FALSE, allfrag = FALSE, cores = 1)

    result <- msPurity:::.pa_filter_resolve(source, param, fallback)

    expect_identical(result$assignments$parameter_set_names,
        list("cid", "hcd", "etd", "hcd", character()))
    expect_identical(result$metadata$activationType,
        c(labels[1:3], "HCD", "CID/HCD"))
    expect_identical(source$activationType, labels)
    expect_identical(source$activationTypeAccession, accessions)
    expect_identical(as.list(result$parameter_objects$hcd)$select,
        list(activationType = "MS:1000422"))
})

test_that("structured activation precedes the limited scan-filter fallback", {
    source <- data.frame(msLevel = 2L,
        activationType = c("CID", NA, "vendor method", NA, NA),
        activationTypeAccession = c("MS:1000133", "MS:1000598", NA,
            "MS:9999999", NA),
        filterString = "FTMS Full ms2 205.0971@hcd30 [50-216]")

    result <- msPurity:::.pa_filter_metadata(source)

    expect_identical(result$activationType,
        c("CID", "ETD", "vendor method", NA_character_, "HCD"))
    expect_identical(result$scanFilterAnalyzer, rep("FTMS", 5L))
    expect_true(all(is.na(result$collisionEnergyUnit)))
})

test_that("equivalent activation aliases agree and known conflicts fail", {
    source <- data.frame(msLevel = 2L,
        activationType = "beam-type collision-induced dissociation",
        activationMethod = "HCD", activationTypeAccession = "MS:1000422")
    expect_identical(msPurity:::.pa_filter_metadata(source)$activationType,
        source$activationType)

    source$activationTypeAccession <- "MS:1000133"
    expect_error(msPurity:::.pa_filter_metadata(source),
        "Conflicting activation label and activationTypeAccession")
    source$activationTypeAccession <- "MS:1000422"
    source$activationMethod <- "CID"
    expect_error(msPurity:::.pa_filter_metadata(source),
        "Conflicting known acquisition aliases")
})

test_that("combined and unmapped activation identifiers are not collapsed", {
    source <- data.frame(msLevel = 2L,
        activationType = c("ETD+HCD", "custom activation", NA),
        activationTypeAccession = c("MS:1000598", "MS:9999999",
            "MS:1000598;MS:1000422"),
        filterString = "FTMS Full ms2 205.0971@cid30 [50-216]")
    expect_identical(msPurity:::.pa_filter_metadata(source)$activationType,
        c("ETD+HCD", "custom activation", NA_character_))
    expect_identical(msPurity:::.pa_activation_type(c("ETD+HCD", "other")),
        c("ETD+HCD", "other"))
})
