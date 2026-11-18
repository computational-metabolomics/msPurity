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

test_that("an activation recorded as a PSI-MS accession reads as its method", {
    meta <- function(labels, accessions = rep(NA_character_, length(labels)),
                     filters = NULL) {
        d <- data.frame(msLevel = 2L, activationType = labels,
                        activationTypeAccession = accessions)
        if (!is.null(filters)) d$filterString <- filters
        msPurity:::.pa_filter_metadata(d)
    }
    ## A Parquet study records the accession; mzML input names the method
    ## through its filter string. Both read as the same method.
    study <- meta(c("MS:1000422", "MS:1000133", NA))
    expect_identical(study$activationType, c("HCD", "CID", NA))
    expect_identical(attr(study, "activationTypeAccession"),
                     c("MS:1000422", "MS:1000133", NA))
    expect_identical(attr(study, "basis")$activationType,
                     c("recorded", "recorded", NA))
    mzml <- meta(rep(NA_character_, 2L), filters = c(
        "FTMS + p ESI d Full ms2 205.0971@hcd30.00 [50.0000-215.0000]",
        "ITMS + c ESI d Full ms2 205.0971@cid35.00 [50.0000-215.0000]"))
    expect_identical(mzml$activationType, c("HCD", "CID"))
    expect_identical(attr(mzml, "basis")$activationType,
                     rep("filter_string", 2L))
    ## A recorded accession is kept even where a filter string disagrees.
    both <- meta("MS:1000422", filters =
        "ITMS + c ESI d Full ms2 205.0971@cid35.00 [50.0000-215.0000]")
    expect_identical(both$activationType, "HCD")
    expect_error(meta("MS:1000422", "MS:1000133"),
                 "Conflicting activation label")
    ## An accession msPurity does not map stays as recorded.
    ecd <- meta("MS:1000250")
    expect_identical(ecd$activationType, "MS:1000250")
    expect_identical(attr(ecd, "basis")$activationType, "recorded")
})
