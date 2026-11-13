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
