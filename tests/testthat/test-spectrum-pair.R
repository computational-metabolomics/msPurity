.pair_spectra <- function(level = 2L, precursor = 200,
    mz = c(50, 100), intensity = c(20, 100)) {
    d <- S4Vectors::DataFrame(msLevel = as.integer(level), precursorMz = precursor,
        acquisitionNum = seq_along(level), dataOrigin = "run1", polarity = 1L,
        activationType = "CID", massAnalyzer = "Orbitrap", scanFilterAnalyzer = "FTMS",
        collisionEnergy = 20, collisionEnergyUnit = "eV")
    d$mz <- rep(list(mz), length(level))
    d$intensity <- rep(list(intensity), length(level))
    Spectra::Spectra(d, source = Spectra::MsBackendMemory())
}

test_that("pair settings and indices are validated before reading peaks", {
    defaults <- as.list(SpectralMatchingParam())
    expect_equal(unlist(defaults), c(q_ppmPrec = 5, l_ppmPrec = 5, q_daPrec = 0,
        l_daPrec = 0, q_ppmProd = 10, l_ppmProd = 10, q_daProd = 0, l_daProd = 0,
        collisionEnergyTolerance = 0))
    for (field in names(defaults))
        for (bad in list(-1, Inf, NA_real_, NaN, numeric(), c(0, 1), "2"))
            expect_error(do.call(SpectralMatchingParam, stats::setNames(list(bad), field)))
    sp <- Spectra::addProcessing(.pair_spectra(), function(x, ...) stop("peaks read"))
    for (index in list(0L, 2L, 1.5, NA_integer_, integer(), c(1L, 1L), TRUE))
        expect_error(scoreSpectrumPair(sp, sp, queryIndex = index), "queryIndex")
    expect_error(scoreSpectrumPair(sp, sp, referenceIndex = 0L), "referenceIndex")
    p <- SpectralMatchingParam()
    p@l_daProd <- -1
    expect_error(scoreSpectrumPair(sp, sp, p), "l_daProd")
    expect_error(scoreSpectrumPair(sp, sp, param = list()), "SpectralMatchingParam")
    expect_error(scoreSpectrumPair(1, sp), "query must")
})

test_that("the public MS2 scorer preserves all six established scores", {
    q <- .pair_spectra(mz = 100.1234, intensity = 10)
    l <- .pair_spectra(mz = c(50.0617, 100.1234), intensity = c(10, 10))
    z <- scoreSpectrumPair(q, l)
    expect_identical(z$reason, "compatible")
    expect_true(z$eligible)
    expect_equal(z$scores, c(dpc = 4 / sqrt(17), rdpc = 4 / sqrt(17),
        cdpc = 8 / (3 * sqrt(17)), mcount = 1, allcount = 2, mpercent = 0.5),
        tolerance = 1e-14)
    expect_true(all(z$checks$outcome == "match"))
    l <- Spectra::filterMzRange(l, mz = c(49, 51))
    expect_equal(scoreSpectrumPair(q, l)$scores[["mcount"]], 0)
})

test_that("precursor and product intervals combine both sides and include endpoints", {
    # 100 * 1250 ppm = 0.125. Query half-width is 0.25, reference is 0.125.
    p <- SpectralMatchingParam(q_ppmPrec = 1250, l_ppmPrec = 0,
        q_daPrec = 0.125, l_daPrec = 0.125, q_ppmProd = 1250, l_ppmProd = 0,
        q_daProd = 0.125, l_daProd = 0.125)
    q <- .pair_spectra(precursor = 100, mz = 100, intensity = 1)
    pair <- function(prec, prod) scoreSpectrumPair(q,
        .pair_spectra(precursor = prec, mz = prod, intensity = 1), p)
    for (sign in c(-1, 1)) {
        edge <- 100 + sign * 0.375
        expect_equal(pair(edge, edge)$scores[["mcount"]], 1)
        expect_equal(pair(edge, edge + sign * 1e-11)$scores[["mcount"]], 1)
        expect_equal(pair(edge, edge + sign * 1e-9)$scores[["mcount"]], 0)
        expect_identical(pair(edge + sign * 1e-11, edge)$reason, "precursor_mismatch")
    }
    zero <- SpectralMatchingParam(q_ppmPrec = 0, l_ppmPrec = 0, q_ppmProd = 0,
        l_ppmProd = 0)
    expect_identical(scoreSpectrumPair(q, q, zero)$reason, "compatible")
    expect_false(scoreSpectrumPair(q, .pair_spectra(precursor = 100.00000001), zero)$eligible)
    l <- .pair_spectra(precursor = 101, mz = 101, intensity = 1)
    # At one percent the larger mass has a wider ppm window.
    q$precursorMz <- 99.995
    parms <- function(qppm, lppm) SpectralMatchingParam(q_ppmPrec = qppm,
        l_ppmPrec = lppm)
    expect_false(scoreSpectrumPair(q, l, parms(10000, 0))$eligible)
    expect_true(scoreSpectrumPair(q, l, parms(0, 10000))$eligible)
})

test_that("required fields and known conflicts reject pairs before reading peaks", {
    q <- Spectra::addProcessing(.pair_spectra(), function(x, ...) stop("peaks read"))
    compare <- function(field, value) {
        l <- .pair_spectra()
        l[[field]] <- value
        z <- scoreSpectrumPair(q, l)
        expect_false(z$eligible)
        expect_true(all(is.na(z$scores)))
        z$reason
    }
    expect_identical(compare("msLevel", 3L), "ms_level_mismatch")
    expect_identical(compare("msLevel", NA_integer_), "ms_level_missing")
    expect_identical(compare("msLevel", 1L), "ms_level_unsupported")
    expect_identical(compare("msLevel", 5L), "ms_level_unsupported")
    for (value in c(NA_real_, 0, -5))
        expect_identical(compare("precursorMz", value), "missing_precursor")
    for (field in c("activationType", "massAnalyzer", "scanFilterAnalyzer",
                    "collisionEnergyUnit", "polarity", "collisionEnergy")) {
        v <- switch(field, activationType = "HCD", massAnalyzer = "quadrupole",
            scanFilterAnalyzer = "ITMS", collisionEnergyUnit = "NCE",
            polarity = 0L, collisionEnergy = 30)
        expect_identical(compare(field, v), paste0(field, "_mismatch"))
    }
})

test_that("partial metadata and mixed energy scales retain uncertainty", {
    q <- .pair_spectra()
    l <- q
    l$activationType <- "MS:1000133"
    l$massAnalyzer <- "MS:1000484"
    l$collisionEnergyUnit <- "UO:0000266"
    l$collisionEnergy <- 21
    p <- SpectralMatchingParam(collisionEnergyTolerance = 1)
    expect_identical(scoreSpectrumPair(q, l, p)$reason, "compatible")
    l$collisionEnergyUnit <- NA_character_
    z <- scoreSpectrumPair(q, l)
    expect_identical(z$reason, "compatible_partial")
    expect_identical(z$checks$outcome[z$checks$field == "collisionEnergy"], "unit_unverified")
    l$msPurityAcquisitionValues <- list(list(activationType = c("ETD", "CID"),
        collisionEnergyUnit = c("eV", "NCE"), collisionEnergy = c(50, 100)))
    z <- scoreSpectrumPair(q, l)
    expect_identical(z$reason, "compatible_partial")
    expect_identical(z$checks$outcome[z$checks$field == "activationType"], "partial")
    l$msPurityAcquisitionValues <- list(list(activationType = c("ETD", "HCD")))
    expect_identical(scoreSpectrumPair(q, l)$reason, "activationType_mismatch")
    l <- q
    l$polarity <- 0L
    l$x_mspurity_acquisition_values <- '{"polarity":null}'
    l$x_mspurity_acquisition_status <- '{"polarity":"unknown"}'
    expect_identical(scoreSpectrumPair(q, l)$reason, "compatible_partial")
    l <- .pair_spectra(mz = numeric(), intensity = numeric())
    z <- scoreSpectrumPair(q, l)
    expect_true(z$eligible)
    expect_identical(z$reason, "no_usable_peaks")
    expect_true(all(is.na(z$scores)))
    expect_identical(scoreSpectrumPair(q, .pair_spectra(intensity = c(0, 0)))$reason,
        "no_usable_peaks")
})

test_that("MS3 and MS4 walk the supplied run's ancestry and disclose absent stages", {
    sp <- .pair_spectra(rep(2:4, 2), c(200, 150, 80, 210, 150, 80))
    sp$dataOrigin <- rep(c("a", "b"), each = 3L)
    sp$acquisitionNum <- rep(10:12, 2L)
    sp$precScanNum <- rep(c(NA_integer_, 10L, 11L), 2L)
    for (i in 2:3) {
        z <- scoreSpectrumPair(sp, sp, queryIndex = i, referenceIndex = i + 3L)
        expect_identical(z$reason, "ancestor_precursor_mismatch")
        expect_identical(z$checks$stage[z$checks$outcome == "mismatch"], 2L)
    }
    l <- sp[6:4]
    l$precursorMz <- c(80, 150, 200)
    expect_identical(scoreSpectrumPair(sp, l, queryIndex = 3)$reason, "compatible")
    l$activationType <- c("CID", "HCD", "CID")
    z <- scoreSpectrumPair(sp, l, queryIndex = 3)
    expect_identical(z$reason, "ancestor_activationType_mismatch")
    z <- scoreSpectrumPair(sp, l[1], queryIndex = 3)
    expect_identical(z$reason, "compatible_partial")
    expect_true(all(z$checks$outcome[z$checks$field == "precursorMz" &
        z$checks$stage < 4L] == "missing"))
    expect_identical(z$checks$outcome[z$checks$field == "ancestorStage"],
        c("missing", "missing"))
    # Neither run needs to use the other's ID convention.
    l$spectrumId <- c("leaf", "middle", "root")
    l$spectrumRef <- c("middle", "root", NA_character_)
    l$precScanNum <- rep(NA_integer_, 3L)
    expect_identical(scoreSpectrumPair(sp, l, queryIndex = 3)$reason,
        "ancestor_activationType_mismatch")
})

test_that("measured MS4 in the bundled tryptophan run resolves both ancestors", {
    f <- .msn_standard_file("tryptophan", "MSn")
    sp <- Spectra::Spectra(f, source = Spectra::MsBackendMzR(),
        BPPARAM = BiocParallel::SerialParam())
    rows <- which(Spectra::msLevel(sp) == 4L)
    i <- rows[which(lengths(Spectra::mz(sp[rows])) > 0L)[1L]]
    z <- scoreSpectrumPair(sp, sp, queryIndex = i, referenceIndex = i)
    mass <- z$checks[z$checks$field == "precursorMz", ]
    expect_identical(mass$stage, 4:2)
    expect_identical(mass$outcome, rep("match", 3))
    expect_true(z$eligible)
    expect_gt(z$scores[["mcount"]], 0)
})

test_that("consensus path masses and acquisition conditions are compared by stage", {
    q <- .pair_spectra(4L, precursor = 80)
    q$msPurityPrecursorPath <- list(data.frame(msLevel = 3:4, precursorMz = c(150, 80)))
    q$msPurityAncestorConditions <- list(data.frame(msLevel = 2:3,
        activationType = c("CID", "HCD")))
    l <- .pair_spectra(4L, precursor = 80)
    l$x_mspurity_precursor_path <- '{"msLevel":[4,3],"precursorMz":[80,150]}'
    l$x_mspurity_ancestor_conditions <- '{"msLevel":[3,2],"activationType":["HCD","CID"]}'
    z <- scoreSpectrumPair(q, l)
    expect_identical(z$reason, "compatible_partial")
    masses <- z$checks[z$checks$field == "precursorMz", ]
    expect_equal(masses$outcome[match(2:4, masses$stage)], c("missing", "match", "match"))
    l$x_mspurity_precursor_path <- '{"msLevel":[4,3],"precursorMz":[80,145]}'
    expect_identical(scoreSpectrumPair(q, l)$reason, "ancestor_precursor_mismatch")
    l$x_mspurity_ancestor_conditions <- '{"msLevel":2,"activationType":"ETD"}'
    l$x_mspurity_precursor_path <- NA_character_
    z <- scoreSpectrumPair(q, l)
    expect_identical(z$reason, "ancestor_activationType_mismatch")
    expect_identical(z$checks$stage[z$checks$outcome == "mismatch"], 2L)
})

test_that("ambiguous scan numbers do not supply an ancestor", {
    sp <- .pair_spectra(c(2L, 2L, 3L), c(200, 210, 150))
    sp$acquisitionNum <- c(10L, 11L, 12L)
    sp$scanIndex <- c(9L, 10L, 11L)
    sp$precScanNum <- c(NA_integer_, NA_integer_, 10L)
    z <- scoreSpectrumPair(sp, sp, queryIndex = 3, referenceIndex = 3)
    expect_identical(z$reason, "compatible_partial")
    expect_identical(z$checks$outcome[z$checks$field == "precursorMz" &
        z$checks$stage == 2L], "missing")
})

test_that("a controlled analyzer identifier resolves different source labels", {
    q <- .pair_spectra()
    l <- q
    q$massAnalyzer <- "time-of-flight"
    l$massAnalyzer <- "TOF"
    q$massAnalyzerAccession <- "MS:1000084"
    l$massAnalyzerAccession <- "MS:1000084"
    expect_identical(scoreSpectrumPair(q, l)$reason, "compatible")
    l$massAnalyzerAccession <- "MS:1000484"
    expect_identical(scoreSpectrumPair(q, l)$reason, "massAnalyzer_mismatch")
    l$massAnalyzerAccession <- NA_character_
    q$massAnalyzer <- NA_character_
    expect_identical(scoreSpectrumPair(q, l)$reason, "compatible_partial")
})
