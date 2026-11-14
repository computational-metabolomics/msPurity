context("D05 acquisition metadata and original peak positions")

.d05_metadata_fixture <- function(duplicate = FALSE,
    unit_accessions = rep("UO:0000266", 2)) {
    d <- S4Vectors::DataFrame(msLevel = c(2L, 2L), scanIndex = 1:2,
        acquisitionNum = 101:102, rtime = c(10, 11), precursorMz = 205,
        dataOrigin = "memory", polarity = 1L, massAnalyzer = "Orbitrap",
        massAnalyzerAccession = "MS:1000484",
        scanFilterAnalyzer = "FTMS",
        activationType = c("HCD", "CID"),
        activationTypeAccession = c("MS:1000422", "MS:1000133"),
        collisionEnergy = c(20, 40),
        collisionEnergyUnit = "eV",
        collisionEnergyUnitAccession = unit_accessions,
        spectrumId = c("scan=101", "scan=102"))
    d$mz <- list(c(100, 150), c(100, 150))
    d$intensity <- list(c(100, 50), c(2, 1))
    sp <- Spectra::Spectra(d, source = Spectra::MsBackendMemory())
    tab <- data.frame(pid = 1:2, fileid = 1L, seqNum = 1:2,
        acquisitionNum = 101:102, retentionTime = c(10, 11),
        precursorMZ = 205, precursorScanNum = NA_integer_, inPurity = 0.9)
    pa <- methods::new("purityA", fileList = c(memory = "memory"), cores = 1)
    pa@spectra <- msPurity:::.pa_spectra_new(sp, c(1L, 1L), tab)
    pa@grped_df <- data.frame(grpid = 1L, pid = 1:2, sample = 1L,
        cid = 1:2, inPurity = 0.9, precurMtchMZ = 205,
        retentionTime = c(10, 11))
    if (duplicate) pa@grped_df <- rbind(pa@grped_df, pa@grped_df[1, ])
    matrices <- as.list(Spectra::peaksData(sp))
    linked <- matrices[pa@grped_df$pid]
    pa@fragSpectra <- msPurity:::.pa_frag_from_grouped(pa, list(`1` = linked))
    msPurity:::.pa_sync_legacy(pa, puritydf = tab,
        grped_ms2 = list(`1` = linked))
}

test_that("source acquisition dimensions survive grouping and filtering", {
    pa <- .d05_metadata_fixture()
    out <- filterFragSpectra(pa, plim = 0, snr = 0)
    sp <- groupedSpectra(out, legacy = FALSE)
    expect_identical(sp$activationType, c("HCD", "CID"))
    expect_equal(sp$collisionEnergy, c(20, 40))
    expect_identical(sp$massAnalyzer, rep("Orbitrap", 2))
    expect_identical(sp$collisionEnergyUnit, rep("eV", 2))
    expect_identical(sp$activationTypeAccession,
        c("MS:1000422", "MS:1000133"))
    expect_identical(sp$massAnalyzerAccession, rep("MS:1000484", 2))
    expect_identical(sp$spectrumId, c("scan=101", "scan=102"))
    expect_equal(Spectra::polarity(sp), c(1, 1))
})

test_that("averages summarize used measurements and preserve repeated links", {
    pa <- .d05_metadata_fixture(duplicate = TRUE,
        unit_accessions = c("UO:0000266", NA_character_))
    pa <- filterFragSpectra(pa, plim = 0, snr = 0)
    mixed <- averageAllFragSpectra(pa, minfrac = 0)
    sp <- averagedSpectra(mixed, legacy = FALSE)
    expect_identical(sp$activationType, "mixed")
    expect_true(is.na(sp$collisionEnergy))
    expect_equal(sp$msPurityAcquisitionValues[[1]]$collisionEnergy, c(20, 40))
    expect_identical(sp$msPurityAcquisitionStatus[[1]]$activationType, "mixed")
    expect_identical(sp$msPurityAcquisitionConfiguration, "mixed")
    expect_true(is.na(sp$activationTypeAccession))
    expect_identical(sp$massAnalyzerAccession, "MS:1000484")
    expect_identical(sp$msPurityAcquisitionValues[[1]]$activationTypeAccession,
        c("MS:1000422", "MS:1000133"))
    expect_identical(sp$msPurityAcquisitionStatus[[1]]$activationTypeAccession,
        "mixed")
    expect_true(is.na(sp$collisionEnergyUnitAccession))
    expect_identical(
        sp$msPurityAcquisitionValues[[1]]$collisionEnergyUnitAccession,
        c("UO:0000266", NA_character_))
    expect_identical(
        sp$msPurityAcquisitionStatus[[1]]$collisionEnergyUnitAccession,
        "mixed")
    expect_equal(nrow(sp$msPurityContributors[[1]]), 3)
    expect_equal(sp$msPurityContributors[[1]]$pid, c(1, 2, 1))
    expect_equal(Spectra::polarity(sp), 1)
    out <- averageAllFragSpectra(filterFragSpectra(
        .d05_metadata_fixture(unit_accessions = c("UO:0000266", NA_character_)),
        ilim = 10, plim = 0, snr = 0), minfrac = 0)
    used <- averagedSpectra(out, legacy = FALSE)
    expect_identical(used$activationType, "HCD")
    expect_identical(used$activationTypeAccession, "MS:1000422")
    expect_identical(used$collisionEnergyUnitAccession, "UO:0000266")
    expect_equal(used$collisionEnergy, 20)
    expect_equal(used$msPurityContributors[[1]]$pid, 1)
})

test_that("unknown acquisition information stays unknown", {
    pa <- .d05_metadata_fixture()
    pa@fragSpectra$massAnalyzer <- rep(NA_character_, 2)
    pa@fragSpectra$activationType <- rep(NA_character_, 2)
    pa@fragSpectra$massAnalyzerAccession <- rep(NA_character_, 2)
    pa@fragSpectra$activationTypeAccession <- rep(NA_character_, 2)
    sp <- averagedSpectra(averageAllFragSpectra(pa, minfrac = 0),
        legacy = FALSE)
    expect_true(is.na(sp$massAnalyzer))
    expect_identical(sp$msPurityAcquisitionStatus[[1]]$massAnalyzer, "unknown")
    expect_true(is.na(sp$activationTypeAccession))
    expect_true(is.na(sp$massAnalyzerAccession))
    expect_identical(sp$msPurityAcquisitionStatus[[1]]$activationTypeAccession,
        "unknown")
    expect_identical(sp$msPurityAcquisitionStatus[[1]]$massAnalyzerAccession,
        "unknown")
})

test_that("average summaries resolve supported acquisition aliases", {
    pa <- .d05_metadata_fixture()
    pa@fragSpectra$massAnalyzer <- NULL
    pa@fragSpectra$activationType <- NULL
    pa@fragSpectra$analyzer <- c("Orbitrap", "Orbitrap")
    pa@fragSpectra$activationMethod <- c("CID", "CID")
    pa@fragSpectra$activationTypeAccession <- rep("MS:1000133", 2)
    pa@fragSpectra$collisionEnergy <- c(20, 20)

    sp <- averagedSpectra(averageAllFragSpectra(pa, minfrac = 0),
        legacy = FALSE)

    expect_identical(sp$massAnalyzer, "Orbitrap")
    expect_identical(sp$activationType, "CID")
    expect_identical(sp$msPurityAcquisitionConfiguration, "uniform")
    expect_identical(sp$analyzer, "Orbitrap")
    expect_identical(sp$activationMethod, "CID")
})

test_that("used occurrences follow retained all-member average clusters", {
    pa <- .d05_metadata_fixture(duplicate = TRUE)
    pa@fragSpectra$mz <- list(c(100, 150), c(200, 250))

    out <- averageAllFragSpectra(pa, minnum = 2, minfrac = 0, rmp = TRUE)
    sp <- averagedSpectra(out, legacy = FALSE)

    expect_equal(as.numeric(Spectra::mz(sp)[[1]]), c(100, 150))
    expect_equal(sp$msPurityContributors[[1]]$featureLinkIndex, c(1, 3))
    expect_equal(sp$msPurityContributors[[1]]$pid, c(1, 1))
    expect_identical(sp$activationType, "HCD")
    expect_equal(sp$collisionEnergy, 20)
})

test_that("inter-file summaries follow passed intra-file clusters", {
    pa <- .d05_metadata_fixture(duplicate = TRUE)
    pa@fragSpectra$mz <- list(c(100, 150), c(200, 250))
    intra <- averageIntraFragSpectra(pa, minnum = 2, minfrac = 0,
        rmp = FALSE)
    before <- averagedSpectra(intra, legacy = FALSE)
    expect_equal(before$msPurityContributors[[1]]$pid, c(1, 2, 1))

    out <- averageInterFragSpectra(intra, minfrac = 0)
    sp <- averagedSpectra(out, legacy = FALSE)
    sp <- sp[sp$av_level == "av_inter"]

    expect_equal(sp$msPurityContributors[[1]]$featureLinkIndex, c(1, 3))
    expect_equal(sp$msPurityContributors[[1]]$pid, c(1, 1))
    expect_identical(sp$activationType, "HCD")
    expect_equal(sp$collisionEnergy, 20)
    expect_equal(as.numeric(Spectra::mz(sp)[[1]]), c(100, 150))
})

test_that("queued subsets retain original peak positions", {
    pa <- .d05_metadata_fixture()
    pa@fragSpectra <- Spectra::filterMzRange(pa@fragSpectra, mz = c(125, 175))
    sp <- groupedSpectra(filterFragSpectra(pa,
            param = FilterFragSpectraParam(plim = 0, snr = 0)), legacy = FALSE)
    peaks <- Spectra::peaksData(sp, columns = c("mz", "original_peak_index"))
    expect_equal(peaks[[1]][, "mz"], 150)
    expect_equal(peaks[[1]][, "original_peak_index"], 2)
    expect_equal(peaks[[2]][, "original_peak_index"], 2)
    expect_equal(lengths(Spectra::mz(pa@fragSpectra)), c(1L, 1L))
    empty <- Spectra::filterMzRange(pa@fragSpectra, mz = c(500, 600))
    input <- msPurity:::.pa_filter_source_peaks(empty)
    expect_equal(lengths(input$original_positions), c(0L, 0L))
})

test_that("untracked ambiguous or value-changing queues fail explicitly", {
    d <- S4Vectors::DataFrame(msLevel = 2L)
    d$mz <- list(c(100, 100))
    d$intensity <- list(c(20, 20))
    sp <- Spectra::Spectra(d, source = Spectra::MsBackendMemory())
    one <- Spectra::addProcessing(sp, function(x) x[1, , drop = FALSE])
    expect_error(msPurity:::.pa_filter_source_peaks(one), "ambiguous")
    changed <- Spectra::addProcessing(sp, function(x) {
        x[, "intensity"] <- x[, "intensity"] * 2
        x
    })
    expect_error(msPurity:::.pa_filter_source_peaks(changed),
        "cannot be recovered")
    expect_equal(msPurity:::.pa_filter_source_peaks(sp)$original_positions[[1]],
        1:2)
})

test_that("legacy unavailable-average rows retain their declared MS2 level", {
    pa <- filterFragSpectra(.d05_metadata_fixture(), ilim = 1000,
        plim = 0, snr = 0)
    old_shape <- list(`1` = list(av_intra = list(`1` = NULL),
            av_inter = NULL, av_all = NULL))
    previous <- pa@avSpectra
    pa@avSpectra <- msPurity:::.pa_av_to_spectra(old_shape, pa@grped_df)
    pa <- msPurity:::.pa_av_metadata(pa, "intra", previous)
    sp <- averagedSpectra(pa, legacy = FALSE)
    expect_true(all(sp$is_null))
    expect_equal(Spectra::msLevel(sp), rep(2L, length(sp)))
    expect_identical(sp$msPurityAcquisitionStatus[[1]]$msLevel, "unknown")
})
