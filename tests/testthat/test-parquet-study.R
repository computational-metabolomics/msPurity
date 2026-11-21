context("Raw MS1-MS4 studies converted by mzMLToParquet()")

test_that("a selected or processed study input numbers its own peaks", {
    skip_if_no_study()
    files <- .d02_copy(c("kynurenine", "tryptophan"))
    study <- .msn_study(files)
    sp <- .msn_study_spectra(study)
    trp <- Spectra::filterDataOrigin(sp, normalizePath(files[2L]))
    ref <- purityA(files[2L], cores = 1L, plotP = FALSE, msLevel = 2:4)
    pa <- purityA(trp, cores = 1L, plotP = FALSE, msLevel = 2:4)
    for (v in c("msPurityRawPeakOffset", "msPurityRawPeakCount", "inPurity",
                "msPurityStagePurity"))
        expect_identical(pa@spectra[[v]], ref@spectra[[v]], info = v)
    queued <- function(x) Spectra::filterIntensity(x, intensity = c(1e4, Inf))
    mz <- Spectra::Spectra(files[2L], source = Spectra::MsBackendMzR())
    ref <- purityA(queued(mz), cores = 1L, plotP = FALSE, msLevel = 2:4)
    pa <- purityA(queued(trp), cores = 1L, plotP = FALSE, msLevel = 2:4)
    for (v in c("msPurityRawPeakOffset", "msPurityRawPeakCount", "inPurity"))
        expect_identical(pa@spectra[[v]], ref@spectra[[v]], info = v)
    expect_equal(as.list(Spectra::peaksData(
                     msPurity:::.pa_raw_scans(pa, pa@spectra$pid[1:5]))),
                 as.list(Spectra::peaksData(
                     msPurity:::.pa_raw_scans(ref, ref@spectra$pid[1:5]))))
})

test_that("a moved or changed study is reported, not read from mzML", {
    skip_if_no_study()
    files <- .d02_copy()
    study <- .copy_dataset(.msn_study(.msn_standard_file("tryptophan",
                                                         "MSn")))
    pa <- purityA(.msn_study_spectra(study), cores = 1L, plotP = FALSE)
    pids <- pa@spectra$pid[1:3]
    expect_s4_class(msPurity:::.pa_raw_scans(pa, pids), "Spectra")
    moved <- paste0(study, "-moved")
    file.rename(study, moved)
    MsBackendParquet::invalidateDatasetCache(study)
    expect_parquet_error(msPurity:::.pa_raw_scans(pa, pids), "reference")
    expect_parquet_error(getScanPeaks(pa), "reference")
    file.rename(moved, study)
    MsBackendParquet::invalidateDatasetCache(study)
    m <- .manifest(study)
    m$runs[[1L]]$ingested_at <- m$runs[[1L]]$ingested_at + 1L
    .write_raw_manifest(study, m)
    expect_parquet_error(msPurity:::.pa_raw_scans(pa, pids), "stale")
    expect_parquet_error(getScanPeaks(pa), "stale")
})

test_that("datasets written before MSn export stay valid and readable", {
    skip_if_no_arrow()
    dir <- tempfile("v1.0-")
    utils::unzip(testthat::test_path("fixtures", "parquet-v1.0.zip"),
                 exdir = dir)
    p <- file.path(dir, "d02_ms2.parquet")
    expect_identical(.manifest(p)$results$version, "1.0.0")
    expect_false(any(validateParquet(p)$level == "MUST"))
    sc <- .read_table(p, "x_mspurity_scan")
    info <- processedSpectrumInfo(p)
    sp <- readProcessedSpectra(p)
    expect_s4_class(sp@backend, "MsBackendParquet")
    expect_equal(info$peak_count, unname(lengths(Spectra::mz(sp))))
    expect_true(all(is.na(info$used_measurements)))
    empty <- processedSpectrumInfo(p, featureId = -1L)
    expect_identical(names(empty), names(info))
    expect_identical(nrow(empty), 0L)
    expect_length(readProcessedSpectra(p, featureId = -1L), 0L)
    expect_equal(nrow(readProcessedProvenance(p, featureId = -1L)$mergeMembers), 0L)
    local_mocked_bindings(peaksData = function(...) stop("peak arrays were requested"),
                          .package = "Spectra")
    expect_equal(processedSpectrumInfo(p)$spectrum_id_, info$spectrum_id_)
    expect_equal(readProcessedProvenance(p)$outputs$spectrum_id_, info$spectrum_id_)
    expect_null(sc$ms_level)
    expect_identical(msPurity:::.mzs_ms2_scans(sc), sc)
})
