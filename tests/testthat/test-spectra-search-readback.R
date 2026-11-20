test_that("offline measured and consensus MS2-MS4 searches survive Parquet readback", {
    skip_if_no_arrow()
    skip_if_not_installed("xcms")
    f <- .msn_standard_file("tryptophan", "MSn")
    study <- .msn_study(f)
    measured <- Spectra::Spectra(f, source = Spectra::MsBackendMzR(),
        BPPARAM = BiocParallel::SerialParam())
    reopened <- .msn_study_spectra(study)

    # Select an MS4 scan with two ancestors and signal at all three stages.
    num <- measured$acquisitionNum
    parent <- match(measured$precScanNum, num)
    has_signal <- lengths(Spectra::mz(measured)) > 0L
    leaf <- which(measured$msLevel == 4L & has_signal & has_signal[parent] &
        has_signal[parent[parent]])[1L]
    chain <- c(leaf, parent[leaf], parent[parent[leaf]])
    expect_identical(measured$msLevel[chain], 4:2)
    q <- measured[chain]
    r <- reopened[match(num[chain], reopened$acquisitionNum)]
    scores <- c("dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent")
    compare <- function(q, r) {
        result <- spectralMatching(q, q)
        direct <- result$matchedResults
        expect_true(all(direct$sharedMeasurement & direct$sharedSourceFile))
        expect_identical(direct$msLevel, 4:2)
        expect_true(all(direct$mcount > 0L))
        expect_false(anyNA(direct[scores]))
        for (reopened in list(spectralMatching(q, r), spectralMatching(r, r))) {
            hit <- reopened$matchedResults
            expect_identical(reopened$queryStatus$status, result$queryStatus$status)
            expect_identical(reopened$queryStatus$nEligible, result$queryStatus$nEligible)
            expect_identical(hit$referenceType, direct$referenceType)
            expect_true(all(hit$sharedMeasurement & hit$sharedSourceFile))
            expect_identical(hit$msLevel, direct$msLevel)
            expect_equal(hit[scores], direct[scores], tolerance = 1e-12)
        }
    }
    compare(q, r)
    # The native identifiers and acquisition numbers are repeated in a
    # different source file. Reverse the reopened reference selection.
    other <- Spectra::setBackend(q, Spectra::MsBackendMemory())
    other$dataOrigin <- rep(file.path(dirname(f), "replicate", basename(f)), length(q))
    separate <- spectralMatching(other, r[3:1])$matchedResults
    expect_identical(separate$referenceIndex, 3:1)
    expect_identical(separate$sharedMeasurement, rep(FALSE, 3))
    expect_identical(separate$sharedSourceFile, rep(FALSE, 3))

    xe <- .msn_standard_feature("tryptophan", "MSn")$xe
    pa <- purityA(f, cores = 1L, plotP = FALSE, msLevel = 2:4)
    pa <- frag4feature(pa, xe, ppm = 10, convert2RawRT = FALSE, msLevel = 2:4)
    pa <- filterFragSpectra(pa, allfrag = TRUE, msLevel = 2:4,
        param = FilterFragSpectraParam(plim = NA_real_, snr = 0))
    pa <- averageAllFragSpectra(pa, param = AverageFragSpectraParam(
        minfrac = 0, minnum = 0, snr = 0))
    live <- averagedSpectra(pa, legacy = FALSE)
    available <- which(lengths(Spectra::mz(live)) > 0L & is.finite(live$precursorMz))
    q <- live[available[match(4:2, live$msLevel[available])]]
    out <- tempfile()
    dir.create(out)
    path <- suppressMessages(createDatabase(pa, xe, outDir = out,
        dbName = "results.parquet", format = "parquet", study = study))
    stored <- readProcessedSpectra(path, msLevel = 2:4, method = "all")
    r <- stored[match(q$spectrumId, stored$x_mspurity_spectrum_key)]
    compare(q, r)
    comparison <- spectralMatching(measured, r)
    mixed <- comparison$matchedResults
    expect_true(all(mixed$queryType == "measured" & mixed$referenceType == "consensus"))
    expect_true(all(mixed$sharedSourceFile))
    for (j in seq_along(r)) {
        inputs <- averagePeakMembers(q, q$spectrumId[j], depth = "original")
        h <- mixed[mixed$referenceIndex == j, ]
        expect_identical(h$sharedMeasurement, h$query_native_id %in% inputs$member)
        sources <- subset(comparison$sourceMeasurements, side == "reference" & inputIndex == j)
        expect_setequal(sources$native_id, inputs$member)
        expect_equal(nrow(sources), length(unique(inputs$member)))
    }
    expect_true(any(mixed$sharedMeasurement))
    expect_setequal(unique(mixed$msLevel), 2:4)
    expect_true(all(mixed$query_native_id == measured$spectrumId[mixed$queryIndex]))
    expect_identical(mixed$reference_native_id, r$x_mspurity_spectrum_key[mixed$referenceIndex])
    # A results path supplies experimental references with their original keys.
    found <- spectralMatching(q, path)
    expect_identical(found$referenceRole, "experimental")
    hit <- found$matchedResults
    self <- hit[hit$query_native_id == hit$reference_native_id, ]
    expect_true(all(self$sharedMeasurement & self$sharedSourceFile))
    expect_identical(self$query_native_id, q$spectrumId)
    expect_equal(self[scores], spectralMatching(q, r)$matchedResults[scores],
        ignore_attr = TRUE)
})

test_that("a Parquet library accession does not identify an experimental measurement", {
    skip_if_no_arrow()
    d <- S4Vectors::DataFrame(msLevel = 2L, precursorMz = 200, polarity = 1L,
        spectrumId = "scan=1", accession = "entry-1", dataOrigin = "/experiment/raw.mzML")
    d$mz <- list(c(50, 80))
    d$intensity <- list(c(100, 20))
    q <- Spectra::Spectra(d, source = Spectra::MsBackendMemory())
    path <- tempfile()
    convertLibraryToParquet(q, path, version = "fixture")
    res <- spectralMatching(q, path)
    expect_identical(res$queryStatus$status, "scored")
    expect_identical(res$matchedResults$referenceType, "library")
    expect_identical(res$matchedResults$library_accession, "entry-1")
    expect_identical(res$matchedResults$sharedMeasurement, NA)
    expect_identical(res$matchedResults$sharedSourceFile, NA)
})
