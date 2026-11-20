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
        direct <- spectralMatching(q, q)$matchedResults
        expect_identical(direct$msLevel, 4:2)
        expect_true(all(direct$mcount > 0L))
        expect_false(anyNA(direct[scores]))
        for (hit in list(spectralMatching(q, r)$matchedResults,
                         spectralMatching(r, r)$matchedResults)) {
            expect_identical(hit$msLevel, direct$msLevel)
            expect_equal(hit[scores], direct[scores], tolerance = 1e-12)
        }
    }
    compare(q, r)

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
    mixed <- spectralMatching(measured, r)$matchedResults
    expect_setequal(unique(mixed$msLevel), 2:4)
    expect_true(all(mixed$query_native_id == measured$spectrumId[mixed$queryIndex]))
    expect_identical(mixed$reference_native_id, r$x_mspurity_spectrum_key[mixed$referenceIndex])
    # A results path supplies experimental references with their original keys.
    found <- spectralMatching(q, path)
    expect_identical(found$referenceRole, "experimental")
    hit <- found$matchedResults
    self <- hit[hit$query_native_id == hit$reference_native_id, ]
    expect_identical(self$query_native_id, q$spectrumId)
    expect_equal(self[scores], spectralMatching(q, r)$matchedResults[scores],
        ignore_attr = TRUE)
})
