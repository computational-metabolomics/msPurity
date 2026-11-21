context("purityA Spectra slots, accessors and updateObject")

purity_ref <- .parity_ref("purityA")$purityA$puritydf
msms <- unname(.lcmsms_paths())

test_that("purityA() stores the scans and their purity as Spectra", {
    pa <- suppressMessages(purityA(msms, offsets = c(0.5, 0.5)))
    sp <- purityTable(pa, legacy = FALSE)
    expect_s4_class(sp, "Spectra")
    expect_length(sp, nrow(purity_ref))
    expect_identical(sp$pid, purity_ref$pid)
    expect_equal(sp$inPurity, purity_ref$inPurity, tolerance = 1e-12)
    expect_true(all(Spectra::msLevel(sp) == 2L))
    expect_identical(Spectra::scanIndex(sp), purity_ref$seqNum)
    expect_equal(purityTable(pa), purity_ref, tolerance = 1e-12)
    expect_true(validObject(pa))
})

test_that("purityA() accepts Spectra and MsExperiment input", {
    sp <- Spectra::Spectra(msms, source = Spectra::MsBackendMzR(),
        BPPARAM = BiocParallel::SerialParam())
    pa <- suppressMessages(purityA(sp, offsets = c(0.5, 0.5)))
    expect_equal(purityTable(pa), purity_ref, tolerance = 1e-12)
    expect_identical(unname(basename(pa@fileList)), basename(msms))

    me <- MsExperiment::readMsExperiment(msms, BPPARAM = BiocParallel::SerialParam())
    pa <- suppressMessages(purityA(me, offsets = c(0.5, 0.5)))
    expect_equal(purityTable(pa), purity_ref, tolerance = 1e-12)
})

test_that("purityA() reads a Parquet study dataset through MsBackendParquet", {
    skip_if_no_arrow()
    sp <- Spectra::Spectra(Spectra::backendInitialize(
        MsBackendParquet::MsBackendParquet(), path = .study_dataset()))
    pa <- suppressMessages(purityA(sp, offsets = c(0.5, 0.5)))
    expect_equal(purityTable(pa), purity_ref, tolerance = 1e-12)
})

test_that("legacy slots can be switched off", {
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    pa <- suppressMessages(purityA(msms, offsets = c(0.5, 0.5)))
    expect_identical(nrow(pa@puritydf), 0L)
    expect_equal(purityTable(pa), purity_ref, tolerance = 1e-12)
})

test_that("updateObject() converts every saved purityA object", {
    files <- list.files(.historical_fixture("purityA"),
                        pattern = "_pa(_OLD)?\\.rds$", full.names = TRUE)
    expect_gt(length(files), 10)
    for (f in files) {
        old <- readRDS(f)
        expect_false(msPurity:::.pa_is_current(old))
        pa <- updateObject(old)
        expect_true(msPurity:::.pa_is_current(pa), label = basename(f))
        expect_true(validObject(pa))
        expect_identical(purityTable(pa), old@puritydf, label = basename(f))
        expect_equal(.parity_strip(groupedSpectra(pa)),
                     .parity_strip(.label_ms2(old@grped_ms2)),
                     label = basename(f))
        expect_equal(.parity_strip(averagedSpectra(pa)),
                     .parity_strip(old@av_spectra), label = basename(f))
        expect_identical(updateObject(pa), pa)
    }
})
