context("purityA Spectra slots, accessors and updateObject")

purity_ref <- .parity_ref("purityA")$purityA$puritydf
msms <- unname(.lcmsms_paths())

test_that("purityA() stores the scans and their purity as Spectra", {
    pa <- suppressMessages(purityA(msms))
    sp <- purityTable(pa, legacy = FALSE)
    expect_s4_class(sp, "Spectra")
    expect_length(sp, nrow(purity_ref))
    expect_identical(sp$pid, purity_ref$pid)
    expect_identical(sp$inPurity, purity_ref$inPurity)
    expect_true(all(Spectra::msLevel(sp) == 2L))
    expect_identical(Spectra::scanIndex(sp), purity_ref$seqNum)
    expect_identical(purityTable(pa), purity_ref)
    expect_true(validObject(pa))
})

test_that("purityA() accepts Spectra and MsExperiment input", {
    sp <- Spectra::Spectra(msms, source = Spectra::MsBackendMzR())
    pa <- suppressMessages(purityA(sp))
    expect_identical(purityTable(pa), purity_ref)
    expect_identical(unname(basename(pa@fileList)), basename(msms))

    me <- MsExperiment::readMsExperiment(msms)
    pa <- suppressMessages(purityA(me))
    expect_identical(purityTable(pa), purity_ref)
})

test_that("purityA() reads an mzStack study dataset through MsBackendParquet", {
    skip_if_no_mzstack()
    sp <- Spectra::Spectra(Spectra::backendInitialize(
        MsBackendParquet::MsBackendParquet(), path = .study_dataset()))
    pa <- suppressMessages(purityA(sp))
    expect_identical(purityTable(pa), purity_ref)
})

test_that("legacy slots can be switched off", {
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    pa <- suppressMessages(purityA(msms))
    expect_identical(nrow(pa@puritydf), 0L)
    expect_identical(purityTable(pa), purity_ref)
})

test_that("updateObject() converts every saved purityA object", {
    files <- list.files(system.file("extdata", "tests", "purityA",
                                    package = "msPurity"),
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

test_that("accessors return Spectra on request", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
                              "9_averageAllFragSpectra_with_filter_pa.rds",
                              package = "msPurity"))
    g <- groupedSpectra(pa, legacy = FALSE)
    expect_s4_class(g, "Spectra")
    expect_setequal(g$pid, unique(pa@grped_df$pid))
    expect_false(anyDuplicated(g$pid) > 0)
    expect_true(all(c("snr", "ra", "pass_flag") %in% Spectra::peaksVariables(g)))

    av <- averagedSpectra(pa, legacy = FALSE)
    expect_s4_class(av, "Spectra")
    expect_setequal(unique(av$av_level), c("av_intra", "av_inter", "av_all"))
    expect_true(all(c("frac", "count", "pass_flag") %in%
                        Spectra::peaksVariables(av)))
    a <- av[av$grpid == "187" & av$av_level == "av_all" & !av$is_null]
    expect_equal(Spectra::mz(a)[[1]],
                 pa@av_spectra[["187"]]$av_all$mz)

    expect_length(allFragSpectra(pa, legacy = FALSE), 0)
    expect_identical(nrow(allFragSpectra(pa)), 0L)
})

test_that("validity catches links to scans that do not exist", {
    pa <- updateObject(readRDS(system.file(
        "extdata", "tests", "purityA", "2_frag4feature_pa.rds",
        package = "msPurity")))
    pa@grped_df$pid[1] <- 1e6L
    expect_error(validObject(pa), "grped_df refers to pid values")
})
