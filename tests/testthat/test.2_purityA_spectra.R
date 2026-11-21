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

test_that("a historical purityA stage upgrades without changing its purity table", {
    legacy <- .golden_pa("9_averageAllFragSpectra_with_filter_pa_OLD.rds")
    expect_false(msPurity:::.pa_is_current(legacy))
    current <- updateObject(legacy)
    expect_true(validObject(current))
    expect_identical(purityTable(current), legacy@puritydf)
    expect_identical(updateObject(current), current)
})
