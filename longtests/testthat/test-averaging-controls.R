test_that("intra, inter and all averaging honour call-level controls", {
    pa <- updateObject(.golden_pa("2_frag4feature_pa.rds"))
    groups <- head(unique(as.character(pa@grped_df$grpid)), 2L)
    pa@grped_df <- pa@grped_df[pa@grped_df$grpid %in% groups, ]
    pa@fragSpectra <- pa@fragSpectra[pa@fragSpectra$pid %in% pa@grped_df$pid]
    pa@grped_ms2 <- msPurity:::.pa_grouped_legacy(pa)
    pa@cores <- 2
    serial <- parallel <- pa
    for (level in c("Intra", "Inter", "All")) {
        step <- paste0("average", level, "FragSpectra")
        serial <- do.call(step, list(serial, cores = 1, blockSize = 100L))
        parallel <- do.call(step, list(parallel, cores = 2))
        expect_equal(averagedSpectra(serial), averagedSpectra(parallel), info = step)
        sp <- averagedSpectra(serial, legacy = FALSE)
        pp <- averagedSpectra(parallel, legacy = FALSE)
        expect_identical(sp$msPurityClusterContributors,
            pp$msPurityClusterContributors, info = step)
        parameters <- tail(processingParameters(serial, step), 1L)[[1L]]
        expect_identical(parameters[c("cores", "effective_cores", "blockSize")],
            list(cores = 1, effective_cores = 1L, blockSize = 100L), info = step)
        expect_identical(tail(processingParameters(parallel, step), 1L)[[1L]]$effective_cores, 2L, info = step)
        expect_identical(serial@cores, 2)
    }
    inherited <- averageAllFragSpectra(pa)
    expect_identical(tail(processingParameters(inherited,
        "averageAllFragSpectra"), 1L)[[1L]]$effective_cores, 2L)
    pa@cores <- 0
    expect_identical(tail(processingParameters(averageAllFragSpectra(pa),
        "averageAllFragSpectra"), 1L)[[1L]]$effective_cores, 1L)
})
