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

test_that("DIMS averaging and grouping forward execution controls", {
    path <- tempfile(fileext = ".csv")
    on.exit(unlink(path), add = TRUE)
    write.csv(data.frame(mz = rep(c(100, 200), 2), i = c(10, 20, 30, 40),
        scanid = rep(1:2, each = 2)), path, row.names = FALSE)
    files <- data.frame(filepth = rep(path, 2), name = c("a", "b"),
        sampleType = "sample", class = "sample")
    pd <- purityD(files, mzML = FALSE, cores = 2)
    single <- averageSpectraSingle(path, csvFile = TRUE, snthr = 0,
        cores = 1, blockSize = 2L)
    expect_equal(single$mz, c(100, 200))
    serial <- averageSpectra(pd, snthr = 0, cores = 1, blockSize = 2L)
    parallel <- averageSpectra(pd, snthr = 0)
    expect_equal(averagedPeaks(serial), averagedPeaks(parallel))
    expect_identical(serial@avParam[c("cores", "blockSize")],
        list(cores = 1, blockSize = 2L))
    expect_identical(parallel@avParam$cores, 2)
    expect_identical(serial@cores, 2)
    serial <- groupPeaks(serial, cores = 1, blockSize = 2L)
    parallel <- groupPeaks(parallel, cores = 2, blockSize = 2L)
    expect_equal(serial@groupedPeaks, parallel@groupedPeaks)
    expect_identical(serial@cores, 2)
})
