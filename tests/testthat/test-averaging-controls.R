test_that("DIMS averaging and grouping forward execution controls", .with_dims_deprecation({
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
}))
