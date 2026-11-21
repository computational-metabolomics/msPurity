test_that("frag4feature() and createDatabase() accept XcmsExperiment", {
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    skip_if_not_installed("MsExperiment")
    skip_if_not_installed("MSnbase")
    suppressPackageStartupMessages({
        library(xcms)
        library(MsExperiment)
    })
    files <- unname(.lcmsms_paths())
    cwp <- xcms::CentWaveParam(snthresh = 10, noise = 5e5, ppm = 10,
                               peakwidth = c(3, 30))
    pdp <- xcms::PeakDensityParam(sampleGroups = c(1, 1), minFraction = 0,
                                  bw = 30)
    xn <- MSnbase::readMSData(files, mode = "onDisk", msLevel. = 1)
    serial <- BiocParallel::SerialParam()
    xn <- xcms::groupChromPeaks(xcms::findChromPeaks(xn, param = cwp,
                                                     BPPARAM = serial),
                                param = pdp)
    xe <- MsExperiment::readMsExperiment(files, BPPARAM = serial)
    xe <- xcms::findChromPeaks(xe, param = cwp, msLevel = 1L,
                               BPPARAM = serial)
    xe <- xcms::groupChromPeaks(xe, param = pdp, msLevel = 1L)
    expect_s4_class(xe, "XcmsExperiment")
    pa <- .golden_pa("1_purityA_pa.rds")
    pn <- suppressMessages(frag4feature(pa, xn))
    pe <- suppressMessages(frag4feature(pa, xe))
    expect_identical(pe@grped_df, pn@grped_df)
    td <- tempfile("xe-")
    dir.create(td)
    dn <- suppressWarnings(suppressMessages(createDatabase(
        pn, xn, outDir = td, dbName = "n.sqlite")))
    de <- suppressWarnings(suppressMessages(createDatabase(
        pe, xe, outDir = td, dbName = "e.sqlite")))
    expect_identical(.sqlite_dump(de), .sqlite_dump(dn))
    skip_if_no_arrow()
    p <- suppressMessages(createDatabase(pe, xe, outDir = td, dbName = "e",
                                         format = "parquet"))
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("both routes write the same tables from the same results", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    td <- tempfile("routes-")
    dir.create(td)
    db <- suppressWarnings(suppressMessages(createDatabase(
        .golden_pa(), .golden_xcms("msms_only_xcmsnexp.rds"), outDir = td,
        dbName = "x.sqlite")))
    conv <- file.path(td, "conv")
    convertSqliteToParquet(db, conv)
    live <- .results_dataset()
    for (t in c("feature", "chromatographic_peak",
                "chromatographic_peak_feature", "abundance", "assay",
                "sample", "x_mspurity_scan")) {
        a <- .read_table(conv, t)
        b <- .read_table(live, t)
        a$activity <- b$activity <- NULL
        expect_identical(b[order(b[[1]]), names(a)], a[order(a[[1]]), ],
                         info = t)
    }
})
