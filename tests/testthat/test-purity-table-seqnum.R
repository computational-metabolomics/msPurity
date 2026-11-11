context("Purity spectrum positions and export")

.seqnum_file <- function() {
    system.file("extdata", "msn-standards",
        "msPurity_L_tryptophan_6305_MSn.mzML.gz", package = "msPurity",
        mustWork = TRUE)
}

test_that("repeated metadata selectors return one unchanged spectrum field", {
    tab <- data.frame(pid = 1L, fileid = 1L, seqNum = 19L,
        retentionTime = 10, precursorMZ = 205, precursorScanNum = 1L)
    pa <- methods::new("purityA")
    pa@spectra <- msPurity:::.pa_spectra_from_table(tab, "example.mzML")
    pa@spectra@metadata$puritydfColumns <- c(names(tab), "seqNum")
    path <- tempfile(fileext = ".rds")
    saveRDS(pa, path)
    reopened <- readRDS(path)
    before <- serialize(reopened, NULL)

    out <- purityTable(reopened)

    expect_identical(names(out), names(tab))
    expect_equal(out, tab)
    expect_identical(serialize(reopened, NULL), before)
})

test_that("purity construction refuses reordered initial spectrum keys", {
    original <- plyr::ddply
    local_mocked_bindings(ddply = function(.data, .variables, .fun, ...) {
        out <- original(.data, .variables, .fun, ...)
        if (identical(.fun, msPurity:::get_init_purity))
            out$seqNum <- rev(out$seqNum)
        out
    }, .package = "plyr")

    expect_error(assessPuritySingle(.seqnum_file(), cores = 1L,
            plotP = FALSE),
        "Initial purity results do not match the spectrum positions",
        fixed = TRUE)
})

.seqnum_processed <- function(xe) {
    pa <- purityA(.seqnum_file(), cores = 1L, plotP = FALSE)
    pa <- frag4feature(pa, xe, ppm = 10, convert2RawRT = FALSE)
    pa <- filterFragSpectra(pa)
    averageAllFragSpectra(pa, cores = 1L)
}

test_that("fresh and reopened MSn results export without duplicate positions", {
    skip_if_no_arrow()
    skip_if_not_installed("xcms")
    old <- options(msPurity.legacySlots = TRUE)
    on.exit(options(old))
    serial <- BiocParallel::SerialParam()
    xe <- xcms::findChromPeaks(MsExperiment::readMsExperiment(.seqnum_file()),
        param = xcms::CentWaveParam(ppm = 10, mzdiff = 0.001,
            peakwidth = c(3, 30), snthresh = 10, prefilter = c(3, 100),
            fitgauss = FALSE, noise = 0, verboseColumns = FALSE,
            firstBaselineCheck = TRUE), msLevel = 1L, BPPARAM = serial)
    xe <- xcms::groupChromPeaks(xe,
        param = xcms::PeakDensityParam(sampleGroups = 1L,
            binSize = 0.005, bw = 0.5, minFraction = 0.5,
            maxFeatures = 50))
    results <- list()
    for (legacy in c(TRUE, FALSE)) {
        options(msPurity.legacySlots = legacy)
        pa <- .seqnum_processed(xe)
        case <- paste("legacy slots:", legacy)
        expect_identical(sum(names(purityTable(pa)) == "seqNum"), 1L,
            info = case)
        expect_identical(sum(pa@spectra@metadata$puritydfColumns ==
                    "seqNum"), 1L, info = case)
        historical <- pa
        historical@spectra@metadata$puritydfColumns <- c(
            historical@spectra@metadata$puritydfColumns, "seqNum")
        path <- tempfile(fileext = ".rds")
        saveRDS(historical, path)
        reopened <- readRDS(path)
        expect_identical(purityTable(reopened), purityTable(pa), info = case)
        results[[as.character(legacy)]] <- list(purity = purityTable(pa),
            grouped = groupedSpectra(pa, legacy = TRUE),
            peaks = allFragSpectra(pa, legacy = TRUE),
            averages = averagedSpectra(pa, legacy = TRUE))

        for (kind in c("fresh", "reopened")) {
            selected <- if (kind == "fresh") pa else reopened
            label <- paste(case, kind)
            output <- tempfile("seqnum-export-")
            dir.create(output)
            sqlite <- .with_deprecations(createDatabase(selected,
                    xe, outDir = output, dbName = "query.sqlite",
                    format = "sqlite"))$value
            con <- DBI::dbConnect(RSQLite::SQLite(), sqlite)
            fields <- DBI::dbListFields(con, "s_peak_meta")
            DBI::dbDisconnect(con)
            expect_identical(grep("^seqNum", fields, value = TRUE),
                "seqNum", info = label)
            parquet <- createDatabase(selected, xe,
                outDir = output, dbName = "query-parquet", format = "parquet")
            report <- validateParquet(parquet)
            expect_false(any(report$level == "MUST"), info = label)
        }
    }
    expect_identical(results[["TRUE"]], results[["FALSE"]])
})
