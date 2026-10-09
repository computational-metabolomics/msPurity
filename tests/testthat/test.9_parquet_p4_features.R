context("Parquet: features, chromatographic peaks, abundance and links; XcmsExperiment")

test_that("the feature layer is complete and referentially sound", {
    p <- .results_dataset()
    tab <- function(n) .read_table(p, n)
    ft <- tab("feature")
    pk <- tab("chromatographic_peak")
    cpf <- tab("chromatographic_peak_feature")
    ab <- tab("abundance")
    as <- tab("assay")
    sf <- tab("spectrum_feature")
    expect_identical(nrow(ft), 521L)
    expect_identical(nrow(pk), 800L)
    expect_identical(nrow(cpf), 799L)
    expect_true(all(cpf$feature_id_ %in% ft$feature_id_))
    expect_true(all(cpf$chromatographic_peak_id_ %in% pk$chromatographic_peak_id_))
    expect_true(all(ab$feature_id_ %in% ft$feature_id_))
    expect_true(all(ab$assay_id_ %in% as$assay_id_))
    ## Dense: one row per feature and assay (R-043).
    expect_identical(nrow(ab), nrow(ft) * nrow(as))
    expect_identical(.manifest(p)$results$tables$abundance$fill_policy, "dense")
    expect_true(all(pk$run_id %in% as$run_id))
    expect_true(all(sf$feature_id_ %in% ft$feature_id_))
    ## R-045: chromatographic_peak links are run-scoped, carry their peak.
    lk <- sf[sf$link_mode == "chromatographic_peak", ]
    expect_identical(nrow(lk), 75L)
    expect_true(all(lk$run_scoped))
    expect_false(anyNA(lk$chromatographic_peak_id_))
    expect_true(all(pk$run_id[match(lk$chromatographic_peak_id_,
                                    pk$chromatographic_peak_id_)] ==
                    lk$ms2_run_id))
    ## Averaged spectra are associated with their feature.
    self <- sf[sf$ms2_source == "self", ]
    expect_identical(nrow(self), 89L)
    expect_true(all(self$link_mode == "manual"))
    ## The live route records what each scan was matched at.
    expect_false(all(is.na(lk$precursor_mz_error)))
    expect_false(all(is.na(lk$x_mspurity_precursor_mz_error_exact)))
})

test_that("features keep their values exactly, and abundances their samples (R-098)", {
    p <- .results_dataset()
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    fd <- xcms::featureDefinitions(x)
    fv <- xcms::featureValues(x)
    ft <- .read_table(p, "feature")
    ft <- ft[order(ft$feature_id_), ]
    expect_identical(ft$exp_mass_to_charge, unname(fd$mzmed))
    expect_identical(ft$retention_time_in_seconds, unname(fd$rtmed))
    expect_identical(ft$feature_name, rownames(fd))
    ab <- .read_table(p, "abundance")
    for (j in 1:2) {
        a <- ab[ab$assay_id_ == j, ]
        a <- a[order(a$feature_id_), ]
        expect_identical(a$value, unname(fv[, j]))
    }
})

test_that("SQLite -> Parquet -> SQLite is exact under semantic equality", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    td <- tempfile("rt-")
    dir.create(td)
    dbs <- c(system.file("extdata", "tests", "db",
                         c("createDatabase_example.sqlite",
                           "createDatabase_example_OLD.sqlite"),
                         package = "msPurity"),
             suppressWarnings(suppressMessages(createDatabase(
                 .golden_pa(), .golden_xcms("msms_only_xset.rds"),
                 outDir = td, dbName = "xset.sqlite"))))
    for (db in dbs) {
        p <- tempfile("rt-")
        convertSqliteToParquet(db, p)
        a <- .mzs(".mzs_canonical_sqlite")(db)
        b <- .mzs(".mzs_canonical_parquet")(p)
        expect_identical(names(a), names(b))
        for (n in names(a))
            expect_identical(b[[n]], a[[n]], info = paste(basename(db), n))
    }
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

test_that("feature-width links are not run-scoped (R-045)", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    pa <- frag4feature(.golden_pa("1_purityA_pa.rds"), x, useGroup = TRUE)
    td <- tempfile("group-")
    dir.create(td)
    p <- suppressMessages(createDatabase(pa, x, outDir = td, dbName = "g",
                                         format = "parquet"))
    sf <- .read_table(p, "spectrum_feature")
    expect_true(nrow(sf) > 0)
    expect_true(all(sf$link_mode == "feature_width"))
    expect_false(any(sf$run_scoped))
    expect_true(all(is.na(sf$chromatographic_peak_id_)))
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("a per-sample column count that disagrees with the files is refused (R-098)", {
    skip_if_no_arrow()
    db <- tempfile(fileext = ".sqlite")
    file.copy(system.file("extdata", "tests", "db",
                          "createDatabase_example.sqlite",
                          package = "msPurity"), db)
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    DBI::dbExecute(con, "ALTER TABLE c_peak_groups ADD COLUMN junk TEXT")
    DBI::dbExecute(con, "UPDATE c_peak_groups SET junk = 'x' || grpid")
    DBI::dbDisconnect(con)
    expect_parquet_error(convertSqliteToParquet(db, tempfile("bad-")),
                         "semantic")
})

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
    xn <- xcms::groupChromPeaks(xcms::findChromPeaks(xn, param = cwp),
                                param = pdp)
    xe <- MsExperiment::readMsExperiment(files)
    xe <- xcms::findChromPeaks(xe, param = cwp, msLevel = 1L)
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
