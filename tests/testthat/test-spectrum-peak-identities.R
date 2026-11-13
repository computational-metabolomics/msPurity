context("Spectrum peak identities and average membership")

.identity_pa <- function(empty = FALSE) {
    raw <- lapply(c("first", "second"), function(file) {
        d <- S4Vectors::DataFrame(msLevel = c(1L, 2L, 3L, 4L, 2L, 2L),
            scanIndex = seq_len(6), acquisitionNum = 201:206,
            rtime = as.numeric(seq_len(6)), dataOrigin = file)
        d$mz <- lapply(seq_len(6), function(i) {
            if (empty || i == 2L) numeric() else as.numeric(100L + i)
        })
        d$intensity <- lapply(d$mz, function(x) rep(20, length(x)))
        Spectra::Spectra(d, source = Spectra::MsBackendMemory())
    })
    names(raw) <- c("first", "second")
    tab <- data.frame(pid = c(101:103, 201:203),
        fileid = rep(1:2, each = 3), seqNum = rep(c(2L, 5L, 6L), 2),
        acquisitionNum = rep(c(202L, 205L, 206L), 2),
        retentionTime = rep(c(2, 5, 6), 2), precursorMZ = 205,
        inPurity = c(1, 0.5, 1, 1, 0.5, 1),
        precursorScanNum = rep(201L, 6))
    pa <- methods::new("purityA")
    pa@fileList <- c(first = "first", second = "second")
    pa@puritydf <- tab
    pa@spectra <- msPurity:::.pa_spectra_from_table(tab, pa@fileList)
    list(pa = pa, raw = raw)
}

.identity_with_legacy <- function(value, code) {
    old <- options(msPurity.legacySlots = value)
    on.exit(options(old))
    force(code)
}

test_that("raw identities survive empty and higher-level spectra", {
    fixture <- .identity_pa()
    local_mocked_bindings(.msp_read = function(files, ...) {
        fixture$raw[[as.character(files)]]
    }, .package = "msPurity")
    out <- msPurity:::getScanPeaks(fixture$pa)
    expect_equal(out$pid, c(102, 103, 202, 203))
    expect_equal(out$mz, c(105, 106, 105, 106))
    expect_equal(as.character(out$scan), c("5", "6", "5", "6"))
    expect_equal(out$fileid, c(1, 1, 2, 2))
    expect_equal(out$sid, c(3, 4, 7, 8))
    fixture$pa@puritydf <- data.frame()
    expect_equal(msPurity:::getScanPeaks(fixture$pa)$pid, out$pid)
})

test_that("empty peak inputs retain a typed schema", {
    fixture <- .identity_pa(empty = TRUE)
    local_mocked_bindings(.msp_read = function(files, ...) {
        fixture$raw[[as.character(files)]]
    }, .package = "msPurity")
    expected <- c(sid = "integer", pid = "integer", fileid = "integer",
        scan = "factor", mz = "numeric", i = "numeric", type = "character")
    inputs <- list(empty_arrays = fixture$pa,
        no_assessed_spectra = methods::new("purityA"))
    for (case in names(inputs)) {
        out <- msPurity:::getScanPeaks(inputs[[case]])
        expect_equal(nrow(out), 0L, info = case)
        expect_identical(vapply(out, class, character(1L)), expected,
            info = case)
    }
})

test_that("both allfrag filtering routes retain source peak row identifiers", {
    fixture <- .identity_pa()
    local_mocked_bindings(.msp_read = function(files, ...) {
        do.call(c, unname(fixture$raw[as.character(files)]))
    }, .package = "msPurity")
    for (legacy in c(TRUE, FALSE)) for (remove in c(TRUE, FALSE)) {
        case <- paste("legacy:", legacy, "remove:", remove)
        result <- .identity_with_legacy(legacy, {
            list(
                scalar = filterFragSpectra(fixture$pa, plim = 0.8,
                    snr = 0, allfrag = TRUE, rmp = remove),
                parameter = filterFragSpectra(fixture$pa, allfrag = TRUE,
                    rmp = remove,
                    param = FilterFragSpectraParam(plim = 0.8, snr = 0)))
        })
        peaks <- lapply(result, function(pa) {
            as.list(Spectra::peaksData(allFragSpectra(pa, legacy = FALSE),
                columns = Spectra::peaksVariables(pa@fragSpectra)))
        })
        expect_identical(peaks$scalar, peaks$parameter, info = case)
        observed <- do.call(rbind, peaks$parameter)
        expect_identical(unname(observed[, "sid"]),
            if (remove) c(4, 8) else c(3, 4, 7, 8), info = case)
    }
})

test_that("real MSn peaks retain assessed spectrum positions", {
    path <- system.file("extdata", "msn-standards",
        "msPurity_L_tryptophan_6305_MSn.mzML.gz", package = "msPurity",
        mustWork = TRUE)
    pa <- purityA(path, cores = 1L, plotP = FALSE)
    out <- msPurity:::getScanPeaks(pa)
    scans <- purityTable(pa)
    raw <- Spectra::Spectra(path, source = Spectra::MsBackendMzR())
    peaks <- Spectra::peaksData(raw)

    expect_true(all(Spectra::msLevel(raw)[scans$seqNum] == 2L))
    expect_equal(nrow(out), sum(vapply(peaks[scans$seqNum], nrow,
                integer(1L))))
    expect_identical(anyDuplicated(out$sid), 0L)
    for (i in seq_len(nrow(scans))) {
        case <- paste("pid:", scans$pid[i], "spectrum position:",
            scans$seqNum[i])
        observed <- out[out$pid == scans$pid[i], , drop = FALSE]
        expected <- peaks[[scans$seqNum[i]]]
        expect_identical(unname(observed$mz), unname(expected[, "mz"]),
            info = case)
        expect_identical(unname(observed$i), unname(expected[, "intensity"]),
            info = case)
    }
})

.membership_av <- function(type = "all", method = type) {
    data.frame(pid = 7L, spectrum_type = type, method = method,
        run_id = "av_all", spectrum_id_ = 1, native_id = "average-1",
        data_origin = "mspurity:test", retention_time = 300,
        polarity = "positive", precursor_mz = 205, inPurity = 1,
        grpid = 1L)
}

.membership_src <- function(type = c("all", "scan")) {
    list(s_peaks = data.frame(pid = c(7L, 7L), mz = c(100, 999),
            i = c(20, 500), type = type, cl = c(1L, NA_integer_)))
}

.identity_collision_db <- function() {
    db <- tempfile(fileext = ".sqlite")
    original <- system.file("extdata", "tests", "db",
        "createDatabase_example.sqlite", package = "msPurity")
    stopifnot(file.copy(original, db))
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    on.exit(DBI::dbDisconnect(con))
    pid <- DBI::dbGetQuery(con, paste("SELECT pid FROM s_peak_meta",
            "WHERE spectrum_type = 'all' ORDER BY pid LIMIT 1"))$pid
    stopifnot(length(pid) == 1L)
    peaks <- DBI::dbReadTable(con, "s_peaks")
    average <- which(peaks$pid == pid & peaks$type == "all")[1L]
    raw <- which(peaks$type == "scan")[1L]
    peaks$mz[average] <- 100
    peaks$i[average] <- 20
    peaks$mz[raw] <- 999
    peaks$i[raw] <- 500
    peaks$pid[raw] <- pid
    keep <- peaks$pid != pid | peaks$type == "scan" |
        seq_len(nrow(peaks)) == average
    DBI::dbExecute(con, "DELETE FROM s_peaks")
    DBI::dbWriteTable(con, "s_peaks", peaks[keep, ], append = TRUE)
    list(path = db, average_pid = pid)
}

.identity_exported_average <- function(dataset, source_pid) {
    manifest <- msPurity:::.mzs_read_manifest(dataset)
    ids <- msPurity:::.mzs_read_table(dataset, manifest, "source_identifier")
    id <- as.numeric(ids$target_key[ids$target_table == "spectra" &
                ids$source_column == "pid" & ids$source_value == source_pid])
    spectra <- lapply(manifest$runs, function(run) {
        as.data.frame(arrow::read_parquet(file.path(
            msPurity:::.mzs_run_dir(dataset, run), "part-0.parquet")))
    })
    selected <- lapply(spectra, function(x) {
        x[x$spectrum_id_ == id, , drop = FALSE]
    })
    list(id = id, spectra = Filter(function(x) nrow(x) > 0L, selected))
}

.identity_scan_library <- function(precursor_mz, peaks) {
    path <- tempfile(fileext = ".sqlite")
    original <- system.file("extdata", "tests", "library",
        "mini_library.sqlite", package = "msPurity")
    stopifnot(file.copy(original, path))
    con <- DBI::dbConnect(RSQLite::SQLite(), path)
    on.exit(DBI::dbDisconnect(con))
    meta <- DBI::dbGetQuery(con,
        "SELECT * FROM library_spectra_meta ORDER BY id LIMIT 1")
    meta$precursor_mz <- precursor_mz
    meta$accession <- "synthetic-original-scan"
    DBI::dbExecute(con, "DELETE FROM library_spectra")
    DBI::dbExecute(con, "DELETE FROM library_spectra_meta")
    DBI::dbWriteTable(con, "library_spectra_meta", meta, append = TRUE)
    library_peaks <- data.frame(id = seq_len(nrow(peaks)),
        mz = peaks[, 1L], i = peaks[, 2L], other = NA_character_,
        library_spectra_meta_id = meta$id)
    DBI::dbWriteTable(con, "library_spectra", library_peaks, append = TRUE)
    list(path = path, pid = meta$id)
}

test_that("raw pid collisions cannot enter averaged arrays", {
    out <- msPurity:::.mzs_derived_runs(.membership_src(), .membership_av())
    expect_equal(as.numeric(out$runs[[1L]]$peaks$mz[[1L]]), 100)
    expect_equal(as.numeric(out$runs[[1L]]$peaks$intensity[[1L]]), 20)
})

test_that("average membership rejects unknown or contradictory peak types", {
    for (type in c(NA_character_, "", "unsupported", "inter")) {
        src <- .membership_src(c(type, "scan"))
        expect_error(msPurity:::.mzs_derived_runs(src, .membership_av()),
            class = "msPurity_parquet_semantic",
            info = paste("peak type:", deparse(type)))
    }
})

test_that("stored average types retain their peak membership", {
    for (type in c("all", "intra", "inter")) {
        src <- .membership_src(c(type, "scan"))
        src$s_peaks$count <- c(2L, 999L)
        out <- msPurity:::.mzs_derived_runs(src, .membership_av(type))
        expect_equal(as.numeric(out$runs[[1L]]$peaks$mz[[1L]]), 100,
            info = type)
        expect_equal(as.integer(out$runs[[1L]]$peaks$contributor_count[[1L]]),
            2L, info = type)
    }
})

test_that("empty averages retain empty peak arrays", {
    src <- .membership_src()
    src$s_peaks <- src$s_peaks[FALSE, , drop = FALSE]
    out <- msPurity:::.mzs_derived_runs(src, .membership_av())
    expect_length(out$runs[[1L]]$peaks$mz[[1L]], 0L)
})

test_that("historical all metadata retains its intra averaging method", {
    av <- .membership_av("all", "intra")
    av$source_key <- "source"
    av$source_run <- "run"
    out <- msPurity:::.mzs_derived_runs(.membership_src(), av)
    expect_equal(as.numeric(out$runs[[1L]]$peaks$mz[[1L]]), 100)
    expect_identical(out$results_runs$av_all$method, "intra")
})

test_that("conversion uses positions, not coincident acquisition numbers", {
    src <- list(s_peak_meta = data.frame(pid = c(11L, 12L, 21L),
            fileid = c(1L, 1L, 2L), seqNum = c(2L, 4L, 2L),
            acquisitionNum = c(4L, 102L, 4L),
            spectrum_type = c("scan", NA, "scan")))
    sp <- data.frame(pid = c(1L, 2L, 3L), fileid = c(1L, 1L, 2L),
        scan = c(2L, 4L, 2L))
    expect_identical(msPurity:::.mzs_scan_peak_pid(src, sp),
        c(11L, 12L, 21L))
    src$s_peak_meta$seqNum[1L] <- NA
    expect_error(msPurity:::.mzs_scan_peak_pid(src, sp),
        class = "msPurity_parquet_semantic")
})

test_that("all-fragment flags agree with legacy slots on and off", {
    fixture <- .identity_pa()
    fixture$pa@grped_df <- data.frame(grpid = integer(), pid = integer(),
        inPurity = numeric())
    local_mocked_bindings(.msp_read = function(files, ...) {
        do.call(c, unname(fixture$raw[as.character(files)]))
    }, .package = "msPurity")
    for (legacy in c(TRUE, FALSE)) {
        .identity_with_legacy(legacy, {
            for (remove in c(TRUE, FALSE)) {
                case <- paste("legacySlots =", legacy, ", remove =", remove)
                out <- filterFragSpectra(fixture$pa, allfrag = TRUE,
                    plim = 0.8, snr = 0, rmp = remove)
                peaks <- allFragSpectra(out, legacy = TRUE)
                expect_equal(peaks$pid, if (remove) c(103, 203)
                    else c(102, 103, 202, 203), info = case)
                expect_identical(peaks$purity_pass_flag,
                    if (remove) c(TRUE, TRUE) else c(FALSE, TRUE, FALSE, TRUE),
                    info = case)
                expect_identical(peaks$pass_flag, peaks$purity_pass_flag,
                    info = case)
                if (legacy)
                    expect_identical(out@all_frag_scans, peaks, info = case)
                else expect_equal(nrow(out@all_frag_scans), 0L, info = case)
            }
        })
    }
})

test_that("both SQLite writers retain peak ownership", {
    skip_if_not_installed("msPurityData")
    writers <- .golden_writers()[c("createDatabase_xcmsnexp",
            "createDatabase_xset", "create_database")]
    for (name in names(writers)) {
        td <- tempfile("identity-writer-")
        dir.create(td)
        db <- .with_deprecations(writers[[name]](td))$value
        con <- DBI::dbConnect(RSQLite::SQLite(), db)
        peaks <- DBI::dbReadTable(con, "s_peaks")
        scans <- DBI::dbReadTable(con, "s_peak_meta")
        DBI::dbDisconnect(con)
        raw <- if ("type" %in% names(peaks))
            peaks[peaks$type == "scan", ] else peaks
        idx <- match(paste(raw$fileid, raw$scan),
            paste(scans$fileid, scans$seqNum))
        expect_false(anyNA(idx), info = name)
        expect_equal(raw$pid, scans$pid[idx], info = name)
    }
})

test_that("historical conversion warns without changing its SQLite input", {
    skip_if_not_installed("arrow")
    db <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite",
        package = "msPurity", mustWork = TRUE)
    before <- tools::md5sum(db)
    out <- tempfile("historical-identities-")
    expect_warning(convertSqliteToParquet(db, out),
        "Stored scan peak identities disagree",
        class = "msPurity_parquet_identity_warning")
    expect_identical(tools::md5sum(db), before)
})

test_that("SQLite collision conversion retains only the actual average peak", {
    skip_if_not_installed("arrow")
    fixture <- .identity_collision_db()
    expect_length(fixture$average_pid, 1L)
    out <- tempfile("collision-")
    expect_warning(convertSqliteToParquet(fixture$path, out),
        "Stored scan peak identities disagree",
        class = "msPurity_parquet_identity_warning")
    result <- .identity_exported_average(out, fixture$average_pid)
    expect_length(result$id, 1L)
    expect_length(result$spectra, 1L)
    target <- result$spectra[[1L]]
    expect_equal(as.numeric(target$mz[[1L]]), 100)
    expect_equal(as.numeric(target$intensity[[1L]]), 20)
})

test_that("missing average metadata cannot silently drop recorded peaks", {
    for (type in c(NA_character_, "", "unsupported", "scan")) {
        src <- .membership_src()
        src$s_peak_meta <- data.frame(pid = 7L, spectrum_type = type)
        expect_error(msPurity:::.mzs_derived_index(src, data.frame(), "test"),
            class = "msPurity_parquet_semantic",
            info = paste("metadata type:", deparse(type)))
    }
})

test_that("scan matching after an empty spectrum uses its original peaks", {
    skip_if_not_installed("msPurityData")
    pa <- .golden_pa()
    sp <- msPurity:::.msp_read(.lcmsms_paths()[[1L]])
    header <- msPurity:::.msp_header(sp)
    peaks <- msPurity:::.msp_peaks(sp)
    size <- vapply(peaks, nrow, integer(1L))
    empty <- which(header$msLevel == 2L & size == 0L)[1L]
    position <- which(header$msLevel == 2L & size > 0L &
            seq_along(size) > empty)[1L]
    expect_true(is.finite(position))
    scans <- purityTable(pa)
    pid <- scans$pid[scans$fileid == 1L & scans$seqNum == position]
    expect_length(pid, 1L)
    td <- tempfile("scan-matching-")
    dir.create(td)
    query <- .with_deprecations(
        .golden_writers()$createDatabase_xcmsnexp(td))$value
    original <- peaks[[position]]
    reference <- .identity_scan_library(
        scans$precursorMZ[match(pid, scans$pid)], original)
    matches <- spectralMatching(query, reference$path, q_pids = pid,
        l_pids = reference$pid, q_spectraTypes = "scan", q_pol = NA,
        l_pol = NA, q_spectraFilter = FALSE, cores = 1L)$matchedResults
    expect_equal(nrow(matches), 1L)
    expect_equal(as.numeric(matches$dpc), 1, tolerance = 1e-14)
    expect_equal(as.numeric(matches$mcount), nrow(original))
})

test_that("averaged peak ids remain valid when the raw table is empty", {
    skip_if_not_installed("msPurityData")
    local_mocked_bindings(getScanPeaks = function(pa) {
        data.frame(sid = integer(), pid = integer(), fileid = integer(),
            scan = factor(), mz = numeric(), i = numeric(), type = character())
    }, .package = "msPurity")
    td <- tempfile("empty-raw-")
    dir.create(td)
    db <- .with_deprecations(
        .golden_writers()$createDatabase_xcmsnexp(td))$value
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    peaks <- DBI::dbReadTable(con, "s_peaks")
    DBI::dbDisconnect(con)
    expect_gt(nrow(peaks), 0L)
    expect_false(anyNA(peaks$sid))
    expect_false(anyDuplicated(peaks$sid) > 0L)
    expect_equal(min(peaks$sid), 1)
})

test_that("SQLite flags use the correct spectrum with either slot setting", {
    skip_if_not_installed("msPurityData")
    for (legacy in c(TRUE, FALSE)) {
        .identity_with_legacy(legacy, {
            for (remove in c(TRUE, FALSE)) {
                case <- paste("legacySlots =", legacy, ", remove =", remove)
                pa <- filterFragSpectra(.golden_pa("2_frag4feature_pa.rds"),
                    plim = 0.7, snr = 3, allfrag = TRUE, rmp = remove)
                td <- tempfile("flagged-writer-")
                dir.create(td)
                db <- .with_deprecations(createDatabase(pa,
                        .golden_xcms("msms_only_xcmsnexp.rds"), outDir = td,
                        dbName = "out.sqlite"))$value
                con <- DBI::dbConnect(RSQLite::SQLite(), db)
                peaks <- DBI::dbReadTable(con, "s_peaks")
                scans <- DBI::dbReadTable(con, "s_peak_meta")
                DBI::dbDisconnect(con)
                idx <- match(paste(peaks$fileid, peaks$scan),
                    paste(scans$fileid, scans$seqNum))
                expect_false(anyNA(idx), info = case)
                expect_equal(peaks$pid, scans$pid[idx], info = case)
                expect_identical(as.numeric(peaks$purity_pass_flag) == 1,
                    scans$inPurity[idx] > 0.7, info = case)
                if (remove)
                    expect_true(all(peaks$pass_flag == 1), info = case)
            }
        })
    }
})

test_that("empty filtered exports do not reread or restore raw peaks", {
    skip_if_not_installed("msPurityData")
    pa <- filterFragSpectra(.golden_pa("2_frag4feature_pa.rds"),
        plim = 1, allfrag = TRUE, rmp = TRUE)
    td <- tempfile("all-removed-")
    dir.create(td)
    db <- with_mocked_bindings(
        .with_deprecations(createDatabase(pa,
                .golden_xcms("msms_only_xcmsnexp.rds"), outDir = td,
                dbName = "out.sqlite"))$value,
        .msp_peaks = function(...) stop("Raw peak arrays must not be read."),
        .package = "msPurity")
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM s_peaks")$n, 0)
    expect_gt(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM s_peak_meta")$n, 0)
    removed_fields <- DBI::dbListFields(con, "s_peaks")
    removed_types <- DBI::dbGetQuery(con, "PRAGMA table_info(s_peaks)")$type
    DBI::dbDisconnect(con)

    retained <- filterFragSpectra(.golden_pa("2_frag4feature_pa.rds"),
        plim = 1, allfrag = TRUE, rmp = FALSE)
    db <- .with_deprecations(createDatabase(retained,
            .golden_xcms("msms_only_xcmsnexp.rds"), outDir = td,
            dbName = "retained.sqlite"))$value
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    expect_identical(removed_fields, DBI::dbListFields(con, "s_peaks"))
    expect_identical(removed_types,
        DBI::dbGetQuery(con, "PRAGMA table_info(s_peaks)")$type)
    DBI::dbDisconnect(con)
})
