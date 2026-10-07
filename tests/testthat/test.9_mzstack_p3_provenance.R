context("mzStack: provenance, converter tables and the SQLite converter")

.example_db <- function(name = "createDatabase_example.sqlite", dir = "db")
    system.file("extdata", "tests", dir, name, package = "msPurity")

.converted <- local({
    cache <- list()
    function(name = "createDatabase_example.sqlite") {
        skip_if_no_mzstack()
        if (is.null(cache[[name]]) || !dir.exists(cache[[name]])) {
            p <- tempfile("conv-")
            convertSqliteToMzstack(.example_db(name), p)
            cache[[name]] <<- p
        }
        cache[[name]]
    }
})

test_that("every row names the activity that wrote it (R-072, D-067)", {
    for (p in list(.converted(), .results_dataset())) {
        m <- .manifest(p)
        ids <- vapply(m$provenance, function(a) a$id, "")
        for (t in names(m$results$tables)) {
            df <- .read_table(p, t)
            expect_true(all(df$activity %in% ids), info = t)
        }
        expect_identical(.mzs(".mzs_schema_check")(m, "manifest.schema.json"),
                         character())
    }
})

test_that("parameters round-trip exactly (D-070)", {
    skip_if_no_mzstack()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa()
    pa@av_all_params$ppm <- 0.1 + 0.2
    pa@filter_frag_params$snr <- 0.010000000000000002
    out <- tempfile("exact-")
    dir.create(out)
    p <- suppressMessages(createDatabase(
        pa, .golden_xcms("msms_only_xcmsnexp.rds"), outDir = out,
        dbName = "x", format = "mzstack"))
    tp <- .read_table(p, "tool_provenance")
    prm <- function(step) jsonlite::fromJSON(
        tp$parameters[tp$step_name == step], simplifyVector = FALSE)
    expect_identical(prm("averageAllFragSpectra")$ppm, 0.1 + 0.2)
    expect_identical(prm("filterFragSpectra")$snr, 0.010000000000000002)
    expect_true(all(tp$parameters_completeness %in%
                    c("complete", "partial", "absent")))
    m <- .manifest(p)
    expect_identical(m$provenance[[1]]$parameters$pass_criteria$averaged,
                     .mzs(".MZS_PASS_CRITERIA")$averaged)
})

test_that("frag4feature() records its arguments for provenance", {
    skip_if_no_mzstack()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- frag4feature(.golden_pa("1_purityA_pa.rds"),
                       .golden_xcms("msms_only_xcmsnexp.rds"), ppm = 7)
    expect_identical(pa@params$frag4feature$ppm, 7)
    out <- tempfile("f4f-")
    dir.create(out)
    p <- suppressMessages(createDatabase(
        pa, .golden_xcms("msms_only_xcmsnexp.rds"), outDir = out,
        dbName = "x", format = "mzstack"))
    tp <- .read_table(p, "tool_provenance")
    expect_identical(tp$parameters_completeness[tp$step_name ==
                                                "frag4feature"], "complete")
    ## The purityA object was saved before the params slot existed.
    expect_identical(tp$parameters_completeness[tp$step_name == "purityA"],
                     "partial")
})

test_that("a database keeps no parameters, and says so (R-090)", {
    p <- .converted()
    tp <- .read_table(p, "tool_provenance")
    expect_true(all(tp$parameters_completeness == "absent"))
    ll <- .read_table(p, "loss_ledger")
    expect_true("processing parameters" %in% ll$source_construct)
    cv <- .read_table(p, "conversion")
    expect_identical(cv$route, "mspurity:sqlite-file")
    expect_identical(cv$source_checksum,
                     unname(tools::sha256sum(.example_db())))
    expect_identical(cv$source_uri, .example_db())
})

test_that("the loss ledger uses only the seven dispositions (R-083)", {
    for (p in list(.converted(), .results_dataset())) {
        ll <- .read_table(p, "loss_ledger")
        expect_true(all(ll$disposition %in% .mzs(".MZS_DISPOSITIONS")))
        ok <- !is.na(ll$affected_rows)
        expect_equal(as.numeric(ll$affected_rows[ok]),
                     round(as.numeric(ll$affected_rows[ok])))
        cov <- jsonlite::fromJSON(.read_table(p, "conversion")$coverage_manifest)
        expect_true(all(unlist(cov) %in% .mzs(".MZS_COVERAGE_STATES")))
    }
})

test_that("every wholly-null column is explained (R-084)", {
    for (p in list(.converted(), .results_dataset())) {
        m <- .manifest(p)
        ll <- .read_table(p, "loss_ledger")
        for (t in setdiff(names(m$results$tables), "loss_ledger")) {
            df <- .read_table(p, t)
            for (c in names(df)) {
                if (!all(is.na(df[[c]])) || !nrow(df))
                    next
                expect_true(paste0(t, ".", c) %in% ll$target_ref,
                            info = paste0(t, ".", c))
            }
        }
    }
})

test_that("copied values are exactly the database's (B-075)", {
    p <- .converted()
    con <- DBI::dbConnect(RSQLite::SQLite(), .example_db())
    on.exit(DBI::dbDisconnect(con))
    meta <- DBI::dbGetQuery(con, "SELECT * FROM s_peak_meta")
    peaks <- DBI::dbGetQuery(con, "SELECT * FROM s_peaks")
    sc <- .read_table(p, "x_mspurity_scan")
    si <- .read_table(p, "source_identifier")
    si <- si[si$target_table == "x_mspurity_scan", ]
    pid <- as.integer(si$source_value[match(sc$scan_annotation_id_,
                                            si$target_key)])
    scans <- meta
    k <- match(pid, scans$pid)
    expect_false(anyNA(k))
    expect_true(all(scans$spectrum_type[k] == "scan"))
    expect_identical(sc$precursor_mz, as.numeric(scans$precursorMZ[k]))
    expect_identical(sc$in_purity, as.numeric(scans$inPurity[k]))
    expect_identical(sc$a_purity, as.numeric(scans$aPurity[k]))
    ## Averaged peaks.
    sp <- .read_run(p, "av_all")
    grp <- as.integer(sub("^msPurity:av_all:", "", sp$id))
    for (i in seq_along(grp)) {
        pid <- meta$pid[meta$spectrum_type == "all" & meta$grpid == grp[i]]
        pk <- peaks[peaks$pid == pid, ]
        pk <- pk[order(pk$mz, pk$cl), ]
        expect_identical(sp$mz[[i]], as.numeric(pk$mz))
        expect_identical(sp$intensity[[i]], as.numeric(pk$i))
        expect_identical(sp$x_mspurity_ra[[i]], as.numeric(pk$ra))
    }
    ## Each derived spectrum and scan carries its source identifier.
    si <- .read_table(p, "source_identifier")
    expect_true(all(si$source_value_type %in% c("positional", "name",
                                                 "path")))
    expect_setequal(si$target_key[si$target_table == "spectra"],
                    seq_len(sum(vapply(.manifest(p)$runs,
                                       function(r) r$n_spectra, 1L))))
})

test_that("the OLD (xcmsSet) layout converts too", {
    p <- .converted("createDatabase_example_OLD.sqlite")
    v <- validateMzstack(p)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("unknown or unmapped tables are refused, not dropped (R-081)", {
    skip_if_no_mzstack()
    expect_mzstack_error(convertSqliteToMzstack(
        .example_db("create_database_example.sqlite"), tempfile("old-")),
        "unsupported")
    db <- tempfile(fileext = ".sqlite")
    file.copy(.example_db(), db)
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    DBI::dbExecute(con, "CREATE TABLE extra (a INTEGER)")
    DBI::dbDisconnect(con)
    expect_mzstack_error(convertSqliteToMzstack(db, tempfile("x-")),
                         "unsupported")
    expect_mzstack_error(convertSqliteToMzstack(.example_db("metab_compound_subset.sqlite"),
                                                tempfile("y-")), "format")
})

test_that("a basename collision needs a fileMap (R-096, R-097)", {
    skip_if_no_mzstack()
    db <- tempfile(fileext = ".sqlite")
    file.copy(.example_db(), db)
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    DBI::dbExecute(con, "UPDATE fileinfo SET filename = 'QC.mzML',
                         filepth = '/day' || fileid || '/QC.mzML'")
    DBI::dbDisconnect(con)
    expect_mzstack_error(convertSqliteToMzstack(db, tempfile("c-")),
                         "semantic")
    p <- tempfile("c-")
    convertSqliteToMzstack(db, p, fileMap = c("/day1/QC.mzML" = "QC_day1",
                                              "/day2/QC.mzML" = "QC_day2"))
    m <- .manifest(p)
    expect_setequal(vapply(m$sources, function(s) s$key, ""),
                    c("QC_day1", "QC_day2"))
    expect_identical(m$provenance[[1]]$parameters$fileMap$`/day1/QC.mzML`,
                     "QC_day1")
})

test_that("the capability set is schema-valid and declares the converter", {
    skip_if_no_mzstack()
    cap <- mzstackCapabilities()
    expect_identical(.mzs(".mzs_schema_check")(cap, "capabilities.schema.json"),
                     character())
    routes <- vapply(cap$results$converter$source_formats, function(f) f$route, "")
    expect_true(all(c("mspurity:sqlite-file", "mspurity:live-object") %in% routes))
})
