# Canonical dumps of the SQLite databases msPurity writes.
#
# The SQLite output is frozen while it is deprecated: the same inputs must
# keep producing the same tables, column types and values. Byte equality is
# not usable, because the `source` table records a timestamp and the
# package version, and paths depend on where msPurityData is installed. A
# dump records, per table, its CREATE statement and every row in rowid
# order, with those machine- and time-dependent values masked.

.golden_mask <- function(x, roots) {
    if (!is.character(x))
        return(x)
    for (nm in names(roots))
        x <- gsub(roots[[nm]], nm, x, fixed = TRUE)
    x
}

.golden_roots <- function() {
    roots <- c(
        "<msPurityData>" = system.file(package = "msPurityData"),
        "<msPurity>" = system.file(package = "msPurity"),
        "<tmp>" = normalizePath(tempdir(), mustWork = FALSE),
        "<tmp>" = tempdir())
    roots[nzchar(roots)]
}

.sqlite_dump <- function(db) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    on.exit(DBI::dbDisconnect(con))
    master <- DBI::dbGetQuery(
        con, "SELECT name, sql FROM sqlite_master WHERE type = 'table' ORDER BY name")
    roots <- .golden_roots()
    out <- lapply(seq_len(nrow(master)), function(i) {
        nm <- master$name[i]
        # Columns typed by their first row hold mixed types, which RSQLite
        # reports on every read; the values are what is compared.
        rows <- suppressWarnings(DBI::dbGetQuery(
            con, sprintf('SELECT * FROM "%s" ORDER BY rowid', nm)))
        if (nm == "source") {
            for (c in intersect(c("name", "parsing_software"), names(rows)))
                rows[[c]] <- sub("[0-9][-0-9:. ]*[0-9]$", "<masked>", rows[[c]])
        }
        rows[] <- lapply(rows, .golden_mask, roots = roots)
        list(sql = master$sql[i], rows = rows)
    })
    names(out) <- master$name
    out
}

.golden_dir <- function() {
    d <- testthat::test_path("fixtures", "sqlite-golden")
    if (dir.exists(d)) d else
        file.path("tests", "testthat", "fixtures", "sqlite-golden")
}

.lcmsms_paths <- function() {
    p <- list.files(system.file("extdata", "lcms", "mzML",
                                package = "msPurityData"),
                    full.names = TRUE, pattern = "MSMS")
    c(LCMSMS_1 = p[basename(p) == "LCMSMS_1.mzML"],
      LCMSMS_2 = p[basename(p) == "LCMSMS_2.mzML"])
}

.golden_pa <- function(stage = "9_averageAllFragSpectra_with_filter_pa.rds") {
    pa <- readRDS(system.file("extdata", "tests", "purityA", stage,
                              package = "msPurity"))
    p <- .lcmsms_paths()
    pa@fileList[1] <- p[["LCMSMS_1"]]
    pa@fileList[2] <- p[["LCMSMS_2"]]
    pa
}

.golden_xcms <- function(fn) {
    x <- readRDS(system.file("extdata", "tests", "xcms", fn,
                             package = "msPurity"))
    p <- .lcmsms_paths()
    if (is(x, "XCMSnExp")) {
        x@phenoData@data[1, ] <- p[["LCMSMS_1"]]
        x@phenoData@data[2, ] <- p[["LCMSMS_2"]]
        x@processingData@files[1] <- p[["LCMSMS_1"]]
        x@processingData@files[2] <- p[["LCMSMS_2"]]
    } else if (is(x, "xcmsSet")) {
        x@filepaths[1] <- p[["LCMSMS_1"]]
        x@filepaths[2] <- p[["LCMSMS_2"]]
    }
    x
}

# Each writer the frozen SQLite output must keep reproducing. Every element
# is a function of a fresh output directory returning the database written.
.golden_writers <- function() {
    list(
        createDatabase_xcmsnexp = function(td) {
            createDatabase(pa = .golden_pa(),
                           xcmsObj = .golden_xcms("msms_only_xcmsnexp.rds"),
                           outDir = td, dbName = "out.sqlite")
        },
        createDatabase_xset = function(td) {
            createDatabase(pa = .golden_pa(),
                           xcmsObj = .golden_xcms("msms_only_xset.rds"),
                           outDir = td, dbName = "out.sqlite")
        },
        frag4feature_createDb = function(td) {
            pa <- frag4feature(pa = .golden_pa("1_purityA_pa.rds"),
                               xcmsObj = .golden_xcms("msms_only_xcmsnexp.rds"),
                               createDb = TRUE, outDir = td,
                               dbName = "out.sqlite")
            pa@db_path
        },
        create_database = function(td) {
            pa <- .golden_pa("9_averageAllFragSpectra_with_filter_pa_OLD.rds")
            create_database(pa, .golden_xcms("msms_only_xset_OLD.rds"),
                            out_dir = td, db_name = "out.sqlite")
        },
        spectralMatching_qvq = function(td) {
            q <- system.file("extdata", "tests", "db",
                             "createDatabase_example.sqlite",
                             package = "msPurity")
            out <- file.path(td, "out.sqlite")
            spectralMatching(q_dbPth = q, l_dbPth = q,
                             q_xcmsGroups = c(89, 410),
                             q_spectraTypes = "av_all", q_spectraFilter = TRUE,
                             q_pol = NA, l_xcmsGroups = c(89, 410),
                             l_spectraTypes = "av_all", l_pol = NA,
                             l_spectraFilter = TRUE, cores = 1,
                             updateDb = TRUE, copyDb = TRUE,
                             usePrecursors = TRUE, outPth = out)
            out
        },
        combineAnnotations = function(td) {
            ext <- function(...) system.file("extdata", "tests", ...,
                                             package = "msPurity")
            out <- file.path(td, "out.sqlite")
            file.copy(ext("sm", "spectralMatching_result.sqlite"), out)
            combineAnnotations(
                out, compoundDbPth = ext("db", "metab_compound_subset.sqlite"),
                metfrag_resultPth = ext("external_annotations", "metfrag.tsv"),
                sirius_csi_resultPth = ext("external_annotations",
                                           "sirus_csifingerid.tsv"),
                probmetab_resultPth = ext("external_annotations",
                                          "probmetab.tsv"),
                ms1_lookup_resultPth = ext("external_annotations", "beams.tsv"),
                weights = list(sm = 0.3, metfrag = 0.2,
                               sirius_csifingerid = 0.2, probmetab = 0,
                               ms1_lookup = 0.05, biosim = 0.25),
                ms1_lookup_dbSource = "hmdb",
                ms1_lookup_checkAdducts = FALSE,
                ms1_lookup_keepAdducts = c("[M+H]+", "[M-H]-"))
            out
        },
        purityX_saveEIC = function(td) {
            out <- file.path(td, "out.sqlite")
            wd <- setwd(td)
            on.exit(setwd(wd))
            purityX(.golden_xcms("msms_only_xset_OLD.rds"), saveEIC = TRUE,
                    sqlitePth = out, plotP = FALSE, xgroups = c(1, 2, 3))
            out
        })
}

# Run one writer in a fresh directory and dump what it wrote.
.golden_run <- function(name) {
    td <- tempfile(paste0("golden-", name, "-"))
    dir.create(td)
    db <- suppressWarnings(suppressMessages(
        utils::capture.output(res <- .golden_writers()[[name]](td))))
    .sqlite_dump(res)
}
