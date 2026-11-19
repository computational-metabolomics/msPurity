# Canonical dumps of the Parquet results datasets msPurity writes.
#
# A dump holds the manifest, every results table, every native run and the
# other JSON files of a dataset, with the values that change between runs
# masked: the dataset uid, timestamps, the package and R versions, and the
# machine's paths. Two datasets are equal when their dumps are identical.

.parquet_mask <- function(x, uid) {
    if (inherits(x, "POSIXt") || inherits(x, "Date"))
        return(rep(NA_character_, length(x)))
    if (!is.character(x))
        return(x)
    if (length(uid) && !is.na(uid))
        x <- gsub(uid, "<uid>", x, fixed = TRUE)
    x <- gsub("[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9:.]+Z", "<time>", x)
    x <- .golden_mask(x, .golden_roots())
    x <- gsub("(study|results|copy|golden|mzs)-[0-9a-f]+", "\\1-<id>", x)
    x <- gsub("msPurity::createDatabase [0-9.]+", "msPurity::createDatabase <version>", x)
    x
}

.parquet_canon <- function(path) {
    m <- jsonlite::fromJSON(file.path(path, "mzStack.json"),
                            simplifyVector = FALSE)
    uid <- m$uid
    mask <- function(x) .parquet_mask(x, uid)
    for (i in seq_along(m$provenance)) {
        m$provenance[[i]]$agent$version <- "<version>"
        m$provenance[[i]]$environment <- NULL
    }
    m <- rapply(m, mask, how = "replace", classes = "character")
    tables <- lapply(names(m$results$tables), function(t) {
        df <- .read_table(path, t)
        if (!is.null(df)) df[] <- lapply(df, mask)
        df
    })
    names(tables) <- names(m$results$tables)
    runs <- lapply(m$runs, function(r) {
        df <- .read_run(path, r$run_id)
        df[] <- lapply(df, function(v) if (is.list(v)) v else mask(v))
        df
    })
    names(runs) <- vapply(m$runs, `[[`, "", "run_id")
    json <- setdiff(list.files(path, "\\.json$", recursive = TRUE),
                    "mzStack.json")
    other <- lapply(json, function(f) mask(readLines(file.path(path, f),
                                                      warn = FALSE)))
    names(other) <- json
    list(manifest = m, tables = tables, runs = runs, json = other)
}

.parquet_golden_dir <- function() {
    d <- testthat::test_path("fixtures", "parquet-golden")
    if (dir.exists(d)) d else
        file.path("tests", "testthat", "fixtures", "parquet-golden")
}

# A purityA of the D02 tryptophan MSn acquisition processed through the MS2
# workflow, and its XCMS object.
.d02_ms2_workflow <- local({
    cache <- NULL
    function() {
        if (!is.null(cache)) return(cache)
        serial <- BiocParallel::SerialParam()
        f <- .msn_standard_file("tryptophan", "MSn")
        xe <- xcms::findChromPeaks(
            MsExperiment::readMsExperiment(f, BPPARAM = serial),
            param = xcms::CentWaveParam(ppm = 10, mzdiff = 0.001,
                peakwidth = c(3, 30), snthresh = 10, prefilter = c(3, 100),
                fitgauss = FALSE, noise = 0, verboseColumns = FALSE,
                firstBaselineCheck = TRUE), msLevel = 1L, BPPARAM = serial)
        xe <- xcms::groupChromPeaks(xe, param = xcms::PeakDensityParam(
            sampleGroups = 1L, binSize = 0.005, bw = 0.5, minFraction = 0.5,
            maxFeatures = 50))
        pa <- purityA(f, cores = 1L, plotP = FALSE)
        pa <- frag4feature(pa, xe, ppm = 10, convert2RawRT = FALSE)
        pa <- filterFragSpectra(pa)
        pa <- averageAllFragSpectra(pa, cores = 1L)
        cache <<- list(pa = pa, xe = xe)
        cache
    }
})

# Parquet writers whose output is frozen in fixtures/parquet-golden. Each is a
# function of a fresh output directory returning the dataset written.
.parquet_golden_writers <- function() {
    xe <- function() .golden_xcms("msms_only_xcmsnexp.rds")
    write <- function(pa, x, td, ...)
        createDatabase(pa, x, outDir = td, dbName = "out.parquet",
                       format = "parquet", ...)
    list(
        xcmsnexp = function(td) write(.golden_pa(), xe(), td),
        xset = function(td) write(.golden_pa(),
                                  .golden_xcms("msms_only_xset.rds"), td),
        study = function(td) write(.golden_pa(), xe(), td,
                                   study = .study_dataset()),
        metadata = function(td) write(.golden_pa(), xe(), td,
            metadata = list(polarity = "positive", instrument = "Q-Exactive",
                            instrument_type = "Orbitrap")),
        frag4feature = function(td)
            write(.golden_pa("2_frag4feature_pa.rds"), xe(), td),
        filterFragSpectra = function(td)
            write(.golden_pa("3_filterFragSpectra_pa.rds"), xe(), td),
        old = function(td)
            write(.golden_pa("9_averageAllFragSpectra_with_filter_pa_OLD.rds"),
                  .golden_xcms("msms_only_xset_OLD.rds"), td),
        useGroup = function(td) {
            pa <- frag4feature(.golden_pa("1_purityA_pa.rds"), xe(),
                               useGroup = TRUE)
            write(pa, xe(), td)
        },
        sqlite = function(td) {
            db <- system.file("extdata", "tests", "db",
                              "createDatabase_example.sqlite",
                              package = "msPurity")
            .without_identity_warning(convertSqliteToParquet(
                db, file.path(td, "out.parquet")))
        },
        d02_ms2 = function(td) {
            w <- .d02_ms2_workflow()
            write(w$pa, w$xe, td)
        })
}

# Run one writer in a fresh directory and dump what it wrote.
.parquet_golden_run <- function(name, writers = .parquet_golden_writers()) {
    td <- tempfile(paste0("golden-", name, "-"))
    dir.create(td)
    utils::capture.output(res <- suppressWarnings(suppressMessages(
        writers[[name]](td))))
    .parquet_canon(res)
}

# Whether JSON text `got` holds every value of JSON text `ref`: an object may
# gain keys, all else must be equal. Processing steps gain parameters.
.json_contains <- function(got, ref) {
    if (is.na(ref) || is.na(got)) return(identical(is.na(got), is.na(ref)))
    within <- function(g, r) {
        if (!is.list(r)) return(identical(g, r))
        if (!is.list(g)) return(FALSE)
        if (is.null(names(r)))
            return(length(g) == length(r) && all(mapply(within, g, r)))
        all(names(r) %in% names(g)) &&
            all(vapply(names(r), function(k) within(g[[k]], r[[k]]), TRUE))
    }
    within(jsonlite::fromJSON(got, simplifyVector = FALSE),
           jsonlite::fromJSON(ref, simplifyVector = FALSE))
}

# Whether dump `got` holds everything frozen in `ref` unchanged: every native
# run, the manifest's runs and sources, every table, column and row, with
# rows of the loss ledger and source identifiers kept as a subset (their ids
# renumbered). Added tables, columns, ledger rows and step parameters are
# allowed; the versions of msPurity and of the results layout, and the
# manifest's table bookkeeping, may change.
.parquet_golden_contains <- function(got, ref) {
    msg <- character()
    fail <- function(...) msg <<- c(msg, paste0(...))
    if (!identical(got$runs, ref$runs)) fail("native runs differ")
    for (k in c("runs", "sources"))
        if (!identical(got$manifest[[k]], ref$manifest[[k]]))
            fail("manifest ", k, " differ")
    if (!identical(got$manifest$results$runs, ref$manifest$results$runs))
        fail("manifest results runs differ")
    row_key <- function(d, drop) do.call(paste, c(lapply(
        d[setdiff(names(d), drop)], function(v) ifelse(is.na(v), "<NA>",
                                                      as.character(v))),
        sep = "\r"))
    for (t in names(ref$tables)) {
        r <- ref$tables[[t]]
        g <- got$tables[[t]]
        if (is.null(r)) next
        if (is.null(g)) { fail("table ", t, " is gone"); next }
        gone <- setdiff(names(r), names(g))
        if (length(gone)) fail(t, ": columns gone: ", paste(gone, collapse = ", "))
        if (t %in% c("loss_ledger", "source_identifier")) {
            id <- if (t == "loss_ledger") "loss_ledger_id_"
                  else "source_identifier_id_"
            missing <- setdiff(row_key(r, id), row_key(g[names(r)], id))
            if (length(missing))
                fail(t, ": ", length(missing), " row(s) gone, e.g. ",
                     gsub("\r", " | ", missing[1L]))
            next
        }
        if (nrow(g) != nrow(r)) { fail(t, ": rows ", nrow(r), " -> ", nrow(g)); next }
        cols <- intersect(names(r), names(g))
        ## The versions of msPurity and the results layout that wrote them.
        if (t == "conversion")
            cols <- setdiff(cols, c("results_version", "converter_version",
                                    "source_format_version"))
        if (t == "tool_provenance") cols <- setdiff(cols, "tool_version")
        for (c in cols) {
            same <- if (t == "tool_provenance" && c == "parameters")
                all(mapply(.json_contains, g[[c]], r[[c]]))
            else identical(g[[c]], r[[c]])
            if (!same) fail(t, ".", c, " changed")
        }
    }
    if (length(msg)) msg else TRUE
}
