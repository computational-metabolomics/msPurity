# Helpers for the Parquet backend tests.

skip_if_no_arrow <- function() {
    for (p in c("arrow", "jsonlite"))
        testthat::skip_if_not_installed(p)
    testthat::skip_if_not(arrow::codec_is_available("zstd") &&
                          arrow::codec_is_available("snappy"),
                          "arrow lacks the zstd or snappy codec")
}

skip_if_no_mspuritydata <- function() {
    testthat::skip_if_not_installed("msPurityData")
}

.mzs <- function(name) get(name, envir = asNamespace("msPurity"))

# Run `expr`, returning its value and the deprecation warnings it raised.
.with_deprecations <- function(expr) {
    n <- 0L
    msgs <- character()
    value <- withCallingHandlers(
        expr,
        deprecatedWarning = function(w) {
            n <<- n + 1L
            msgs <<- c(msgs, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
    list(value = value, n = n, messages = msgs)
}

# Expect `expr` to raise a Parquet error of the given class.
expect_parquet_error <- function(expr, class) {
    testthat::expect_error(expr, class = paste0("msPurity_parquet_", class))
}

# A fresh directory path, not yet created.
.fresh_path <- function(prefix = "mzs-") {
    tempfile(prefix)
}

# A study dataset of the two LC-MS/MS files, converted once per session with
# MsBackendParquet. Returns its path; tests copy it before changing it.
.study_dataset <- local({
    path <- NULL
    function() {
        testthat::skip_if_not_installed("MsBackendParquet")
        skip_if_no_mspuritydata()
        if (is.null(path) || !dir.exists(path)) {
            p <- tempfile("study-")
            suppressMessages(MsBackendParquet::mzMLToParquet(
                unname(.lcmsms_paths()), path = p))
            path <<- p
        }
        path
    }
})

# A copy of a dataset directory.
.copy_dataset <- function(from) {
    to <- tempfile("copy-")
    dir.create(to)
    file.copy(list.files(from, full.names = TRUE, all.files = TRUE,
                         no.. = TRUE), to, recursive = TRUE)
    to
}

# The results dataset createDatabase(format = "parquet") writes for the
# standard fixtures (purityA after averaging, the XCMSnExp), with or without
# the study. Written once per session; tests copy it before changing it.
.results_dataset <- local({
    cache <- list()
    function(study = FALSE, stage = "9_averageAllFragSpectra_with_filter_pa.rds") {
        skip_if_no_arrow()
        skip_if_no_mspuritydata()
        testthat::skip_if_not_installed("xcms")
        key <- paste(study, stage)
        if (is.null(cache[[key]]) || !dir.exists(cache[[key]])) {
            suppressPackageStartupMessages(library(xcms))
            out <- tempfile("results-")
            dir.create(out)
            p <- suppressMessages(createDatabase(
                .golden_pa(stage), .golden_xcms("msms_only_xcmsnexp.rds"),
                outDir = out, dbName = "results.parquet", format = "parquet",
                study = if (study) .study_dataset()))
            cache[[key]] <<- p
        }
        cache[[key]]
    }
})

# A results table as a data.frame.
.read_table <- function(path, name) {
    m <- jsonlite::fromJSON(file.path(path, "mzStack.json"),
                            simplifyVector = FALSE)
    e <- m$results$tables[[name]]
    if (is.null(e) || isTRUE(e$omitted))
        return(NULL)
    ds <- arrow::open_dataset(file.path(path, e$path), format = "parquet")
    df <- as.data.frame(dplyr::collect(ds))
    for (c in names(df)) if (is.factor(df[[c]])) df[[c]] <- as.character(df[[c]])
    df
}

# A native run's spectra as a data.frame (list columns kept).
.read_run <- function(path, run_id) {
    as.data.frame(arrow::read_parquet(file.path(
        path, "spectra", paste0("run_id=", run_id), "part-0.parquet")))
}

.manifest <- function(path) {
    jsonlite::fromJSON(file.path(path, "mzStack.json"), simplifyVector = FALSE)
}
