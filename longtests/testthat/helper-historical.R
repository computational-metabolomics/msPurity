.historical_fixture <- function(...) {
    system.file("extdata", "regression-v1", ..., package = "msPurityData",
                mustWork = TRUE)
}

.parity_dir <- function() .historical_fixture("parity")

.golden_dir <- function() .historical_fixture("sqlite-golden")

.parquet_golden_dir <- function() .historical_fixture("parquet-golden")

.parity_ref <- function(name) readRDS(file.path(.parity_dir(), paste0(name, ".rds")))

# Run one writer in a fresh directory and dump what it wrote.
.golden_run <- function(name) {
    td <- tempfile(paste0("golden-", name, "-"))
    dir.create(td)
    db <- suppressWarnings(suppressMessages(
        utils::capture.output(res <- .golden_writers(.historical_fixture)[[name]](td))))
    .sqlite_dump(res)
}
