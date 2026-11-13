# Run existing tests in fresh processes on one runner. The default backend
# runs first, so this is a diagnostic comparison, not an order-balanced benchmark.
test_filter <- "^test[.](0_(parity|spectra_io)|9_parquet_p8_spectra)$"
expected_files <- c("test.0_parity.R", "test.0_spectra_io.R",
                    "test.9_parquet_p8_spectra.R")

time_backend <- function(backend, output) {
    stopifnot(backend %in% c("default", "serial"))
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    Sys.setenv(R_TESTS = "", TESTTHAT_IS_CHECKING = "true")
    if (backend == "serial")
        BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)

    library(testthat)
    library(msPurity)
    param <- BiocParallel::bpparam()
    if (backend == "serial")
        stopifnot(is(param, "SerialParam"))
    writeLines(c(paste("Mode:", backend), capture.output(show(param)),
                 paste("Commit:", Sys.getenv("GITHUB_SHA")),
                 paste("Filter:", test_filter)),
               file.path(output, "backend.txt"))
    write.csv(installed.packages()[, c("Package", "Version", "Built")],
              file.path(output, "packages.csv"), row.names = FALSE)
    capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))

    started <- proc.time()
    results <- testthat::test_package("msPurity", filter = test_filter,
        reporter = "summary", stop_on_failure = FALSE, stop_on_warning = FALSE)
    elapsed <- (proc.time() - started)[["elapsed"]]
    saveRDS(results, file.path(output, "results.rds"))
    timings <- as.data.frame(results)
    timings$result <- NULL
    write.csv(timings, file.path(output, "tests.csv"), row.names = FALSE)
    write.csv(data.frame(backend = backend, test_elapsed = elapsed,
                        tests = nrow(timings), passed = sum(timings$passed),
                        failed = sum(timings$failed), errors = sum(timings$error),
                        skipped = sum(timings$skipped),
                        warnings = sum(timings$warning)),
              file.path(output, "summary.csv"), row.names = FALSE)
    stopifnot(setequal(timings$file, expected_files),
              nrow(timings) > 0L, all(timings$failed == 0L),
              !any(timings$error), !any(timings$skipped),
              all(timings$warning == 0L))
}

compare_backends <- function(output) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    output <- normalizePath(output)
    script <- normalizePath(".github/scripts/time-windows-tests.R")
    modes <- c("default", "serial")
    runs <- lapply(modes, function(backend) {
        log <- file.path(output, paste0(backend, ".log"))
        started <- proc.time()
        run <- tryCatch(callr::rscript(script,
            cmdargs = c(backend, file.path(output, backend)),
            stdout = log, stderr = "2>&1", show = FALSE,
            system_profile = FALSE, user_profile = FALSE,
            fail_on_status = FALSE, timeout = 2400),
            error = function(e) {
                message(backend, ": ", conditionMessage(e))
                list(status = 1L)
            })
        elapsed <- (proc.time() - started)[["elapsed"]]
        lines <- if (file.exists(log)) readLines(log, warn = FALSE) else character()
        cat(lines, sep = "\n")
        data.frame(backend = backend, process_elapsed = elapsed,
                   status = run$status,
                   arrow_startup_messages = sum(grepl(
                       "The tzdb package is not installed.", lines, fixed = TRUE)))
    })
    runs <- do.call(rbind, runs)
    write.csv(runs, file.path(output, "processes.csv"), row.names = FALSE)

    paths <- file.path(output, modes, "tests.csv")
    if (all(file.exists(paths))) {
        tables <- lapply(paths, read.csv, stringsAsFactors = FALSE)
        keys <- c("file", "context", "test")
        comparison <- merge(tables[[1]], tables[[2]], by = keys, all = TRUE,
                            suffixes = c("_default", "_serial"), sort = FALSE)
        comparison$elapsed_saved <- comparison$real_default - comparison$real_serial
        comparison <- comparison[order(comparison$elapsed_saved, decreasing = TRUE), ]
        write.csv(comparison, file.path(output, "comparison.csv"), row.names = FALSE)
        print(comparison[, c("file", "test", "real_default", "real_serial",
                             "elapsed_saved")], row.names = FALSE)
        # Both modes must run the same cases and expectations, including skips.
        counts <- c("nb", "failed", "skipped", "error", "warning", "passed")
        stopifnot(!anyDuplicated(tables[[1]][keys]),
                  !anyDuplicated(tables[[2]][keys]))
        for (column in counts)
            stopifnot(identical(comparison[[paste0(column, "_default")]],
                                comparison[[paste0(column, "_serial")]]))
    }
    stopifnot(all(runs$status == 0L), all(file.exists(paths)))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 2L) {
    time_backend(args[1], args[2])
} else if (!length(args)) {
    output <- Sys.getenv("TIMING_OUTPUT")
    if (!nzchar(output))
        stop("Set TIMING_OUTPUT to the diagnostic output directory.")
    compare_backends(output)
} else {
    stop("Usage: time-windows-tests.R [default|serial OUTPUT_DIRECTORY]")
}
