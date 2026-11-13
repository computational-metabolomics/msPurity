# Measure the complete focused test, including worker startup and cleanup.
script_started <- proc.time()
output <- Sys.getenv("TIMING_OUTPUT")
if (!nzchar(output))
    stop("Set TIMING_OUTPUT to the diagnostic output directory.")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(R_TESTS = "", TESTTHAT_IS_CHECKING = "true")
BiocParallel::register(BiocParallel::SerialParam(), default = TRUE)
library(testthat)
library(msPurity)
before <- BiocParallel::bpparam()
writeLines(c(paste("Commit:", Sys.getenv("GITHUB_SHA")),
             paste("OS:", Sys.info()[["sysname"]]),
             capture.output(show(before))), file.path(output, "backend.txt"))
write.csv(installed.packages()[, c("Package", "Version", "Built")],
          file.path(output, "packages.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))

started <- proc.time()
results <- testthat::test_package("msPurity", filter = "^parallel-reader$",
    reporter = "summary", stop_on_failure = FALSE, stop_on_warning = FALSE)
elapsed <- (proc.time() - started)[["elapsed"]]
saveRDS(results, file.path(output, "results.rds"))
timings <- as.data.frame(results)
timings$result <- NULL
write.csv(timings, file.path(output, "tests.csv"), row.names = FALSE)
restored <- identical(BiocParallel::bpparam(), before)
summary <- data.frame(test_elapsed = elapsed,
    script_elapsed = (proc.time() - script_started)[["elapsed"]],
    tests = nrow(timings), passed = sum(timings$passed),
    failed = sum(timings$failed), errors = sum(timings$error),
    skipped = sum(timings$skipped), warnings = sum(timings$warning),
    backend_restored = restored)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
print(summary, row.names = FALSE)
stopifnot(restored, nrow(timings) == 1L,
    identical(unique(timings$file), "test-parallel-reader.R"),
    all(timings$failed == 0L), !any(timings$error), !any(timings$skipped),
    all(timings$warning == 0L))
