library(testthat)
library(msPurity)

# R CMD check --test-dir=longtests runs from a copy of this directory.
# Direct invocation from the source checkout is useful for CI and maintainers.
suite <- if (dir.exists("longtests/testthat")) "longtests/testthat" else "testthat"
if (!nzchar(system.file("extdata", "regression-v1", "SHA256SUMS",
                        package = "msPurityData")))
    stop("The long tests require regression-v1 from msPurityData ",
         "(first available in 1.41.2).")
options(timeout = max(1200, getOption("timeout")))
Sys.setenv(R_TESTS = "")
test_dir(suite, package = "msPurity", stop_on_failure = TRUE)
