skip_if_no_study <- function() {
    skip_if_no_arrow()
    testthat::skip_if_not_installed("mzR")
    testthat::skip_if(utils::packageVersion("MsBackendParquet") < "0.99.12",
                      "MsBackendParquet converts MSn mzML from 0.99.12")
}

# A private copy of D02 files, so a test can remove them after conversion.
.d02_copy <- function(compounds = "tryptophan", kind = "MSn") {
    dir <- tempfile("d02-")
    dir.create(dir)
    files <- vapply(compounds, function(cmp) {
        f <- file.path(dir, basename(.msn_standard_file(cmp, kind)))
        file.copy(.msn_standard_file(cmp, kind), f)
        f
    }, "")
    unname(files)
}
