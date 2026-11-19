test_that("createMSP writes grouped m/z and RT for scans and averages", {
    skip_if_not_installed("xcms")
    prepared <- .msn_standard_feature("tryptophan", "MSn")
    feature <- xcms::featureDefinitions(prepared$xe)[prepared$feature, ]
    pa <- purityA(.msn_standard_file("tryptophan", "MSn"))
    pa <- frag4feature(pa, prepared$xe, ppm = 10, useGroup = TRUE)
    pa <- averageIntraFragSpectra(pa)
    pa <- averageInterFragSpectra(pa)
    pa <- averageAllFragSpectra(pa)

    file <- tempfile(fileext = ".msp")
    expected <- c(paste("PRECURSORMZ:", feature$mzmed),
                  paste("RETENTIONTIME:", feature$rtmed))
    records <- c(all = 17L, max = 1L, av_intra = 1L, av_inter = 1L,
                 av_all = 1L)
    for (method in names(records)) {
        createMSP(pa, file, xcms_groupids = prepared$feature, method = method,
                  msp_schema = "mona", filter = FALSE)
        headers <- grep("^(PRECURSORMZ|RETENTIONTIME):", readLines(file),
                        value = TRUE)
        expect_equal(headers, rep(expected, records[[method]]), label = method)
    }
})

test_that("createMSP closes its output connection after a write error", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
        "9_averageAllFragSpectra_with_filter_pa.rds", package = "msPurity"))
    connection <- NULL
    local_mocked_bindings(write.msp = function(precmz, rtmed, grpid, fileid,
            spectra, metadata, metadata_cols, ofile, ...) {
        connection <<- ofile
        stop("MSP write interrupted")
    }, .package = "msPurity")
    expect_error(createMSP(pa, tempfile(), xcms_groupids = 187, method = "max"),
                 "MSP write interrupted")
    expect_false(as.integer(connection) %in% getAllConnections())
})
