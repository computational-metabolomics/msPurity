context("Purity spectrum positions and export")

.seqnum_file <- function() {
    system.file("extdata", "msn-standards",
        "msPurity_L_tryptophan_6305_MSn.mzML.gz", package = "msPurity",
        mustWork = TRUE)
}

test_that("repeated metadata selectors return one unchanged spectrum field", {
    tab <- data.frame(pid = 1L, fileid = 1L, seqNum = 19L,
        retentionTime = 10, precursorMZ = 205, precursorScanNum = 1L)
    pa <- methods::new("purityA")
    pa@spectra <- msPurity:::.pa_spectra_from_table(tab, "example.mzML")
    pa@spectra@metadata$puritydfColumns <- c(names(tab), "seqNum")
    path <- tempfile(fileext = ".rds")
    saveRDS(pa, path)
    reopened <- readRDS(path)
    before <- serialize(reopened, NULL)

    out <- purityTable(reopened)

    expect_identical(names(out), names(tab))
    expect_equal(out, tab)
    expect_identical(serialize(reopened, NULL), before)
})

test_that("purity construction refuses reordered initial spectrum keys", {
    original <- plyr::ddply
    local_mocked_bindings(ddply = function(.data, .variables, .fun, ...) {
        out <- original(.data, .variables, .fun, ...)
        if (identical(.fun, msPurity:::get_init_purity))
            out$seqNum <- rev(out$seqNum)
        out
    }, .package = "plyr")

    expect_error(assessPuritySingle(.seqnum_file(), cores = 1L,
            plotP = FALSE),
        "Initial purity results do not match the spectrum positions",
        fixed = TRUE)
})
