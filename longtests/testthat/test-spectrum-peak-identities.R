test_that("both SQLite writers retain peak ownership", {
    skip_if_not_installed("msPurityData")
    writers <- .golden_writers()[c("createDatabase_xcmsnexp",
            "createDatabase_xset", "create_database")]
    for (name in names(writers)) {
        td <- tempfile("identity-writer-")
        dir.create(td)
        db <- .with_deprecations(writers[[name]](td))$value
        con <- DBI::dbConnect(RSQLite::SQLite(), db)
        peaks <- DBI::dbReadTable(con, "s_peaks")
        scans <- DBI::dbReadTable(con, "s_peak_meta")
        DBI::dbDisconnect(con)
        raw <- if ("type" %in% names(peaks))
            peaks[peaks$type == "scan", ] else peaks
        idx <- match(paste(raw$fileid, raw$scan),
            paste(scans$fileid, scans$seqNum))
        expect_false(anyNA(idx), info = name)
        expect_equal(raw$pid, scans$pid[idx], info = name)
    }
})
