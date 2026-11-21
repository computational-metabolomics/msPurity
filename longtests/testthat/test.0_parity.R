context("Parity with the results before the move to Spectra containers")

# The references in msPurityData hold the legacy-shaped results of each
# workflow before the Spectra transition. Only raw peak identities and their
# purity/overall flags were corrected to match file positions. Calculated
# values allow floating-point rounding across platforms; identity checks
# remain exact. The frag4feature_group stage compares MS/MS times with feature
# windows on the adjusted retention time scale and uses iMz (intense = TRUE),
# so it holds the links of that comparison. grped_df is compared on the legacy
# columns; frag4feature() appends its newer columns after them.

suppressPackageStartupMessages(library(xcms))

purityA_ref <- .parity_ref("purityA")

test_that("purityA slots and accessors reproduce the reference results", {
    stages <- .parity_stages()
    expect_identical(names(stages), names(purityA_ref))
    for (via in c("slots", "accessors")) {
        res <- .parity_purityA(via, stages)
        for (stage in names(purityA_ref))
            expect_equal(res[[stage]], purityA_ref[[stage]],
                         tolerance = 1e-12, label = paste(stage, via))
    }
})

test_that("parity reference stores file-position spectrum identities", {
    ref <- purityA_ref$filter_allfrag
    afs <- ref$all_frag_scans
    p <- ref$puritydf
    true_pid <- p$pid[match(paste(afs$fileid, as.character(afs$scan)),
                            paste(as.character(p$fileid), p$seqNum))]
    # References retain the independently corrected file-position mapping.
    expect_identical(afs$pid, true_pid)
})

test_that("a scan linked to two features is shared, not dropped", {
    g <- purityA_ref$two_features$grped_df
    pid <- g$pid[as.character(g$grpid) == "187"][1]
    expect_setequal(as.character(g$grpid[g$pid == pid]), c("162", "187"))
    ms2 <- purityA_ref$two_features$grped_ms2
    expect_identical(ms2[["162"]][[length(ms2[["162"]])]][, c("mz", "i")],
                     ms2[["187"]][[1]][, c("mz", "i")])
})

test_that("the purityD workflow reproduces the reference results", {
    ref <- .parity_ref("purityD")
    res <- .parity_purityD()
    for (stage in names(ref))
        expect_equal(res[[stage]], ref[[stage]],
                     tolerance = 1e-12, label = stage)
})
