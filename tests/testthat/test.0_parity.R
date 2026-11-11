context("Parity with the results before the move to Spectra containers")

# The references in fixtures/parity hold the legacy-shaped results of each
# workflow before the Spectra transition. Only raw peak identities and their
# purity/overall flags were corrected to match file positions. Comparisons
# remain exact.

suppressPackageStartupMessages(library(xcms))

purityA_ref <- .parity_ref("purityA")

test_that("purityA workflows reproduce the reference results", {
    res <- .parity_purityA()
    expect_identical(names(res), names(purityA_ref))
    for (stage in names(purityA_ref))
        expect_identical(res[[stage]], purityA_ref[[stage]], label = stage)
})

test_that("the accessors rebuild the reference results from the Spectra slots", {
    res <- .parity_purityA(via = "accessors")
    for (stage in names(purityA_ref))
        expect_identical(res[[stage]], purityA_ref[[stage]], label = stage)
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
        expect_identical(res[[stage]], ref[[stage]], label = stage)
})
