context("Parity with the results before the move to Spectra containers")

# The references in fixtures/parity hold the legacy-shaped results of each
# workflow, written before msPurity read or stored spectra with Spectra.
# Every comparison is exact.

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
    for (stage in names(purityA_ref)) {
        ref <- purityA_ref[[stage]]
        if (is.list(ref) && !is.null(ref$all_frag_scans))
            ref <- .parity_fix_allfrag(ref)
        expect_identical(res[[stage]], ref, label = stage)
    }
})

test_that("allFragSpectra() names the scan of every peak correctly", {
    ref <- purityA_ref$filter_allfrag
    afs <- ref$all_frag_scans
    p <- ref$puritydf
    true_pid <- p$pid[match(paste(afs$fileid, as.character(afs$scan)),
                            paste(as.character(p$fileid), p$seqNum))]
    # The reference table is wrong after the first scan without peaks.
    expect_true(any(afs$pid != true_pid))
    expect_identical(.parity_fix_allfrag(ref)$all_frag_scans$pid, true_pid)
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
