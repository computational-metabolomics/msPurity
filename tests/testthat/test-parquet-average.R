test_that("historical scalar attempts with no produced runs can be exported", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    pa <- .golden_pa("2_frag4feature_pa.rds")
    pa <- filterFragSpectra(pa, ilim = 1e30)
    pa <- averageAllFragSpectra(pa, cores = 1L)
    expect_length(pa@avSpectra, 0L)
    # Historical objects retain the scalar settings and visited feature names.
    pa@params$absentAverages <- NULL
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    td <- tempfile()
    dir.create(td)
    path <- createDatabase(pa, x, outDir = td, dbName = "absent", format = "parquet")
    expect_length(.manifest(path)$runs, 0L)
    expect_null(.read_table(path, "merge_member"))
    ledger <- .read_table(path, "loss_ledger")$source_construct
    expect_setequal(grep("^av_", ledger, value = TRUE),
        paste0("av_all: feature ", unique(as.character(pa@grped_df$grpid))))
    expect_false(any(validateParquet(path)$level == "MUST"))
})

test_that("original measurement counts respect peak paths and repeated feature links", {
    # Two links to the same scan supply parent peak 1. A third measured
    # input supplies peak 2, which never reaches the inter-file average.
    members <- data.frame(merged_spectrum_id_ = c(5, 5, 5, 6),
        member_source = c("study", "study", "study", "self"),
        member_spectrum_id_ = c(10, 10, 11, 5), member_rank = c(1, 2, 3, 1),
        members_complete = TRUE, x_mspurity_contributes = TRUE,
        x_mspurity_feature_link_index = c(21L, 22L, 23L, NA_integer_))
    edges <- data.frame(merged_spectrum_id_ = c(5, 5, 5, 6, 6),
        merged_peak_rank = c(1, 1, 2, 1, 2), member_rank = c(1, 2, 3, 1, 1),
        member_peak_rank = c(1, 1, 1, 1, 1))
    reader <- list(table = function(name, ...) {
        if (name == "merge_member") members else edges
    })
    graph <- msPurity:::.mzs_graph_reader(reader, 6)
    expect_equal(length(graph$used(5)), 3L)
    expect_equal(length(graph$used(6)), 2L)
    expect_identical(unique(members$member_spectrum_id_[graph$usedRows(6)]), 10)
    edges$member_peak_rank[4] <- NA_integer_
    partial <- msPurity:::.mzs_graph_reader(reader, 6)
    expect_true(is.na(partial$used(6)))
    expect_setequal(partial$usedRows(6), c(1L, 2L, NA_integer_))
})
