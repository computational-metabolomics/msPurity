context("Averaged MS2-MS4 spectra in Parquet results")

.average_sets <- function() list(
    hcd = AverageFragSpectraParam(select = list(activationType = "HCD"),
                                  combine = "collisionEnergy", minfrac = 0,
                                  ppm = 10),
    cid = AverageFragSpectraParam(select = list(activationType = "CID"),
                                  combine = "collisionEnergy", minfrac = 0,
                                  ppm = 10))

.average_all_levels <- function(pa) {
    sets <- .average_sets()
    pa <- averageIntraFragSpectra(pa, param = sets)
    pa <- averageInterFragSpectra(pa, param = sets)
    averageAllFragSpectra(pa, param = sets)
}

# Averaged spectra of a results dataset, without their peak list columns.
.av_rows <- function(path) {
    runs <- lapply(.manifest(path)$runs, function(r) {
        d <- .read_run(path, r$run_id)
        d$run_id <- r$run_id
        d[!vapply(d, is.list, TRUE)]
    })
    do.call(msPurity:::.mzs_rbind_fill, runs)
}

test_that("averages of every level are exported with their provenance", {
    skip_if_no_arrow()
    testthat::skip_if_not_installed("mzR")
    skip_if_not_installed("xcms")
    testthat::skip_if(utils::packageVersion("MsBackendParquet") < "0.99.12")
    dir <- tempfile("d02-")
    dir.create(dir)
    f <- file.path(dir, basename(.msn_standard_file("tryptophan", "MSn")))
    file.copy(.msn_standard_file("tryptophan", "MSn"), f)
    study <- tempfile("study-")
    suppressMessages(MsBackendParquet::mzMLToParquet(f, study, verbose = FALSE))
    xe <- .msn_standard_feature("tryptophan", "MSn")$xe
    pa_mz <- .average_all_levels(.msn_workflow(f, xe))
    pa_pq <- .average_all_levels(.msn_workflow(.msn_study_spectra(study), xe))
    out <- tempfile("results-")
    dir.create(out)
    write <- function(pa, name, ...) suppressMessages(createDatabase(
        pa, xe, outDir = out, dbName = name, format = "parquet", ...))
    p_mz <- write(pa_mz, "mz.parquet", study = study)
    p_pq <- write(pa_pq, "pq.parquet")
    p_ext <- write(pa_pq, "ext.parquet", study = NA)
    for (p in c(p_mz, p_pq, p_ext))
        expect_false(any(validateParquet(p)$level == "MUST"), info = p)

    view <- averagedSpectra(pa_pq, legacy = FALSE)
    av <- .av_rows(p_pq)
    expect_identical(as.vector(table(av$ms_level)),
                     as.vector(table(Spectra::msLevel(view))))
    expect_setequal(av$x_mspurity_spectrum_key, view$spectrumId)

    ## merge_member and the peak edges are the recorded provenance.
    mm <- .read_table(p_pq, "merge_member")
    mp <- .read_table(p_pq, "x_mspurity_merge_peak")
    for (i in seq_len(nrow(av))) {
        key <- av$x_mspurity_spectrum_key[i]
        got <- mm[mm$merged_spectrum_id_ == av$spectrum_id_[i], ]
        got <- got[order(got$member_rank), ]
        ref <- averageMembers(pa_pq, key)
        expect_identical(got$member_rank, ref$rank, info = key)
        expect_identical(got$x_mspurity_contributes, ref$contributes,
                         info = key)
        measured <- !is.na(ref$linkId)
        expect_identical(got$member_source[measured] == "study",
                         rep(TRUE, sum(measured)), info = key)
        expect_identical(got$member_native_id[measured],
                         ref$member[measured], info = key)
        derived <- av$x_mspurity_spectrum_key[match(
            got$member_spectrum_id_[!measured], av$spectrum_id_)]
        expect_identical(derived, ref$member[!measured], info = key)
        edges <- mp[mp$merged_spectrum_id_ == av$spectrum_id_[i], ]
        expect_identical(nrow(edges), nrow(averagePeakMembers(pa_pq, key)),
                         info = key)
    }

    ## mzML and Parquet input give the same references; without a study the
    ## same native ids.
    same <- function(p) {
        x <- .read_table(p, "merge_member")
        x[order(x$merged_spectrum_id_, x$member_rank),
          c("merged_spectrum_id_", "member_rank", "member_source",
            "member_spectrum_id_", "member_native_id",
            "x_mspurity_contributes")]
    }
    expect_identical(same(p_mz), same(p_pq))
    expect_identical(.read_table(p_mz, "x_mspurity_merge_peak"), mp)
    ext <- same(p_ext)
    expect_true(all(is.na(ext$member_spectrum_id_[ext$member_source !=
                                                  "self"])))
    expect_identical(ext$member_native_id, same(p_pq)$member_native_id)

    ## Spectral matching takes the MS2 averages only.
    q <- msPurity:::.mzs_query_spectra(p_pq,
        MsBackendParquet::readManifest(p_pq),
        list(q_spectraTypes = c("av_all", "av_intra", "av_inter"),
             q_spectraFilter = FALSE, q_purity = NA, q_pol = NA,
             q_rtrange = NA, q_xcmsGroups = NA), NULL)
    expect_identical(sort(as.numeric(q$meta$query_spectrum_id_)),
                     sort(as.numeric(av$spectrum_id_[av$ms_level == 2L])))
    expect_true(any(av$ms_level > 2L))

    expect_error(suppressWarnings(createDatabase(pa_pq, xe, outDir = out,
        dbName = "x.sqlite")), "AverageFragSpectraParam is not supported")
})

test_that("scalar empty averages retain members and absent attempts have no spectrum", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    pa <- .golden_pa("3_filterFragSpectra_pa.rds")
    pa <- averageIntraFragSpectra(pa, minnum = 1000000L, rmp = TRUE, cores = 1L)
    pa <- averageInterFragSpectra(pa, cores = 1L)
    td <- tempfile()
    dir.create(td)
    old <- options(msPurity.parquet.batch_spectra = 5L, msPurity.parquet.batch_rows = 2048L)
    on.exit(options(old))
    path <- createDatabase(pa, x, outDir = td, dbName = "empty", format = "parquet")
    expect_false(any(validateParquet(path)$level == "MUST"))
    av <- .av_rows(path)
    expect_true(nrow(av) > 0)
    view <- pa@avSpectra
    expect_equal(nrow(av), sum(view$av_level == "av_intra" & !view$is_null))
    manifest <- .manifest(path)
    part_sizes <- function(relative) {
        files <- list.files(file.path(path, relative), pattern = "\\.parquet$",
                            recursive = TRUE, full.names = TRUE)
        vapply(files, function(f) arrow::ParquetFileReader$create(f)$num_rows, numeric(1))
    }
    expect_gt(manifest$results$tables$source_identifier$rows, 2048L)
    for (entry in manifest$results$tables)
        expect_true(all(part_sizes(entry$path) <= 2048L))
    for (run in manifest$runs) {
        expect_true(all(part_sizes(run$path) <= 5L))
        rows <- .read_run(path, run$run_id)
        expect_true(all(vapply(rows$mz, identical, logical(1), numeric())))
        expect_true(all(vapply(rows$intensity, identical, logical(1), numeric())))
    }
    mm <- .read_table(path, "merge_member")
    expect_setequal(mm$merged_spectrum_id_, av$spectrum_id_)
    expect_true("x_mspurity_contributes" %in% names(mm))
    expect_true(all(!mm$x_mspurity_contributes))
    expect_null(.read_table(path, "x_mspurity_merge_peak"))
    ledger <- .read_table(path, "loss_ledger")$source_construct
    expect_setequal(grep("^av_inter: feature ", ledger, value = TRUE),
        paste0("av_inter: feature ", unique(as.character(pa@grped_df$grpid))))
    expect_false(any(startsWith(ledger, "av_all: feature ")))
    historic <- pa
    historic@params$averageGraph <- NULL
    historic@params$absentAverages <- NULL
    historic_path <- createDatabase(historic, x, outDir = td, dbName = "historical",
                                    format = "parquet")
    recovered <- .read_table(historic_path, "merge_member")
    expect_equal(table(recovered$merged_spectrum_id_), table(mm$merged_spectrum_id_))
    expect_true(all(!recovered$members_complete & !recovered$x_mspurity_contributes))
    expect_null(.read_table(historic_path, "x_mspurity_merge_peak"))
    again <- averageIntraFragSpectra(pa, cores = 1L)
    expect_null(again@params$absentAverages$inter)
})

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
