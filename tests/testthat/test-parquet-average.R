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
