context("Parquet: derived spectra, merge_member and createDatabase(format = 'parquet')")

test_that("averaged spectra are native runs organised by provenance scope (R-010, R-011)", {
    p <- .results_dataset()
    m <- .manifest(p)
    expect_identical(m$role, "results")
    ids <- vapply(m$runs, function(r) r$run_id, "")
    expect_identical(sum(startsWith(ids, "av_intra_")), 2L)
    expect_true(all(c("av_inter", "av_all") %in% ids))
    expect_identical(ids[length(ids)], "av_all")
    for (r in m$runs) {
        expect_identical(r$kind, "native")
        expect_identical(r$signal$layout, "list")
        rr <- m$results$runs[[r$run_id]]
        expect_false(is.null(rr))
        if (startsWith(r$run_id, "av_intra_")) {
            expect_identical(rr$scope, "run")
            expect_identical(paste0("av_intra_", rr$source_run_id), r$run_id)
        } else {
            expect_identical(rr$scope, "dataset")
        }
    }
    ## Contiguous spectrum id blocks in manifest order (D-019).
    base <- vapply(m$runs, function(r) r$uid_base, 1L)
    n <- vapply(m$runs, function(r) r$n_spectra, 1L)
    expect_identical(base, as.integer(cumsum(c(1L, head(n, -1L)))))
    expect_identical(sum(n), 89L)
})

test_that("derived spectra have a dataOrigin, and ids that are not scans (R-012, R-014)", {
    p <- .results_dataset()
    m <- .manifest(p)
    for (r in m$runs) {
        sp <- .read_run(p, r$run_id)
        expect_false(anyNA(sp$data_origin))
        if (!startsWith(r$run_id, "av_intra_"))
            expect_true(all(startsWith(sp$data_origin, "mspurity:")))
        expect_false(any(grepl("(^|\\s)scan=\\d+$", sp$id)))
        expect_false(any(grepl("^(\\S+=\\d+)( \\S+=\\d+)*$", sp$id)))
        expect_true(all(startsWith(sp$id, "msPurity:av_")))
        expect_identical(sp$spectrum_id_,
                         as.integer(r$uid_base + sp$spectrum_index))
        expect_false("pass_flag" %in% names(sp))
        expect_true(all(c("sn", "contributor_count", "x_mspurity_ra",
                          "x_mspurity_snr_pass_flag") %in% names(sp)))
        expect_identical(lengths(sp$sn), lengths(sp$mz))
        expect_true(all(sp$ms_level == 2L))
    }
})

test_that("averaged peaks equal those in the purityA object", {
    p <- .results_dataset()
    pa <- .golden_pa()
    sp <- .read_run(p, "av_all")
    grp <- sub("^msPurity:av_all:", "", sp$id)
    for (i in seq_along(grp)) {
        a <- pa@av_spectra[[grp[i]]]$av_all
        a <- a[order(a$mz, a$cl), ]
        expect_identical(sp$mz[[i]], a$mz)
        expect_identical(sp$intensity[[i]], a$i)
        expect_identical(sp$x_mspurity_ra[[i]], a$ra)
        expect_identical(sp$contributor_count[[i]], as.integer(a$count))
    }
})

test_that("every derived spectrum names the spectra it merged (R-047)", {
    p <- .results_dataset()
    mm <- .read_table(p, "merge_member")
    m <- .manifest(p)
    n <- sum(vapply(m$runs, function(r) r$n_spectra, 1L))
    expect_setequal(unique(mm$merged_spectrum_id_), seq_len(n))
    expect_true(all(mm$merged_source == "self"))
    expect_true(all(mm$members_complete))
    inter <- mm[mm$merged_run_id == "av_inter", ]
    expect_true(all(inter$member_source == "self"))
    expect_true(all(startsWith(inter$member_run_id, "av_intra_")))
    scans <- mm[mm$member_source != "self", ]
    expect_true(all(grepl("scan=\\d+$", scans$member_native_id)))
    ## External sources: no spectrum id is invented (R-026).
    expect_true(all(is.na(scans$member_spectrum_id_)))
})

test_that("with a study, scans are referenced by their spectrum ids there", {
    p <- .results_dataset(study = TRUE)
    m <- .manifest(p)
    expect_identical(m$sources[[1]]$resolution, "dataset")
    mm <- .read_table(p, "merge_member")
    scans <- mm[mm$member_source == "study", ]
    expect_false(anyNA(scans$member_spectrum_id_))
    study <- .study_dataset()
    sm <- .mzs(".mzs_read_manifest")(study)
    meta <- .mzs(".mzs_read_spectra_meta")(study, sm, "id")
    k <- match(scans$member_spectrum_id_, meta$spectrum_id_)
    expect_identical(meta$id[k], scans$member_native_id)
    expect_identical(meta$run_id[k], scans$member_run_id)
})

test_that("the dataset validates, and its results index is schema-valid", {
    for (study in c(FALSE, TRUE)) {
        p <- .results_dataset(study = study)
        v <- validateParquet(p)
        expect_identical(v[v$level == "MUST", "message"], character())
        idx <- jsonlite::fromJSON(file.path(p, "results", "column_map.json"),
                                  simplifyVector = FALSE)
        expect_identical(.mzs(".mzs_schema_check")(
            idx, "results_index.schema.json"), character())
    }
})

test_that("MsBackendParquet reads the derived spectra as a Spectra object", {
    skip_if_not_installed("MsBackendParquet")
    skip_if_not_installed("Spectra")
    p <- .results_dataset()
    s <- Spectra::Spectra(MsBackendParquet::backendInitialize(
        MsBackendParquet::MsBackendParquet(), path = p))
    expect_identical(length(s), 89L)
    expect_true(all(Spectra::msLevel(s) == 2L))
    expect_true(all(c("sn", "contributor_count", "x_mspurity_ra") %in%
                    Spectra::peaksVariables(s)))
    sp <- .read_run(p, "av_all")
    sd <- Spectra::spectraData(s, c("spectrumId", "rtime"))
    k <- match(sp$id, sd$spectrumId)
    expect_equal(sd$rtime[k], sp$time * 60)
})

test_that("an interrupted write leaves no dataset behind (B-064)", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    out <- tempfile("fail-")
    dir.create(out)
    old <- options(msPurity.parquet.fail_before_commit = TRUE)
    on.exit(options(old))
    expect_error(suppressMessages(createDatabase(
        .golden_pa(), .golden_xcms("msms_only_xcmsnexp.rds"), outDir = out,
        dbName = "x.parquet", format = "parquet")), "fail_before_commit")
    expect_identical(list.files(out, all.files = TRUE, no.. = TRUE),
                     character())
})

test_that("createDatabase() warns for SQLite, not for Parquet", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    out <- tempfile("warn-")
    dir.create(out)
    r <- .with_deprecations(suppressMessages(createDatabase(
        .golden_pa(), .golden_xcms("msms_only_xcmsnexp.rds"), outDir = out,
        dbName = "x.sqlite")))
    expect_identical(r$n, 1L)
    r <- .with_deprecations(suppressMessages(createDatabase(
        .golden_pa(), .golden_xcms("msms_only_xcmsnexp.rds"), outDir = out,
        dbName = "x.parquet", format = "parquet")))
    expect_identical(r$n, 0L)
})

test_that("createDatabase(format = 'parquet') refuses what it cannot map", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa()
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    out <- tempfile("refuse-")
    dir.create(out)
    expect_parquet_error(createDatabase(pa, x, outDir = out,
                                        grpPeaklist = data.frame(a = 1),
                                        format = "parquet"), "unsupported")
    expect_parquet_error(createDatabase(pa, x, xsa = list(), outDir = out,
                                        format = "parquet"), "unsupported")
    ## Never inside the study it reads (R-069).
    study <- .copy_dataset(.study_dataset())
    expect_parquet_error(createDatabase(pa, x, outDir = study,
                                        dbName = "inside", study = study,
                                        format = "parquet"), "semantic")
    ## An existing destination is kept unless overwrite = TRUE.
    p <- suppressMessages(createDatabase(pa, x, outDir = out, dbName = "d",
                                         format = "parquet"))
    expect_error(createDatabase(pa, x, outDir = out, dbName = "d",
                                format = "parquet"), "already exists")
    p2 <- suppressMessages(createDatabase(pa, x, outDir = out, dbName = "d",
                                          format = "parquet",
                                          overwrite = TRUE))
    expect_identical(p2, p)
    expect_identical(sort(list.files(out, all.files = TRUE, no.. = TRUE)),
                     "d")
})

test_that("two source files sharing a name are refused without a fileMap (R-096)", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    fi <- data.frame(fileid = 1:2, filename = c("QC.mzML", "QC.mzML"),
                     filepth = c("/a/QC.mzML", "/b/QC.mzML"), class = NA,
                     stringsAsFactors = FALSE)
    files <- .mzs(".mzs_files")
    expect_parquet_error(files(fi, NULL, "study", NULL, tempdir()),
                         "semantic")
    f <- files(fi, NULL, "study",
               c("/a/QC.mzML" = "QC_a", "/b/QC.mzML" = "QC_b"), tempdir())
    expect_identical(f$run_id, c("QC_a", "QC_b"))
    ## Minted run ids never come from the file name (R-095).
    f <- files(fi[1, ], NULL, "study", NULL, tempdir())
    expect_match(f$run_id, "^run-[0-9a-f]{12}$")
})
