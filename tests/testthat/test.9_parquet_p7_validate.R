context("Parquet: validation (mzStack-1 §11, mzStack-4 §11)")

.valid_matched <- local({
    path <- NULL
    function() {
        skip_if_no_arrow()
        if (is.null(path) || !dir.exists(path)) {
            ex <- function(...) system.file("extdata", "tests", ...,
                                            package = "msPurity")
            lib <- tempfile("lib-")
            convertLibraryToParquet(ex("library", "mini_library.sqlite"), lib)
            q <- tempfile("q-")
            .without_identity_warning(convertSqliteToParquet(
                ex("db", "createDatabase_example.sqlite"), q))
            suppressMessages(spectralMatching(q, lib, q_spectraTypes = "av_all",
                                              format = "parquet",
                                              updateDb = TRUE, topn = 1))
            path <<- q
        }
        path
    }
})

# A copy of a valid dataset, mutated by `fn(path)`, and its findings.
.mutate <- function(fn, base = .valid_matched()) {
    p <- .copy_dataset(base)
    ## Sources recorded relative to the original location still resolve.
    m <- MsBackendParquet::readManifest(p)
    for (i in seq_along(m$sources))
        if (!is.null(m$sources[[i]]$path) &&
            !.mzs(".mzs_is_absolute")(m$sources[[i]]$path))
            m$sources[[i]]$path <- normalizePath(file.path(
                base, m$sources[[i]]$path))
    .write_raw_manifest(p, m)
    fn(p)
    validateParquet(p)
}

.manifest_edit <- function(p, fn) {
    m <- MsBackendParquet::readManifest(p)
    .write_raw_manifest(p, fn(m))
}

.rewrite_table <- function(p, name, fn) {
    m <- MsBackendParquet::readManifest(p)
    root <- file.path(p, m$results$tables[[name]]$path)
    for (f in list.files(root, pattern = "\\.parquet$", recursive = TRUE,
                         full.names = TRUE)) {
        t <- arrow::read_parquet(f, as_data_frame = FALSE, mmap = FALSE)
        out <- fn(as.data.frame(t), t)
        unlink(f)
        arrow::write_parquet(out, f)
    }
}

.has <- function(v, req) any(v$level == "MUST" & v$requirement == req)

test_that("every dataset kind msPurity writes validates", {
    v <- validateParquet(.valid_matched())
    expect_identical(v$message[v$level == "MUST"], character())
    expect_true(all(v$requirement[v$level == "INFO"] == "D-049"))
})

test_that("structural defects are found (R-037, R-038, R-061, D-005)", {
    expect_true(.has(.mutate(function(p) .manifest_edit(p, function(m) {
        m$results$tables$feature$rows <- 1L
        m
    })), "R-037"))
    expect_true(.has(.mutate(function(p) .manifest_edit(p, function(m) {
        m$results$tables$feature$path <- "results/feature/rev-9"
        m
    })), "R-037"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "feature",
        function(df, t) {
            df <- df[rev(seq_len(nrow(df))), ]
            arrow::arrow_table(df, schema = t$schema)
        })), "R-038"))
    expect_true(.has(.mutate(function(p) .manifest_edit(p, function(m) {
        m$generation <- 0L
        m
    })), "D-005"))
    expect_true(.has(.mutate(function(p) {
        m <- MsBackendParquet::readManifest(p)
        cm <- file.path(p, m$results$column_mapping_path)
        idx <- jsonlite::fromJSON(cm, simplifyVector = FALSE)
        idx$files <- Filter(function(f) f$name != "feature", idx$files)
        writeLines(jsonlite::toJSON(idx, auto_unbox = TRUE), cm)
    }), "R-061"))
})

test_that("reference and integrity defects are found (R-018, R-020, R-024, R-099)", {
    expect_true(.has(.mutate(function(p) .manifest_edit(p, function(m) {
        m$sources[[2]] <- m$sources[[1]]
        m
    })), "R-024"))
    expect_true(.has(.mutate(function(p) .manifest_edit(p, function(m) {
        m$sources[[1]]$key <- "self"
        m
    })), "R-020"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "abundance",
        function(df, t) {
            df$feature_id_[1] <- 99999
            arrow::arrow_table(df, schema = t$schema)
        })), "R-099"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "evidence",
        function(df, t) {
            df$reference_source <- factor("nowhere")
            arrow::arrow_table(df, schema = t$schema)
        })), "R-018"))
})

test_that("derived-spectrum and ranking defects are found (R-047, R-051)", {
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "merge_member",
        function(df, t) arrow::arrow_table(df[df$merged_spectrum_id_ != 1, ],
                                           schema = t$schema))), "R-047"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "evidence",
        function(df, t) {
            df$rank <- df$rank + 1L
            arrow::arrow_table(df, schema = t$schema)
        })), "R-051"))
})

test_that("coverage and ledger defects are found (B-044, R-076, R-083, R-084)", {
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "coverage",
        function(df, t) {
            df$hits_retained <- df$hits_retained + 1L
            arrow::arrow_table(df, schema = t$schema)
        })), "B-044"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "coverage",
        function(df, t) {
            df$truncated <- TRUE
            df$boundary_score <- NA_real_
            arrow::arrow_table(df, schema = t$schema)
        })), "R-076"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "loss_ledger",
        function(df, t) {
            df$disposition <- factor("lost")
            arrow::arrow_table(df, schema = t$schema)
        })), "R-083"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "loss_ledger",
        function(df, t) arrow::arrow_table(
            df[!(df$target_ref %in% "compound.uri"), ], schema = t$schema))),
        "R-084"))
})

test_that("naming and InChIKey defects are found (G-016, R-049)", {
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "compound",
        function(df, t) {
            df$inchikey[1] <- "NOT-A-KEY"
            arrow::arrow_table(df, schema = t$schema)
        })), "R-049"))
    expect_true(.has(.mutate(function(p) .rewrite_table(p, "assay",
        function(df, t) {
            df$x_mspurity_sample_id <- 1L
            arrow::arrow_table(df)
        })), "G-016"))
})

test_that("an unknown CV accession is a MUST violation (R-063)", {
    expect_true(.has(.mutate(function(p) {
        m <- MsBackendParquet::readManifest(p)
        cm <- file.path(p, m$results$column_mapping_path)
        idx <- jsonlite::fromJSON(cm, simplifyVector = FALSE)
        idx$files[[1]]$column_mapping[[1]]$accession <- "MS:9999999"
        writeLines(jsonlite::toJSON(idx, auto_unbox = TRUE), cm)
    }), "R-063"))
})

test_that("an unknown results major still reads as spectra (R-035)", {
    p <- .copy_dataset(.valid_matched())
    .manifest_edit(p, function(m) {
        m$results$version <- "2.0.0"
        m
    })
    m <- MsBackendParquet::readManifest(p)
    expect_true(sum(.mzs(".mzs_runs_frame")(m)$n_spectra) > 0)
    expect_parquet_error(.mzs(".mzs_read_table")(p, m, "feature"),
                         "unsupported")
    expect_parquet_error(suppressMessages(spectralMatching(
        p, p, format = "parquet")), "unsupported")
    expect_true(.has(validateParquet(p), "R-034"))
    s <- Spectra::Spectra(MsBackendParquet::backendInitialize(
        MsBackendParquet::MsBackendParquet(), path = p))
    expect_identical(length(s), sum(.mzs(".mzs_runs_frame")(m)$n_spectra))
})

test_that("a dataset moved after it is written stays valid and readable", {
    p <- .copy_dataset(.results_dataset())
    runs <- .mzs(".mzs_runs_frame")(MsBackendParquet::readManifest(p))
    expect_identical(runs$path, paste0("spectra/run_id=", runs$run_id))
    moved <- paste0(p, "-moved")
    expect_true(file.rename(p, moved))
    v <- validateParquet(moved)
    expect_identical(v$message[v$level == "MUST"], character())
    s <- Spectra::Spectra(MsBackendParquet::backendInitialize(
        MsBackendParquet::MsBackendParquet(), path = moved))
    expect_identical(length(s), sum(runs$n_spectra))
    expect_identical(names(.mzs(".mzs_read_peaks")(moved, c(2, 1))),
                     c("1", "2"))
})

test_that("unresolvable sources are SHOULD findings, never failures (R-100)", {
    p <- .copy_dataset(.results_dataset(study = TRUE))
    missing <- tempfile("missing-study-")
    expect_false(file.exists(missing))
    .manifest_edit(p, function(m) {
        m$sources[[1]]$path <- missing
        m
    })
    v <- validateParquet(p)
    expect_false(.has(v, "R-100"))
    expect_true(any(v$requirement == "R-100" & v$level == "SHOULD"))
})
