test_that("combineAnnotations() ranks the same on Parquet as on SQLite", {
    skip_if_no_arrow()
    p <- tempfile("ca-")
    .without_identity_warning(convertSqliteToParquet(
        .example("sm", "spectralMatching_result.sqlite"), p))
    args <- list(
        compoundDbPth = .example("db", "metab_compound_subset.sqlite"),
        metfrag_resultPth = .example("external_annotations", "metfrag.tsv"),
        sirius_csi_resultPth = .example("external_annotations",
                                        "sirus_csifingerid.tsv"),
        probmetab_resultPth = .example("external_annotations", "probmetab.tsv"),
        ms1_lookup_resultPth = .example("external_annotations", "beams.tsv"),
        weights = list(sm = 0.3, metfrag = 0.2, sirius_csifingerid = 0.2,
                       probmetab = 0, ms1_lookup = 0.05, biosim = 0.25),
        ms1_lookup_dbSource = "hmdb", ms1_lookup_checkAdducts = FALSE,
        ms1_lookup_keepAdducts = c("[M+H]+", "[M-H]-"))
    db <- tempfile(fileext = ".sqlite")
    file.copy(.example("sm", "spectralMatching_result.sqlite"), db)
    r <- suppressWarnings(.with_deprecations(suppressMessages(
        do.call(combineAnnotations, c(list(db), args)))))
    expect_identical(r$n, 1L)
    a <- r$value
    r <- suppressWarnings(.with_deprecations(suppressMessages(
        do.call(combineAnnotations, c(list(p, format = "parquet"), args)))))
    expect_identical(r$n, 0L)
    b <- r$value
    expect_identical(b$rank, a$rank)
    expect_identical(b$wscore, a$wscore)
    expect_identical(b$inchikey, a$inchikey)
    ## Recorded without changing any existing row: compound is superseded.
    m <- .manifest(p)
    expect_identical(m$results$tables$compound$revision, 2L)
    expect_identical(m$results$tables$compound$supersedes, 1L)
    ca <- .read_table(p, "x_mspurity_combined_annotation")
    expect_identical(nrow(ca), nrow(a))
    fa <- .read_table(p, "x_mspurity_feature_annotation")
    expect_true(all(fa$score_kind[fa$tool == "sirius_csifingerid"] ==
                    "rescaled"))
    inputs <- m$provenance[[2]]$inputs
    expect_true(all(vapply(inputs, function(i) nchar(i$sha256), 1L) == 64L))
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("Spectra libraries match converted libraries across backends and batches", {
    skip_if_not_installed("MsBackendMsp")
    p <- .query_dataset()
    msp <- MsBackendMsp::MsBackendMsp()
    lib <- Spectra::Spectra(.example("library", "mini_mona.msp"), source = msp,
        mapping = Spectra::spectraVariableMapping(msp, "mona"))
    path <- tempfile()
    convertLibraryToParquet(lib, path, version = "fixture")
    search <- function(x, ...) .sm(p, x, l_spectraTypes = NA, ...)
    expected <- search(path)$matchedResults
    expect_gt(nrow(expected), 0L)
    disk <- Spectra::setBackend(lib, MsBackendParquet::MsBackendParquet(), path = tempfile())
    old <- options(msPurity.parquet.batch_spectra = 1L)
    on.exit(options(old))
    fields <- c("qpid", "library_accession", "library_entry_name", "library_compound_name",
                "library_precursor_type", "inchikey", "library_precursor_mz", "rank",
                "dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent")
    compare <- function(mr) {
        out <- mr[order(mr$qpid, mr$library_accession), fields]
        rownames(out) <- NULL
        out
    }
    checksum <- paste(rep("ab", 32), collapse = "")
    for (x in list(lib, disk)) {
        actual <- search(x, libraryVersion = "fixture", libraryChecksum = checksum,
                         updateDb = TRUE)$matchedResults
        expect_identical(compare(actual), compare(expected))
    }
    src <- Filter(function(s) identical(s$version, "fixture"), .manifest(p)$sources)
    expect_length(src, 1L)
    expect_identical(src[[1]]$digest, paste0("sha256:", checksum))
    expect_identical(src[[1]]$resolution, "external")
    ev <- .read_table(p, "evidence")
    ev <- ev[ev$reference_source == src[[1]]$key, ]
    expect_gt(nrow(ev), 0L)
    expect_true(all(is.na(ev$reference_spectrum_id_)))
    expect_setequal(ev$reference_native_id, expected$library_accession)
    scores <- .read_table(p, "evidence_score")
    counts <- table(factor(scores$evidence_id_, levels = ev$evidence_id_))
    expect_identical(as.integer(counts), rep(6L, nrow(ev)))
    # A later library can introduce new compounds as well as new evidence.
    new <- lib[which(lib$accession %in% expected$library_accession)[1L]]
    new$inchikey <- "ABCDEFGHIJKLMN-ABCDEFGHIJ-A"
    search(new, updateDb = TRUE)
    updated <- .read_table(p, "evidence")
    latest <- updated[updated$activity == max(updated$activity), ]
    cmp <- .read_table(p, "compound")
    expect_gt(nrow(latest), 0L)
    expect_identical(unique(cmp$inchikey[match(latest$compound_id_, cmp$compound_id_)]),
                     new$inchikey)
    expect_parquet_error(.mzs(".mzs_commit_update")(p, expected_generation = -1L),
                         "conflict")
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character(), info = paste(v$message, collapse = "\n"))
    expect_null(search(lib[integer()])$matchedResults)
    expect_null(search(lib, l_accessions = "absent")$matchedResults)
    expect_null(search(lib, l_pol = "negative")$matchedResults)
    expect_parquet_error(search(lib, libraryChecksum = "bad"), "format")
})
