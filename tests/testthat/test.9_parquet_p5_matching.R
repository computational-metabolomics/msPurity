context("Parquet: evidence, coverage, spectralMatching() and combineAnnotations()")

.qvq_args <- list(q_spectraTypes = "av_all", q_spectraFilter = TRUE,
                  q_pol = NA, l_spectraTypes = "av_all", l_pol = NA,
                  l_spectraFilter = TRUE, cores = 1)

.example <- function(...) system.file("extdata", "tests", ...,
                                      package = "msPurity")

# A results dataset converted from the example database, fresh each call.
.query_dataset <- function() {
    skip_if_no_arrow()
    p <- tempfile("q-")
    .without_identity_warning(convertSqliteToParquet(
        .example("db", "createDatabase_example.sqlite"), p))
    p
}

.sm <- function(q, l, ...) suppressMessages(do.call(
    spectralMatching, c(list(q_dbPth = q, l_dbPth = l, format = "parquet"),
                        utils::modifyList(.qvq_args, list(...)))))

test_that("Parquet matching scores exactly as SQLite matching", {
    p <- .query_dataset()
    db <- .example("db", "createDatabase_example.sqlite")
    a <- suppressMessages(do.call(spectralMatching,
                                  c(list(q_dbPth = db, l_dbPth = db), .qvq_args)))
    b <- .sm(p, p)
    expect_identical(nrow(b$matchedResults), nrow(a$matchedResults))
    si <- .read_table(p, "source_identifier")
    sp <- si[si$target_table == "spectra", ]
    pid <- function(id) as.numeric(sp$source_value[match(id, sp$target_key)])
    am <- a$matchedResults
    bm <- b$matchedResults
    k <- match(paste(am$qpid, am$lpid), paste(pid(bm$qpid), pid(bm$lpid)))
    expect_false(anyNA(k))
    for (s in c("dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent"))
        expect_identical(as.numeric(bm[[s]][k]), as.numeric(am[[s]]), info = s)
})

test_that("updateDb = FALSE writes nothing; TRUE appends evidence (R-070)", {
    p <- .query_dataset()
    before <- tools::md5sum(list.files(p, recursive = TRUE, full.names = TRUE))
    g0 <- .manifest(p)$generation
    .sm(p, p)
    expect_identical(.manifest(p)$generation, g0)
    ## An unknown key survives the update (G-011).
    m <- MsBackendParquet::readManifest(p)
    m$x_vendor <- list("keep", list(a = 1L))
    .write_raw_manifest(p, m)
    .sm(p, p, updateDb = TRUE)
    m <- .manifest(p)
    expect_identical(m$generation, g0 + 1L)
    expect_identical(m$x_vendor, list("keep", list(a = 1L)))
    expect_identical(m$provenance[[2]]$agent$`function`, "spectralMatching")
    ## Existing parts are byte-for-byte unchanged.
    old <- names(before)[!endsWith(names(before), "mzStack.json")]
    expect_identical(tools::md5sum(old), before[old])
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("coverage counts equal the evidence rows of each query (B-044)", {
    p <- .query_dataset()
    .sm(p, p, updateDb = TRUE)
    ev <- .read_table(p, "evidence")
    cv <- .read_table(p, "coverage")
    n <- table(factor(ev$query_spectrum_id_, levels = cv$query_spectrum_id_))
    expect_identical(as.integer(n), cv$hits_retained)
    expect_true(all(cv$candidates_considered >= cv$hits_retained))
    es <- .read_table(p, "evidence_score")
    expect_identical(nrow(es), 6L * nrow(ev))
    expect_true(all(es$score_kind %in% c("raw", "count")))
    expect_true(all(ev$identification_method == "MS:1001031"))
})

test_that("not attempted, nothing matched and cut lists are distinguishable (R-073)", {
    p <- .query_dataset()
    ## A window too narrow to match anything but itself, and a top-1 cut.
    .sm(p, p, updateDb = TRUE, topn = 1, q_ppmPrec = 1000, l_ppmPrec = 1000)
    cv <- .read_table(p, "coverage")
    m <- .manifest(p)
    expect_identical(m$results$tables$evidence$selection$policy, "top_n")
    expect_identical(m$results$tables$evidence$selection$ties, "all")
    cut <- cv[cv$truncated, ]
    expect_true(nrow(cut) > 0)
    expect_false(anyNA(cut$boundary_score))
    expect_true(all(cut$hits_retained >= 1L))
    expect_true(all(is.na(cv$boundary_score[!cv$truncated])))
    ## Queries the activity did not attempt have no coverage row; those it
    ## did all have one, matched or not.
    sp <- .read_run(p, "av_all")
    expect_true(all(cv$query_spectrum_id_ %in% sp$spectrum_id_))
})

test_that("a search without candidates returns and records no matches", {
    p <- .query_dataset()
    lib <- tempfile("lib-")
    suppressMessages(convertLibraryToParquet(
        .example("library", "mini_mona.msp"), lib, format = "msp",
        dialect = "mona"))
    ## No query has a candidate among the negative-mode library spectra.
    none <- .sm(p, lib, l_spectraTypes = NA, l_pol = "negative",
                updateDb = TRUE)
    expect_null(none$matchedResults)
    expect_true(all(.read_table(p, "coverage")$candidates_considered == 0))
    ## Matches written later leave a valid dataset.
    .sm(p, lib, l_spectraTypes = NA, l_pol = "positive", updateDb = TRUE)
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("unknown query polarity is matched unless a polarity is required", {
    p <- .query_dataset()
    expected <- .sm(p, p)$matchedResults
    read_meta <- .mzs(".mzs_read_spectra_meta")
    local_mocked_bindings(.mzs_read_spectra_meta = function(...) {
        meta <- read_meta(...)
        if (!is.null(meta$scan_polarity))
            meta$scan_polarity <- rep(NA_integer_, nrow(meta))
        meta
    }, .package = "msPurity")
    expect_identical(.sm(p, p)$matchedResults, expected)
    expect_null(.sm(p, p, q_pol = "positive")$matchedResults)
})

test_that("an empty selection from a results library returns no matches", {
    p <- .query_dataset()
    expect_null(.sm(p, p, l_xcmsGroups = -1)$matchedResults)
})

test_that("a dataset without a library role is refused clearly", {
    p <- .query_dataset()
    l <- .copy_dataset(p)
    m <- MsBackendParquet::readManifest(l)
    m$role <- NULL
    .write_raw_manifest(l, m)
    expect_parquet_error(.sm(p, l), "unsupported")
})

test_that("ranks follow mzTab-M: ties share a rank (R-051)", {
    rk <- .mzs(".mzs_rank_matches")
    m <- data.frame(qkey = c("a", "a", "a", "a", "b"),
                    dpc = c(0.9, 0.5, 0.9, 0.2, 0.1))
    att <- data.frame(qkey = c("a", "b", "c"), query_source = "self",
                      query_run_id = "r", query_spectrum_id_ = 1:3,
                      query_native_id = NA)
    r <- rk(m, att, topn = 1)
    expect_identical(r$matches$rank, c(1L, 3L, 1L, 4L, 1L))
    expect_identical(r$matches$keep, c(TRUE, FALSE, TRUE, FALSE, TRUE))
    expect_identical(r$coverage$hits_retained, c(2L, 1L, 0L))
    expect_identical(r$coverage$truncated, c(TRUE, FALSE, FALSE))
    expect_identical(r$coverage$boundary_score, c(0.9, NA, NA))
})

test_that("evidence keeps one selection per table", {
    p <- .query_dataset()
    .sm(p, p, updateDb = TRUE)
    expect_parquet_error(.sm(p, p, updateDb = TRUE, topn = 5), "unsupported")
    .sm(p, p, updateDb = TRUE)
    expect_identical(length(.manifest(p)$provenance), 3L)
})

test_that("copyDb = TRUE writes to a fork with its own uid (D-011)", {
    p <- .query_dataset()
    out <- tempfile("fork-")
    .sm(p, p, updateDb = TRUE, copyDb = TRUE, outPth = out)
    expect_null(.manifest(p)$results$tables$evidence)
    m <- .manifest(out)
    expect_false(identical(m$uid, .manifest(p)$uid))
    expect_identical(m$provenance[[2]]$parameters$forked_from,
                     .manifest(p)$uid)
    v <- validateParquet(out)
    expect_identical(v$message[v$level == "MUST"], character())
})

test_that("what the Parquet route cannot apply is refused (B-005)", {
    p <- .query_dataset()
    expect_parquet_error(.sm(p, p, q_pids = 1), "unsupported")
    expect_parquet_error(.sm(p, p, q_raThres = 5), "unsupported")
    expect_parquet_error(suppressMessages(spectralMatching(
        p, format = "parquet")), "unsupported")
    expect_error(suppressMessages(spectralMatching(
        .example("db", "createDatabase_example.sqlite"), topn = 3)),
        "only available")
})

test_that("spectralMatching() warns only when it writes SQLite", {
    skip_if_no_arrow()
    db <- .example("db", "createDatabase_example.sqlite")
    r <- .with_deprecations(suppressMessages(do.call(
        spectralMatching, c(list(q_dbPth = db, l_dbPth = db), .qvq_args))))
    expect_identical(r$n, 0L)
    r <- .with_deprecations(suppressMessages(do.call(
        spectralMatching, c(list(q_dbPth = db, l_dbPth = db, updateDb = TRUE,
                                 copyDb = TRUE, outPth = tempfile(
                                     fileext = ".sqlite")), .qvq_args))))
    expect_identical(r$n, 1L)
})

test_that("a database's matches become evidence, compounds and coverage", {
    skip_if_no_arrow()
    p <- tempfile("sm-")
    .without_identity_warning(convertSqliteToParquet(
        .example("sm", "spectralMatching_result.sqlite"), p))
    ev <- .read_table(p, "evidence")
    expect_identical(nrow(ev), 2L)
    expect_identical(sort(ev$reference_native_id),
                     c("CCMSLIB00000577898", "CE000616"))
    expect_true(all(startsWith(ev$reference_source, "library_")))
    m <- .manifest(p)
    lib <- Filter(function(s) startsWith(s$key, "library_"), m$sources)
    expect_true(all(vapply(lib, function(s) s$resolution, "") == "external"))
    cmp <- .read_table(p, "compound")
    expect_setequal(cmp$inchikey, c("ONIBWKKTOPOVIA-UHFFFAOYSA-N",
                                    "AGPKZVBTJJNPAG-UHFFFAOYSA-N"))
    expect_true(all(ev$compound_id_ %in% cmp$compound_id_))
    expect_true(all(ev$compound_match_level == "exact_structure"))
    cv <- .read_table(p, "coverage")
    expect_identical(nrow(cv), 2L)
    ll <- .read_table(p, "loss_ledger")
    expect_true("spectral matching scope" %in% ll$source_construct)
    v <- validateParquet(p)
    expect_identical(v$message[v$level == "MUST"], character())
})

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

test_that("anonymous Spectra references survive reordering and queued peak processing", {
    p <- .query_dataset()
    lib <- Spectra::Spectra(S4Vectors::DataFrame(
        msLevel = c(2L, 2L), precursorMz = c(123.06, 123.06),
        polarity = c(1L, 1L), name = c("first", "second"),
        spectrumId = c(NA_character_, "native-2"), accession = c("", NA_character_),
        mz = IRanges::NumericList(c(80, 100), c(91, 100)),
        intensity = IRanges::NumericList(c(15, 50), c(20, 25))))
    search <- function(x, ...) .sm(p, x, l_spectraTypes = NA, usePrecursors = FALSE,
                                   ...)$matchedResults
    mr <- search(lib)
    old <- options(msPurity.parquet.batch_spectra = 1L)
    on.exit(options(old))
    expect_identical(search(lib), mr)
    ids <- unique(mr$library_accession)
    expect_length(ids, 2L)
    hash <- ids[startsWith(ids, "sha256:")]
    expect_match(hash, "^sha256:[[:xdigit:]]{64}$")
    expect_true("native-2" %in% ids)
    expect_identical(unique(search(lib[2:1], l_accessions = hash)$library_accession), hash)
    scaled <- Spectra::addProcessing(lib[1], function(x, ...) {
        x[, "intensity"] <- x[, "intensity"] * 2
        x
    })
    expect_false(hash %in% search(scaled)$library_accession)
})

test_that("precursor selection needs only the matching rows of an unordered library", {
    # Expanding the lazy peaks to the converter's text representation would
    # exceed R's single-string limit. Reject any out-of-window peak read.
    p <- .query_dataset()
    n <- 32769L
    lib <- Spectra::Spectra(S4Vectors::DataFrame(
        msLevel = c(2L, 2L), polarity = c(1L, 1L),
        precursorMz = c(1e6, 123.0553), accession = c("far", "near"),
        mz = IRanges::NumericList(pi, c(80, 100)),
        intensity = IRanges::NumericList(exp(1), c(15, 50))))
    lib <- lib[c(rep(1L, n %/% 2L), 2L, rep(1L, n %/% 2L))]
    guard <- new.env()
    guard$searching <- FALSE
    lib <- Spectra::addProcessing(lib, function(x, ...) {
        if (x[1L, "mz"] != pi) return(x)
        if (guard$searching) stop("An out-of-window spectrum was read")
        x[rep(1L, 4097L), , drop = FALSE]
    })
    far <- Spectra::peaksData(lib[1L])[[1L]]
    entry <- paste(format(far[, "mz"], digits = 17), format(far[, "intensity"], digits = 17),
                   collapse = " ")
    expect_gt(as.double(n - 1L) * nchar(entry, type = "bytes"), .Machine$integer.max)
    guard$searching <- TRUE
    mr <- .sm(p, lib, l_spectraTypes = NA)$matchedResults
    expect_gt(nrow(mr), 0L)
    expect_identical(unique(mr$library_accession), "near")
    expect_identical(unique(mr$lpid), as.numeric(n %/% 2L + 1L))
})
