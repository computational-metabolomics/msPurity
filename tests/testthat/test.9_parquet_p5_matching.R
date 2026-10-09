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
    convertSqliteToParquet(.example("db", "createDatabase_example.sqlite"), p)
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
    m <- .mzs(".mzs_read_manifest")(p)
    m$x_vendor <- list("keep", list(a = 1L))
    .mzs(".mzs_write_manifest")(p, m)
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
    convertSqliteToParquet(.example("sm", "spectralMatching_result.sqlite"), p)
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
    convertSqliteToParquet(.example("sm", "spectralMatching_result.sqlite"), p)
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
