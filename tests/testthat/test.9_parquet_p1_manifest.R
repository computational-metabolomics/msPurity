context("Parquet: manifest, sources, references and conditions")

test_that("the nine error classes are distinct conditions (B-011)", {
    skip_if_no_arrow()
    abort <- .mzs(".mzs_abort")
    for (cl in c("format", "archive", "reference", "study", "stale",
                 "unsupported", "semantic", "capability", "resource")) {
        e <- tryCatch(abort(cl, "boom"), error = function(e) e)
        expect_s3_class(e, paste0("mzstack_", cl))
        expect_s3_class(e, "mzstack_error")
        expect_s3_class(e, paste0("msPurity_parquet_", cl))
        expect_s3_class(e, "msPurity_parquet_error")
        expect_identical(conditionMessage(e), "boom")
    }
    expect_error(abort("nonsense", "x"))
    e <- tryCatch(abort("reference", "gone", data = list(source = "s")),
                  error = function(e) e)
    expect_identical(e$source, "s")
})

test_that("the manifest round-trips every double exactly (D-070)", {
    skip_if_no_arrow()
    m <- MsBackendParquet::newManifest()
    tricky <- c(0.010000000000000002, 0.1 + 0.2, 1 / 3, 5e-324, 1e300,
                -0.0, 123456789.123456789)
    m$parameters <- list(tricky = as.list(tricky), one = 0.1)
    td <- .fresh_path()
    MsBackendParquet::writeManifest(td, m)
    back <- MsBackendParquet::readManifest(td)
    expect_identical(vapply(back$parameters$tricky, as.numeric, 0), tricky)
    expect_identical(back$parameters$one, 0.1)
    expect_identical(.mzs(".mzs_schema_check")(back, "manifest.schema.json"),
                     character())
})

test_that("unknown manifest keys survive a commit (G-011, D-065)", {
    skip_if_no_arrow()
    td <- .fresh_path()
    m <- MsBackendParquet::newManifest()
    m$x_vendor <- list(list(a = 1L), list(b = list("only")), "c")
    m$empty_obj <- .mzs(".mzs_object")()
    MsBackendParquet::writeManifest(td, m)
    before <- jsonlite::fromJSON(file.path(td, "mzStack.json"),
                                 simplifyVector = FALSE)
    m2 <- MsBackendParquet::readManifest(td)
    out <- .mzs(".mzs_commit_manifest")(td, m2, m2$generation)
    after <- jsonlite::fromJSON(file.path(td, "mzStack.json"),
                                simplifyVector = FALSE)
    expect_identical(after$x_vendor, before$x_vendor)
    expect_identical(after$empty_obj, before$empty_obj)
    expect_identical(after$generation, before$generation + 1L)
    expect_identical(out$generation, after$generation)
})

test_that("a manifest commit refuses a concurrent change", {
    skip_if_no_arrow()
    td <- .fresh_path()
    MsBackendParquet::writeManifest(td, MsBackendParquet::newManifest())
    m <- MsBackendParquet::readManifest(td)
    .mzs(".mzs_commit_manifest")(td, m, m$generation)
    e <- expect_error(.mzs(".mzs_commit_manifest")(td, m, m$generation),
                      class = "msPurity_parquet_conflict")
    expect_s3_class(e, "mzstack_conflict")
    expect_identical(MsBackendParquet::readManifest(td)$generation,
                     m$generation + 1L)
})

test_that("a manifest commit refuses while another writer holds the lock", {
    skip_if_no_arrow()
    td <- .fresh_path()
    MsBackendParquet::writeManifest(td, MsBackendParquet::newManifest())
    m <- MsBackendParquet::readManifest(td)
    lock <- file.path(td, "mzStack.json.lock")
    dir.create(lock)
    expect_error(.mzs(".mzs_commit_manifest")(td, m, m$generation),
                 class = "msPurity_parquet_conflict")
    expect_identical(MsBackendParquet::readManifest(td)$generation,
                     m$generation)
    unlink(lock, recursive = TRUE)
    expect_identical(
        .mzs(".mzs_commit_manifest")(td, m, m$generation)$generation,
        m$generation + 1L)
})

test_that("non-mzStack manifests are Format errors (G-008, G-009)", {
    skip_if_no_arrow()
    td <- .fresh_path()
    dir.create(td)
    expect_parquet_error(MsBackendParquet::readManifest(td), "format")
    writeLines('{"format": "other", "version": "0.1.0"}',
               file.path(td, "mzStack.json"))
    expect_parquet_error(MsBackendParquet::readManifest(td), "format")
    writeLines('{"format": "mzStack", "version": "9.0.0", "generation": 1,
               "created": "2026-01-01T00:00:00Z", "runs": []}',
               file.path(td, "mzStack.json"))
    expect_parquet_error(MsBackendParquet::readManifest(td), "format")
})

test_that("the schema checker rejects known-bad documents", {
    skip_if_no_arrow()
    check <- .mzs(".mzs_schema_check")
    good <- MsBackendParquet::newManifest()
    expect_identical(check(good, "manifest.schema.json"), character())
    bad <- good
    bad$generation <- 0L
    expect_true(length(check(bad, "manifest.schema.json")) > 0L)
    bad <- good
    bad$format <- "nope"
    expect_true(length(check(bad, "manifest.schema.json")) > 0L)
    bad <- good
    bad$runs <- list(list(run_id = "bad id!", kind = "native"))
    expect_true(length(check(bad, "manifest.schema.json")) > 0L)
    bad <- good
    bad$role <- "results"
    expect_true(length(check(bad, "manifest.schema.json")) > 0L)
    bad <- good
    bad$sources <- list(list(key = "self", resolution = "external"))
    bad$uid <- "0123456789abcdef"
    expect_true(length(check(bad, "manifest.schema.json")) > 0L)
    expect_identical(check(parquetCapabilities(),
                           "capabilities.schema.json"), character())
})

test_that("source keys follow R-024 and never name 'self' (R-020)", {
    skip_if_no_arrow()
    ck <- .mzs(".mzs_check_key")
    expect_silent(ck("study_A"))
    expect_parquet_error(ck("self"), "format")
    expect_parquet_error(ck("no spaces"), "format")
    ext <- .mzs(".mzs_external_source")("run1", location = "file:///x.mzML")
    expect_null(ext$uid)
    expect_true("uid" %in% names(ext))
})

test_that("sources resolve lazily: Reference, Capability, Stale (R-022, B-015)", {
    skip_if_no_arrow()
    study <- .copy_dataset(.study_dataset())
    ds <- .fresh_path("results-")
    dir.create(ds)
    m <- MsBackendParquet::newManifest()
    m$uid <- .mzs(".mzs_uid")()
    m$sources <- list(
        .mzs(".mzs_dataset_source")("study_A", study, ds),
        .mzs(".mzs_external_source")("ext", location = "file:///QC01.mzML"))
    expect_identical(.mzs(".mzs_schema_check")(m, "manifest.schema.json"),
                     character())
    MsBackendParquet::writeManifest(ds, m)
    m <- MsBackendParquet::readManifest(ds)
    open <- .mzs(".mzs_open_source")

    ok <- open(ds, m, "study_A")
    expect_identical(normalizePath(ok$path), normalizePath(study))
    expect_parquet_error(open(ds, m, "ext"), "capability")

    ## A fingerprint edit is Stale; an unrelated generation bump is not.
    sm <- MsBackendParquet::readManifest(study)
    sm <- MsBackendParquet::writeManifest(study, sm)
    expect_silent(open(ds, m, "study_A"))
    sm$runs[[1]]$ingested_at <- sm$generation
    MsBackendParquet::writeManifest(study, sm, bump = FALSE)
    expect_parquet_error(open(ds, m, "study_A"), "stale")

    ## Moved away: Reference, unless the caller says where it went.
    moved <- paste0(study, "-moved")
    file.rename(study, moved)
    expect_parquet_error(open(ds, m, "study_A"), "reference")
    st <- .mzs(".mzs_reference_status")(ds, m, c("study_A", "ext"),
                                        c(NA, NA))
    expect_identical(st, c("unresolved_missing", "unresolved_external"))
})

test_that("validateParquet() reports a directory without a manifest", {
    skip_if_no_arrow()
    td <- .fresh_path()
    dir.create(td)
    v <- validateParquet(td)
    expect_identical(v$requirement, "D-001")
    expect_identical(v$level, "MUST")
})
