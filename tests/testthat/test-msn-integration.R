test_that("MSn association appends options after existing arguments", {
    existing <- c("pa", "xcmsObj", "ppm", "plim", "intense", "convert2RawRT",
        "useGroup", "createDb", "outDir", "dbName", "grpPeaklist", "use_group",
        "out_dir", "create_db", "grp_peaklist", "db_name", "xset")
    generic <- methods::getGeneric("frag4feature")
    method <- methods::selectMethod("frag4feature", "purityA")
    expected <- c(existing, "includeMSn", "parentReference",
        "purityPpm", "purityTolerance")
    expect_identical(names(formals(generic)), expected)
    expect_identical(names(formals(method)), expected)
    expect_identical(formals(method)$includeMSn, FALSE)
    expect_identical(formals(method)$parentReference, "auto")
    expect_false("linkMSnSpectra" %in% getNamespaceExports("msPurity"))
    pa <- methods::new("purityA")
    expect_error(frag4feature(pa, NULL, includeMSn = NA), "includeMSn")
    expect_error(frag4feature(pa, NULL, parentReference = "invalid"),
        "arg.*should be one of")
    expect_error(frag4feature(pa, NULL, purityPpm = -1), "purityPpm")
    expect_error(frag4feature(pa, NULL, purityTolerance = NA),
        "purityTolerance")
    expect_error(frag4feature(pa, NULL, includeMSn = TRUE, createDb = TRUE),
        "deprecated createDb route")
})

test_that("stored ancestry follows acquisition metadata rather than flags", {
    source <- .ancestry_spectra()
    metadata <- msPurity:::.pa_msn_spectra_metadata(source)
    preparation <- msPurity:::.pa_msn_ancestry(metadata, "sourceSpectrumId")
    preparation$source_metadata <- metadata
    pa <- methods::new("purityA")
    pa@spectra <- source
    pa@params$msnLineage <- preparation
    expect_true(msPurity:::.pa_msn_validate_preparation(pa))
    pa@spectra$purity_pass_flag <- c(TRUE, NA, NA)
    expect_true(msPurity:::.pa_msn_validate_preparation(pa))
    reordered <- pa
    reordered@spectra <- pa@spectra[c(3L, 1L, 2L)]
    expect_true(msPurity:::.pa_msn_validate_preparation(reordered))
    changed <- pa
    changed@spectra$precursorMz[2L] <- 151
    expect_error(msPurity:::.pa_msn_validate_preparation(changed),
        "no longer matches acquisition metadata")
    changed <- pa
    changed@spectra$precScanNum[3L] <- 101L
    expect_error(msPurity:::.pa_msn_validate_preparation(changed),
        "no longer matches acquisition metadata")
    changed <- pa
    changed@spectra <- changed@spectra[1:2]
    expect_error(msPurity:::.pa_msn_validate_preparation(changed),
        "no longer matches the source spectra")
})

.msn_association_fixture <- local({
    cache <- list()
    function(compound) {
        if (is.null(cache[[compound]])) {
            prepared <- .msn_standard_feature(compound, "MSn")
            stopifnot(length(prepared$selected) == 1L,
                length(prepared$feature) == 1L)
            prepared$pa <- purityA(.msn_standard_file(compound, "MSn"),
                cores = 1L, plotP = FALSE)
            cache[[compound]] <<- prepared
        }
        cache[[compound]]
    }
})

for (compound in c("tryptophan", "kynurenine")) local({
    name <- compound
    test_that(paste(name, "links descendants before MS2 purity exclusion"), {
        skip_if_not_installed("xcms")
        fixture <- .msn_association_fixture(name)
        ordinary <- frag4feature(fixture$pa, fixture$xe, ppm = 10,
            convert2RawRT = FALSE)
        expect_true(all(Spectra::msLevel(ordinary@fragSpectra) == 2L))
        expect_null(ordinary@params$msnLineage)
        full <- frag4feature(fixture$pa, fixture$xe, ppm = 10,
            includeMSn = TRUE, parentReference = "acquisitionNum",
            convert2RawRT = FALSE)
        selected <- full@grped_df[
            full@grped_df$grpid %in% fixture$feature, , drop = FALSE]
        expected <- if (name == "tryptophan") c(17L, 11L, 6L) else
            c(34L, 21L, 12L)
        expect_identical(tabulate(selected$msLevel - 1L, nbins = 3L), expected)
        expect_true(all(is.na(selected$inPurity[selected$msLevel > 2L])))
        expected <- full@grped_df[full@grped_df$msLevel > 2L, , drop = FALSE]
        rejected <- frag4feature(fixture$pa, fixture$xe, ppm = 10, plim = 1,
            includeMSn = TRUE, parentReference = "acquisitionNum",
            convert2RawRT = FALSE)
        actual <- rejected@grped_df
        rownames(actual) <- rownames(expected) <- NULL
        expect_identical(actual, expected)
        expect_true(all(Spectra::msLevel(rejected@fragSpectra) > 2L))
        history <- processingHistory(rejected, "frag4feature")
        parameters <- history[[length(history)]]$parameters
        expect_true(parameters$includeMSn)
        expect_identical(parameters$resolvedParentReference,
            "acquisitionNum")
        relinked <- frag4feature(rejected, fixture$xe, ppm = 10,
            convert2RawRT = FALSE)
        expect_null(relinked@params$msnLineage)
        expect_true(all(Spectra::msLevel(relinked@fragSpectra) == 2L))
        expect_gt(length(processingHistory(relinked)),
            length(processingHistory(rejected)))
    })
})

test_that("group-width association retains source acquisition times", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    linked <- frag4feature(fixture$pa, fixture$xe, ppm = 10,
        includeMSn = TRUE, parentReference = "acquisitionNum",
        convert2RawRT = FALSE, useGroup = TRUE)
    children <- linked@grped_df[linked@grped_df$msLevel > 2L, , drop = FALSE]
    source <- fixture$pa@spectra
    rows <- match(children$pid, source$pid)
    expect_gt(nrow(children), 0L)
    expect_false(anyNA(rows))
    expect_equal(children$retentionTime, Spectra::rtime(source)[rows])
    expect_true(all(is.na(children$inPurity)))
    expect_identical(linked@f4f_link_type, "group")
})

test_that("empty associations retain diagnostics without reading peaks", {
    skip_if_not_installed("xcms")
    fixture <- .msn_association_fixture("tryptophan")
    source <- .ancestry_spectra()[2:3]
    source$dataOrigin <- rep(unname(fixture$pa@fileList)[1L], 2L)
    source$fileid <- rep(1L, 2L)
    source$seqNum <- source$scanIndex
    pa <- methods::new("purityA", fileList = fixture$pa@fileList, cores = 1)
    pa@spectra <- source
    local_mocked_bindings(.pa_raw_scans = function(...) {
        stop("Empty associations must not request peak arrays")
    }, .package = "msPurity")
    higher_only <- frag4feature(pa, fixture$xe, includeMSn = TRUE,
        parentReference = "sourceSpectrumId", convert2RawRT = FALSE)
    unmatched <- frag4feature(fixture$pa, fixture$xe, ppm = 0,
        includeMSn = TRUE, parentReference = "acquisitionNum",
        convert2RawRT = FALSE)
    cases <- list(higher_only = higher_only, unmatched = unmatched)
    for (name in names(cases)) {
        out <- cases[[name]]
        expect_identical(nrow(out@grped_df), 0L, info = name)
        expect_true(all(c("pid", "grpid", "sample", "inPurity") %in%
                    names(out@grped_df)), info = name)
        expect_identical(length(out@fragSpectra), 0L, info = name)
        expect_identical(length(groupedSpectra(out)), 0L, info = name)
        record <- tail(processingHistory(out, "frag4feature"), 1L)[[1L]]
        expect_identical(record$status, "completed", info = name)
        expect_identical(nrow(record$output$spectra), 0L, info = name)
    }
    expect_identical(nrow(higher_only@params$msnLineage$scans), 2L)
    expect_true(all(is.na(higher_only@params$msnLineage$scans$root_pid)))
    expect_true(any(unmatched@params$msnLineage$scans$path_status == "valid"))
    settings <- tail(processingParameters(unmatched, "frag4feature"), 1L)[[1L]]
    expect_identical(settings$ppm, 0)
})
