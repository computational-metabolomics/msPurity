.d05_settings_spectra <- function(origin = "resolved-memory",
    offsets = c(0.5, 0.5)) {
    target <- 100.1234
    d <- S4Vectors::DataFrame(
        msLevel = c(1L, 2L, 1L), rtime = c(1, 2, 3),
        scanIndex = 1:3, acquisitionNum = 1:3,
        precScanNum = c(NA_integer_, 1L, NA_integer_),
        precursorMz = c(NA_real_, target, NA_real_),
        precursorIntensity = c(NA_real_, 100, NA_real_),
        isolationWindowTargetMz = c(NA_real_, target, NA_real_),
        isolationWindowLowerMz = c(NA_real_, target - offsets[1], NA_real_),
        isolationWindowUpperMz = c(NA_real_, target + offsets[2], NA_real_),
        dataOrigin = rep(origin, 3), centroided = TRUE)
    d$mz <- list(c(99.9000, target, 100.3000),
        c(57.0344, 70.0456), c(99.9000, target, 100.3000))
    d$intensity <- list(c(10, 100, 40), c(10, 20), c(10, 90, 35))
    Spectra::Spectra(d, source = Spectra::MsBackendMemory())
}

test_that("history reuses each file's resolved offsets", {
    first <- .d05_settings_spectra()
    second <- .d05_settings_spectra("second-memory", c(0.4, 0.6))
    calls <- 0L
    resolve <- msPurity:::.msp_isolation_settings
    local_mocked_bindings(
        .msp_isolation_settings = function(sp) {
            calls <<- calls + 1L
            resolve(sp)
        }, .package = "msPurity")

    pa <- purityA(c(first, second), interpol = "none", isotopes = FALSE)
    settings <- processingParameters(pa, "purityA")[[1]]
    resolved <- settings$resolved

    expect_identical(calls, 2L)
    expect_identical(names(resolved), c("1", "2"))
    expect_equal(resolved[[1]]$offsets, c(0.5, 0.5))
    expect_equal(resolved[[2]]$offsets, c(0.4, 0.6))
    expect_identical(vapply(resolved, `[[`, "", "offsets_source"),
        c("1" = "recorded_metadata", "2" = "recorded_metadata"))
    expect_identical(resolved[[2]]$dataOrigin, "second-memory")
    expect_identical(resolved[[2]]$fileid, 2L)
    expect_identical(pa@params$purityA$offsets, NA)
    expect_equal(purityTable(pa)$aPurity, rep(2 / 3, 2))
    expect_null(resolved[[1]]$iwNormFun)
    expect_identical(resolved[[1]]$iwNormFun_source, "not_used")
})

test_that("caller offsets bypass resolution and public output is unchanged", {
    sp <- .d05_settings_spectra()
    callback <- function(x) rep(1, length(x))
    environment(callback) <- baseenv()
    local_mocked_bindings(
        .msp_isolation_settings = function(sp) {
            stop("The explicit offsets must not be resolved again.")
        }, .package = "msPurity")

    pa <- purityA(sp, offsets = c(0.4, 0.6), interpol = "none",
        isotopes = FALSE, iwNorm = TRUE, iwNormFun = callback)
    resolved <- processingParameters(pa, "purityA")[[1]]$resolved[[1]]
    public <- assessPuritySingle(sp, offsets = c(0.4, 0.6),
        interpol = "none", isotopes = FALSE, iwNormFun = callback)
    detailed <- msPurity:::.pa_assess_purity_single(sp,
        offsets = c(0.4, 0.6), interpol = "none", isotopes = FALSE,
        iwNormFun = callback)

    expect_identical(resolved$offsets_source, "caller")
    expect_identical(resolved$offsets_method, "explicit_offsets")
    expect_identical(resolved$iwNormFun, callback)
    expect_identical(resolved$iwNormFun_source, "caller")
    expect_identical(public, detailed$purity)
    expect_s3_class(public, "data.frame")
    expect_false("resolved" %in% names(attributes(public)))
    expect_equal(public$aPurity, 2 / 3)
})

test_that("parallel assessment retains each file's resolved settings", {
    input <- c(.d05_settings_spectra(),
        .d05_settings_spectra("parallel-memory", c(0.4, 0.6)))
    serial <- purityA(input, cores = 1L, interpol = "none", isotopes = FALSE)
    parallel <- purityA(input, cores = 2L, interpol = "none", isotopes = FALSE)

    expect_identical(purityTable(parallel), purityTable(serial))
    expect_identical(
        processingParameters(parallel, "purityA")[[1]]$resolved,
        processingParameters(serial, "purityA")[[1]]$resolved)
})

test_that("applied normalization survives RDS without peak arrays", {
    pa <- purityA(.d05_settings_spectra(), interpol = "none",
        isotopes = FALSE, iwNorm = TRUE)
    parameters <- processingParameters(pa, "purityA")[[1]]
    resolved <- parameters$resolved[[1]]
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))

    saveRDS(parameters, path)
    restored <- readRDS(path)$resolved[[1]]

    expect_identical(resolved$iwNormFun_source, "generated_default")
    expect_true(resolved$iwNorm)
    expect_false(any(c("scans", "sp", "mrdf") %in%
                ls(environment(resolved$iwNormFun))))
    expect_equal(restored$iwNormFun(c(-0.5, 0, 0.5)), c(0, 1, 0))
    expect_identical(restored$offsets, resolved$offsets)
    expect_null(parameters$iwNormFun)
})

test_that("file settings distinguish metadata from the Agilent assumption", {
    path <- tempfile(fileext = ".xml")
    before <- getAllConnections()
    on.exit({
        for (id in setdiff(getAllConnections(), before))
            base::close(getConnection(id))
        unlink(path)
    })
    writeLines(c(
        '<cvParam name="isolation window lower offset" value="0.4"/>',
        '<cvParam name="isolation window upper offset" value="0.6"/>'), path)
    recorded <- msPurity:::.pa_read_isolation_settings(path)
    expect_identical(recorded$offsets_source, "recorded_metadata")
    expect_identical(recorded$offsets_method, "mzml_first_window")
    expect_identical(msPurity:::get_isolation_offsets(path), c(0.4, 0.6))

    writeLines(paste0('<cvParam cvRef="MS" accession="MS:1000490" ',
            'name="Agilent instrument model" value=""/>'), path)
    expect_message(
        approximation <- msPurity:::.pa_read_isolation_settings(path),
        "approximate")
    expect_identical(approximation$offsets, c(0.65, 0.65))
    expect_identical(approximation$offsets_source, "approximation")
    expect_identical(approximation$offsets_method, "agilent_narrow_window")
})

test_that("MS1-only assessment retains NULL and unused settings", {
    sp <- Spectra::filterMsLevel(.d05_settings_spectra(), 1L)
    local_mocked_bindings(
        .msp_isolation_settings = function(sp) {
            stop("Unused offsets must not be resolved.")
        }, .package = "msPurity")

    expect_message(public <- assessPuritySingle(sp), "No MS/MS spectra")
    expect_message(detailed <- msPurity:::.pa_assess_purity_single(sp),
        "No MS/MS spectra")

    expect_null(public)
    expect_null(detailed$purity)
    expect_identical(detailed$resolved$status, "no_ms2")
    expect_identical(detailed$resolved$offsets_source, "not_used")
})

test_that("MS2-only assessment retains unavailable values", {
    sp <- Spectra::filterMsLevel(.d05_settings_spectra(), 2L)
    expect_message(public <- assessPuritySingle(sp), "only MS2 data")
    expect_message(detailed <- msPurity:::.pa_assess_purity_single(sp),
        "only MS2 data")

    expect_identical(public, detailed$purity)
    expect_true(all(is.na(public$inPurity)))
    expect_identical(detailed$resolved$status, "no_ms1")
    expect_identical(detailed$resolved$iwNormFun_source, "not_used")
})
