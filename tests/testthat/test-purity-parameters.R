.purity_param_metadata <- function() {
    data <- S4Vectors::DataFrame(pid = c(7L, 40L, 41L, 8L),
        msLevel = c(2L, 3L, 4L, 2L),
        activationType = c("collision-induced dissociation", "HCD", NA, "ETD"),
        massAnalyzer = "Orbitrap", collisionEnergy = c(20, 30, NA, 20),
        collisionEnergyUnit = c("electronvolt", "normalized collision energy",
            NA, "electronvolt"), polarity = c(1L, 1L, 0L, 1L))
    data$mz <- rep(list(numeric()), 4L)
    data$intensity <- rep(list(numeric()), 4L)
    Spectra::Spectra(data, source = Spectra::MsBackendMemory())
}

.purity_param_fallback <- function() {
    list(mostIntense = FALSE, nearest = TRUE, offsets = NA_real_,
        interpol = "linear", iwNorm = FALSE, iwNormFun = NULL,
        ilim = 0.2, isotopes = FALSE, im = NULL, ppmInterp = 7)
}

test_that("purity parameters retain complete fixed defaults and functions", {
    param <- PurityAParam()
    expect_s4_class(param, "Param")
    expect_true(methods::validObject(param))
    expect_identical(as.list(param), list(select = list(), mostIntense = FALSE,
            nearest = TRUE, offsets = NA_real_, interpol = "linear",
            iwNorm = FALSE, iwNormFun = NULL, ilim = 0.05, isotopes = TRUE,
            im = NULL, ppmInterp = 7))
    selected <- as.list(PurityAParam(select = list(msLevel = 4L)))
    expect_identical(selected[-1L], as.list(param)[-1L])
    expect_identical(as.list(PurityAParam(offsets = NA))$offsets, NA_real_)

    callback <- function(x) 1 - abs(x)
    environment(callback) <- baseenv()
    model <- matrix(c(1, 1.003355, 1.07, 4, 0.1, 1, 12, 1), nrow = 1L)
    complete <- PurityAParam(iwNorm = TRUE, iwNormFun = callback, im = model)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(complete, path)
    restored <- as.list(readRDS(path))
    expect_identical(restored$iwNormFun, callback)
    expect_identical(restored$im, model)
    expect_equal(restored$iwNormFun(c(-0.2, 0, 0.2)), c(0.8, 1, 0.8))
})

test_that("purity parameter validity rejects malformed settings", {
    cases <- list(
        selector = list(select = list(parentPeakPpm = 10)),
        logical = list(mostIntense = NA),
        offsets = list(offsets = c(NA_real_, 0.5)),
        interpolation = list(interpol = "other"),
        fraction = list(ilim = 1.1),
        tolerance = list(ppmInterp = Inf),
        function_value = list(iwNormFun = 1),
        isotope_model = list(im = matrix(1, nrow = 1L, ncol = 7L)))
    errors <- c(selector = "unsupported", logical = "TRUE or FALSE",
        offsets = "Explicit isolation offsets", interpolation = "interpol",
        fraction = "between zero and one", tolerance = "finite",
        function_value = "function or NULL", isotope_model = "eight columns")
    for (name in names(cases))
        expect_error(do.call(PurityAParam, cases[[name]]), errors[[name]],
            info = name)
    expect_error(PurityAParam(parentPeakPpm = 10), "unused")
})

test_that("purity resolution validates complete named parameter objects", {
    source <- .purity_param_metadata()
    resolve <- function(param) {
        msPurity:::.pa_purity_resolve(source, param, .purity_param_fallback())
    }
    expect_error(resolve(list(PurityAParam())), "names")
    expect_error(resolve(list()), "names")
    expect_error(resolve(list(default = list(ilim = 0))),
        "requires PurityAParam")
    invalid <- PurityAParam()
    invalid@ilim <- -1
    expect_error(resolve(invalid), "nonnegative")
})

test_that("purity selectors reuse metadata and level-specific fallback", {
    source <- .purity_param_metadata()
    fallback <- .purity_param_fallback()
    fallback$mostIntense <- TRUE
    fallback$offsets <- c(0.2, 0.3)
    params <- list(
        cid = PurityAParam(select = list(activationType = "MS:1000133",
                collisionEnergy = 20, collisionEnergyUnit = "electronvolt")),
        unknown = PurityAParam(select = list(activationType = NA),
            ilim = 0.1, isotopes = FALSE),
        unused = PurityAParam(select = list(activationType = "HCD",
                collisionEnergyUnit = "electronvolt")))
    local_mocked_bindings(.msp_peaks = function(...) {
        stop("Parameter assignment must not read peaks")
    }, .package = "msPurity")
    resolved <- msPurity:::.pa_purity_resolve(source, params, fallback)

    expect_equal(vapply(resolved$parameters, `[[`, numeric(1L), "ilim"),
        c(0.05, 0, 0.1, 0.2))
    expect_identical(vapply(resolved$parameters, `[[`, logical(1L), "isotopes"),
        c(TRUE, FALSE, FALSE, FALSE))
    expect_identical(vapply(resolved$parameters, `[[`, logical(1L),
            "mostIntense"),
        c(FALSE, TRUE, FALSE, TRUE))
    expect_identical(resolved$parameters[[2L]]$offsets, c(0.2, 0.3))
    expect_identical(resolved$parameters[[1L]]$offsets, NA_real_)
    expect_identical(resolved$assignments$pid, source$pid)
    expect_identical(resolved$assignments$parameter_set_names,
        list("cid", character(), "unknown", character()))
    expect_identical(resolved$assignments$fallback_used,
        c(FALSE, TRUE, FALSE, TRUE))
    expect_identical(resolved$unused_parameter_sets, "unused")
    expect_identical(resolved$parameter_objects, params)
    expect_identical(resolved$parameters[[2L]]$not_applicable,
        c("nearest", "interpol", "iwNorm", "iwNormFun", "ppmInterp"))
    expect_identical(resolved$parameters[[2L]]$interpol, "none")
    expect_true(is.na(resolved$parameters[[2L]]$nearest))
    expect_null(resolved$parameters[[2L]]$iwNormFun)
    expect_true("iwNormFun" %in% names(resolved$parameters[[2L]]))
})

test_that("higher-level overlaps compare only applicable settings", {
    source <- .purity_param_metadata()
    params <- list(first = PurityAParam(),
        second = PurityAParam(nearest = FALSE, interpol = "spline",
            ppmInterp = 20))
    resolved <- msPurity:::.pa_purity_resolve(source[2:3], params,
        .purity_param_fallback())
    expect_length(resolved$parameter_sets, 1L)
    expect_identical(resolved$assignments$parameter_id,
        rep("parameter-0001", 2L))
    expect_identical(resolved$assignments$parameter_set_names,
        rep(list(c("first", "second")), 2L))
    expect_identical(as.list(resolved$parameter_objects$second)$interpol,
        "spline")
    expect_equal(resolved$parameters[[1L]]$ilim, 0.05)
    expect_true(resolved$parameters[[1L]]$isotopes)
    expect_error(msPurity:::.pa_purity_resolve(source, params,
            .purity_param_fallback()), "Conflicting overlapping.*1")
    params$second <- PurityAParam(ilim = 0.1)
    expect_error(msPurity:::.pa_purity_resolve(source[2:3], params,
            .purity_param_fallback()), "Conflicting overlapping.*1")
})

test_that("unsupported weighting is checked only for assigned higher events", {
    source <- .purity_param_metadata()
    callback <- function(x) rep(1, length(x))
    environment(callback) <- baseenv()
    unsupported <- list(enabled = PurityAParam(iwNorm = TRUE),
        callback = PurityAParam(iwNormFun = callback))
    for (name in names(unsupported))
        expect_error(msPurity:::.pa_purity_resolve(source, unsupported[[name]],
                .purity_param_fallback()), "available only for MS2",
            info = name)
    params <- list(ms2 = PurityAParam(select = list(msLevel = 2L),
            iwNorm = TRUE, iwNormFun = callback),
        unused = PurityAParam(select = list(msLevel = 5L), iwNorm = TRUE))
    fallback <- .purity_param_fallback()
    fallback$iwNorm <- TRUE
    fallback$iwNormFun <- callback
    resolved <- msPurity:::.pa_purity_resolve(source, params, fallback)
    expect_true(resolved$parameters[[1L]]$iwNorm)
    expect_identical(resolved$parameters[[1L]]$iwNormFun, callback)
    expect_false(resolved$parameters[[2L]]$iwNorm)
    expect_null(resolved$parameters[[2L]]$iwNormFun)
    expect_identical(resolved$unused_parameter_sets, "unused")

    empty <- msPurity:::.pa_purity_resolve(source[integer()], params, fallback)
    expect_identical(empty$parameters, list())
    expect_identical(empty$parameter_sets, list())
    expect_identical(empty$unused_parameter_sets, names(params))
})
