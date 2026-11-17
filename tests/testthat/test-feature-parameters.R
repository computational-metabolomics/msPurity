test_that("feature parameters are one complete serializable configuration", {
    param <- Frag4FeatureParam()
    expect_s4_class(param, "Param")
    expect_true(methods::validObject(param))
    expect_identical(as.list(param), list(ppm = 5, plim = NA_real_,
            intense = TRUE, convert2RawRT = TRUE, useGroup = FALSE,
            parentReference = "auto", rtTol = 0, rtSource = "precursor",
            groupWindow = "full", includeFilled = FALSE,
            estimateRtMap = FALSE))
    expect_identical(as.list(Frag4FeatureParam(plim = NA))$plim, NA_real_)

    complete <- Frag4FeatureParam(ppm = 10, plim = 0.7, intense = FALSE,
        convert2RawRT = FALSE, useGroup = TRUE, parentReference = "scanIndex",
        rtTol = 2, rtSource = "ms2", groupWindow = "median",
        includeFilled = TRUE, estimateRtMap = TRUE)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(complete, path)
    expect_identical(readRDS(path), complete)
    expect_identical(as.list(readRDS(path)), list(ppm = 10, plim = 0.7,
            intense = FALSE, convert2RawRT = FALSE, useGroup = TRUE,
            parentReference = "scanIndex", rtTol = 2, rtSource = "ms2",
            groupWindow = "median", includeFilled = TRUE,
            estimateRtMap = TRUE))
})

test_that("feature parameter validity rejects invalid or removed settings", {
    cases <- list(negative_ppm = list(ppm = -1),
        infinite_ppm = list(ppm = Inf), multiple_ppm = list(ppm = c(5, 10)),
        negative_plim = list(plim = -0.1), excessive_plim = list(plim = 1.1),
        infinite_plim = list(plim = Inf), intense = list(intense = NA),
        raw_rt = list(convert2RawRT = 1), groups = list(useGroup = logical()),
        reference = list(parentReference = "other"),
        negative_tol = list(rtTol = -1), infinite_tol = list(rtTol = Inf),
        source = list(rtSource = "survey"), window = list(groupWindow = "max"),
        filled = list(includeFilled = NA), estimate = list(estimateRtMap = 1))
    for (name in names(cases))
        expect_error(do.call(Frag4FeatureParam, cases[[name]]),
            names(cases[[name]]), info = name)
    expect_error(Frag4FeatureParam(select = list(msLevel = 4L)), "unused")
    expect_error(Frag4FeatureParam(includeMSn = TRUE), "unused")
    expect_error(Frag4FeatureParam(msLevel = 4L), "unused")
    invalid <- Frag4FeatureParam()
    invalid@ppm <- -1
    expect_error(methods::validObject(invalid), "ppm")
    expect_error(frag4feature(methods::new("purityA"), NULL,
            param = invalid), "ppm")
})

test_that("association accepts one parameter object and rejects conflicts", {
    pa <- methods::new("purityA")
    param <- Frag4FeatureParam()
    cases <- list(unnamed_list = list(param),
        named_list = list(default = param), wrong_class = PurityAParam())
    for (name in names(cases))
        expect_error(frag4feature(pa, NULL, param = cases[[name]]),
            "one Frag4FeatureParam", info = name)
    conflicts <- list(ppm = 10, plim = 0.7, intense = FALSE,
        convert2RawRT = FALSE, useGroup = TRUE, parentReference = "scanIndex",
        rtTol = 2, rtSource = "ms2", groupWindow = "median",
        includeFilled = TRUE, estimateRtMap = TRUE)
    for (name in names(conflicts))
        expect_error(do.call(frag4feature,
                c(list(pa = pa, xcmsObj = NULL, param = param),
                    conflicts[name])), paste0(name, " conflicts"), info = name)
    expect_error(frag4feature(pa, NULL, param = param, use_group = TRUE),
        "use_group conflicts")
    expect_error(frag4feature(pa, NULL,
            param = Frag4FeatureParam(useGroup = TRUE), useGroup = FALSE,
            use_group = TRUE), "useGroup conflicts")
})
