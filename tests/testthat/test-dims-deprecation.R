test_that("DIMS entry points warn before processing inputs", {
    pd <- .with_dims_deprecation(purityD())
    calls <- list(
        purityD = quote(purityD()),
        Getfiles = quote(Getfiles()),
        averageSpectra = quote(averageSpectra(pd)),
        averageSpectraSingle = quote(averageSpectraSingle()),
        filterp = quote(filterp(pd)),
        subtract = quote(msPurity::subtract(pd)),
        subtractMZ = quote(subtractMZ()),
        dimsPredictPurity = quote(dimsPredictPurity(pd)),
        dimsPredictPuritySingle = quote(dimsPredictPuritySingle()),
        groupPeaks = quote(groupPeaks(pd)),
        groupPeaksEx = quote(groupPeaksEx()),
        writeOut = quote(writeOut(pd)),
        get_additional_mzml_meta = quote(get_additional_mzml_meta()),
        averagedPeaks = quote(averagedPeaks(pd)),
        getP = quote(getP(pd)),
        updatePeaks = quote(msPurity:::updatePeaks(pd, list())),
        updateObject = quote(BiocGenerics::updateObject(pd)),
        show = quote(methods::show(pd)))

    for (name in names(calls)) {
        # Catch the first condition before file access or a worker is started.
        warning <- tryCatch(eval(calls[[name]]),
                            deprecatedWarning = function(w) w)
        expect_true(inherits(warning, "deprecatedWarning"), info = name)
        expect_identical(warning$old, name, info = name)
        expect_identical(warning$package, "msPurity", info = name)
        expect_match(conditionMessage(warning),
                     "no further bug fixes or feature development", info = name)
        expect_match(conditionMessage(warning),
                     "no direct replacement in msPurity", info = name)
        expect_match(conditionMessage(warning), "msPurity-deprecated", info = name)
    }
})

test_that("standalone DIMS prediction still calculates purity", {
    path <- tempfile(fileext = ".csv")
    on.exit(unlink(path))
    write.csv(data.frame(scanid = rep(1:2, each = 2),
        mz = rep(c(100, 100.1), 2), i = rep(c(80, 20), 2)),
        path, row.names = FALSE)

    result <- expect_warning(dimsPredictPuritySingle(100, path,
        mzML = FALSE, isotopes = FALSE, ilim = 0),
        "'dimsPredictPuritySingle' is deprecated", class = "deprecatedWarning")

    expect_named(result, c("medianPurity", "meanPurity", "sdPurity",
                          "cvPurity", "sdePurity", "medianPeakNum"))
    expect_equal(unname(unlist(result)), c(0.8, 0.8, 0, 0, 0, 2))
})

test_that("saved DIMS results remain accessible and exportable", {
    files <- data.frame(filepth = "sample.mzML", name = "sample",
                        sampleType = "sample", class = "sample")
    peaks <- data.frame(peakID = 1:2, mz = c(100, 200), i = c(80, 20))
    expected <- list(orig = list(sample = peaks), processed = list(sample = peaks))
    path <- tempfile(fileext = ".rds")
    output <- tempfile("dims-export-")
    on.exit(unlink(c(path, output), recursive = TRUE))

    .with_dims_deprecation({
        pd <- purityD(files)
        pd <- msPurity:::updatePeaks(pd, expected)
        saveRDS(pd, path)
        saved <- readRDS(path)
        current <- BiocGenerics::updateObject(saved)
        expect_equal(averagedPeaks(current), expected)
        expect_equal(getP(current), expected)
        exported <- writeOut(current, output, original = TRUE)
        expect_equal(averagedPeaks(exported), expected)
    })

    for (i in seq_len(nrow(saved@fileList))) {
        sample <- as.character(saved@fileList$name[i])
        for (stage in c("orig", "processed")) {
            suffix <- if (stage == "orig") "original" else "processed"
            path <- file.path(output, paste0(sample, "_", suffix, ".csv"))
            actual <- read.csv(path)
            reference <- expected[[stage]][[i]]
            rownames(reference) <- NULL
            expect_equal(actual, reference, info = paste(sample, stage))
        }
    }
})

test_that("DIMS warning handling does not hide other warnings", {
    expect_warning(.with_dims_deprecation(warning("unrelated warning")),
                   "unrelated warning")
    expect_warning(.with_dims_deprecation(.Deprecated(old = "anotherAPI",
        package = "msPurity")), "anotherAPI", class = "deprecatedWarning")
})

test_that("shared LC-MS interfaces do not emit DIMS deprecations", {
    expect_warning(pa <- methods::new("purityA"), NA)
    expect_warning(BiocGenerics::updateObject(pa), NA)
    expect_warning(purityTable(pa), NA)
    expect_warning(averagedSpectra(pa), NA)
    expect_warning(iwNormGauss(), NA)
    expect_warning(iwNormRcosine(), NA)
    peaks <- cbind(mz = c(100, 100.1), intensity = c(80, 20))
    expect_warning(value <- pcalc(peaks, mzmin = 99.5, mzmax = 100.5,
        mztarget = 100, ppm = 5, iwNorm = FALSE, ilim = 0, isotopes = FALSE), NA)
    expect_equal(unname(value[1]), 0.8)
})
