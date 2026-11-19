.px_calculation_fixture <- function() {
    intensities <- c(10, 40, 1000, 80, 40, 10)
    arrays <- lapply(intensities, function(i) {
        cbind(mz = c(205.0972, 205.2972), intensity = c(i, 20))
    })
    tracking <- data.frame(scan = seq_along(arrays), sample = 1L,
        msLevel = c(1L, 1L, 2L, 1L, 1L, 1L),
        acquisitionNum = 100L + seq_along(arrays),
        retentionTime = seq_along(arrays), rt_corrected = seq_along(arrays))
    target <- data.frame(sample = 1L, id = 1, output_id = 1, c_peak_id = 1,
        mz = 205.0972, mzmin = 205.0970, mzmax = 205.0974,
        rt = 4, rtmin = 1, rtmax = 6, minscan = 1L, maxscan = 6L)
    settings <- list(offsets = c(0.5, 0.5), iwNorm = FALSE,
        iwNormFun = NULL, ilim = 0, isotopes = FALSE, im = NULL)
    list(target = target, prepared = list(scanpeaks = list(arrays),
            tracking = tracking), settings = settings)
}

test_that("purityX selects MS1 before peak-width and purity summaries", {
    x <- .px_calculation_fixture()
    result <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_identical(result$trace$scan, c(1L, 2L, 4L, 5L, 6L))
    expect_equal(result$fwhm, c(2, 5))
    expect_equal(unname(result$metrics[c("purityFWHMmedian", "purityFWmean")]),
        c(2 / 3, mean(c(1 / 3, 2 / 3, 0.8, 2 / 3, 1 / 3))))
    expect_equal(result$trace$isolationIntensity, c(30, 60, 100, 60, 30))
    expect_equal(result$trace$intensity, c(10, 40, 80, 40, 10))

    x$target$scan_positions <- I(list(c(2L, 4L, 5L)))
    subset <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_identical(subset$observations$scan, c(2L, 4L, 5L))

    x$target$scan_positions <- I(list(c(5L, 4L, 2L)))
    reordered <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_identical(reordered$metrics, subset$metrics)
    expect_identical(reordered$observations, subset$observations)
})

test_that("purityX retains full-width estimates without an FWHM", {
    x <- .px_calculation_fixture()
    x$target$minscan <- x$target$maxscan <- 4L
    one <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_true(is.na(one$metrics[["purityFWHMmedian"]]))
    expect_equal(one$metrics[["purityFWmedian"]], 0.8)
    expect_identical(one$fwhm_status, "missing_flank_or_signal")
    x$prepared$scanpeaks[[1L]][[4L]] <- matrix(numeric(), ncol = 2L)
    empty <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_true(all(is.na(empty$metrics)))
    expect_identical(empty$status, "no_usable_purity")
    x$target$minscan <- NA_integer_
    missing <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_identical(missing$status, "invalid_rt_bounds")
    expect_equal(nrow(missing$observations), 0L)
    x$target$observation_status <- "no_ms1_observations"
    absent <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_identical(absent$status, "no_ms1_observations")
})

test_that("purityX handles one peak and asymmetric isolation windows", {
    x <- .px_calculation_fixture()
    x$target$minscan <- x$target$maxscan <- 4L
    x$prepared$scanpeaks[[1L]][[4L]] <- matrix(c(205.0972, 80), nrow = 1L)
    one <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_equal(one$metrics[["purityFWmedian"]], 1)
    x$prepared$scanpeaks[[1L]][[4L]] <- cbind(
        mz = c(204.4972, 205.0972, 205.6972), intensity = c(50, 80, 20))
    x$settings$offsets <- c(0.7, 0.3)
    asymmetric <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_equal(asymmetric$metrics[["purityFWmedian"]], 0.8)
    expect_equal(asymmetric$trace$isolationIntensity, 100)

    x$settings$offsets <- c(0, 0.7)
    boundary <- .px_predict_one(x$target, x$prepared, x$settings)
    expect_equal(with(boundary$trace,
            c(intensity, isolatedTargetIntensity,
                isolationIntensity - isolatedTargetIntensity)), c(80, 0, 50))
    expect_true(is.na(boundary$metrics[["purityFWmedian"]]))
    expect_identical(boundary$status, "no_usable_purity")
})

test_that("purityX preserves peak choices and their actual observations", {
    x <- .px_calculation_fixture()
    targets <- x$target[rep(1L, 3L), ]
    targets$grpid <- 1
    targets$c_peak_id <- c(8, 9, 10)
    targets$sample <- c(1L, 1L, 2L)
    targets$rt <- c(3, 5, 4)
    selected <- .px_choose_targets(targets, single = FALSE)
    expect_equal(selected$targets$c_peak_id, c(8, 10))
    expect_equal(selected$excluded$c_peak_id, 9)
    result <- .px_predict_one(x$target, x$prepared, x$settings)
    refs <- .px_calculation_refs(list(result), selected$excluded)
    expect_identical(refs$observations$scan, c(1L, 2L, 4L, 5L, 6L))
    expect_equal(refs$observations$c_peak_id, rep(1, 5))
    expect_false(any(c("mz", "intensity") %in% names(refs$observations)))
})

test_that("purityX preserves default tables and adds requested means", {
    x <- .px_calculation_fixture()
    result <- .px_predict_one(x$target, x$prepared, x$settings)
    default <- .px_prediction_table(list(result), data.frame(),
        "purityFWHMmedian", single = TRUE)
    expect_identical(names(default), c("id", "purityFWHMmedian",
            "purityFWmedian", "pknmFWHMmedian", "pknmFWmedian"))
    mean <- .px_prediction_table(list(result), data.frame(),
        "purityFWmean", single = TRUE)
    expect_identical(names(mean), c(names(default), "purityFWmean"))
    expect_equal(mean$purityFWmean, mean(c(1 / 3, 2 / 3, .8, 2 / 3, 1 / 3)))
    empty <- .px_prediction_table(list(), data.frame(),
        "purityFWHMmedian", single = TRUE)
    expect_identical(names(empty), names(default))
    expect_equal(nrow(empty), 0L)
})

test_that("purityX summarises peak occurrences in feature order", {
    groups <- data.frame(grpid = c(9L, 2L, 7L), i = 100, mz = 205.0972)
    predictions <- Map(function(group, purity, count) {
        list(target = list(grpid = group),
            metrics = c(purityFWmedian = purity, pknmFWmedian = count))
    }, c(7L, 9L, 7L), c(0.4, 0.9, 0.8), c(2, 1, 4))
    table <- .px_prediction_table(predictions, groups, "purityFWmedian", FALSE)
    expect_equal(table$grpid, groups$grpid)
    expect_equal(table$mean, c(0.9, NA, 0.6))
    expect_equal(table$pknm, c(1, NA, 3))
})
