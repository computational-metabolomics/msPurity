test_that("chunk boundaries keep compatible fragments together", {
    # The former 5000-row split divided the 51 occurrences of m/z 198.
    blocks <- msPurity:::.hc_blocks(rep(100:199, each = 51), 5, 0)
    expect_equal(blocks$first, c(1L, 4999L))
    expect_equal(blocks$last, c(4998L, 5100L))

    # Exercise the full worker with small blocks, avoiding a large matrix.
    plan_blocks <- msPurity:::.hc_blocks
    for (tolerance in c(0, 0.0002)) {
        observed <- with_mocked_bindings(
            msPurity:::clustering(rep(c(100, 200), each = 3), ppm = 5,
                tolerance = tolerance, cores = 1),
            .hc_blocks = function(mz, ppm, tolerance)
                plan_blocks(mz, ppm, tolerance, max_size = 4L),
            .package = "msPurity")
        expect_identical(observed, rep(1:2, each = 3),
            info = paste("tolerance", tolerance))
    }
    # A connected chain must still follow complete-linkage clustering.
    expect_identical(msPurity:::clustering(c(100, 100.0001, 100.0006),
        ppm = 5, cores = 1), c(1L, 1L, 2L))
})

test_that("oversized runs are either one safe cluster or an explicit error", {
    observed <- msPurity:::clustering(
        seq(100, 100.0004, length.out = 5001), ppm = 5, cores = 1)
    expect_identical(observed, rep.int(1L, 5001))
    expect_error(msPurity:::clustering(100 + 0:5000 * 0.0001,
        ppm = 5, cores = 1), "5001 peaks exceeds the 5000-peak clustering limit")
})

test_that("compatible endpoints avoid a distance matrix at every size", {
    for (n in c(2L, 5000L, 5001L)) {
        observed <- with_mocked_bindings(
            msPurity:::clustering(seq(100, 100.0004, length.out = n),
                ppm = 5, cores = 1),
            dist = function(...) stop("Unexpected distance matrix"),
            .package = "stats")
        expect_identical(observed, rep.int(1L, n), info = n)
    }
    expect_identical(with_mocked_bindings(
        msPurity:::clustering(c(100, 200), cores = 2),
        makeCluster = function(...) stop("Unexpected worker pool"),
        .package = "parallel"), 1:2)
})

test_that("parallel blocks preserve results and the registered backend", {
    foreach::setDoPar(function(...) stop("Unrelated backend"),
        info = function(data, item) switch(item, name = "existing", workers = 1L))
    on.exit(foreach::registerDoSEQ(), add = TRUE)
    plan_blocks <- msPurity:::.hc_blocks
    observed <- with_mocked_bindings(
        msPurity:::clustering(c(100, 100.0006, 200, 200.0012),
            ppm = 5, cores = 2),
        .hc_blocks = function(mz, ppm, tolerance)
            plan_blocks(mz, ppm, tolerance, max_size = 2L),
        .package = "msPurity")
    expect_identical(observed, 1:4)
    expect_identical(foreach::getDoParName(), "existing")
})
