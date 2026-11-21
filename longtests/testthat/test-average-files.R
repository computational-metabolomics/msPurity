test_that("uniform groups reproduce intra- and inter-file averages", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
        "3_filterFragSpectra_pa.rds", package = "msPurity"))
    param <- AverageFragSpectraParam(
        combine = c("activationType", "collisionEnergy"))
    scalar <- averagedSpectra(averageInterFragSpectra(
        averageIntraFragSpectra(pa)), legacy = FALSE)
    typed <- averagedSpectra(averageInterFragSpectra(
        averageIntraFragSpectra(pa, param = param), param = param),
        legacy = FALSE)
    key <- function(x) paste(x$av_level, x$grpid, x$sample)
    produced <- function(x) x[lengths(Spectra::mz(x)) > 0L]
    scalar <- produced(scalar)
    typed <- produced(typed)

    expect_setequal(key(typed), key(scalar))
    peaks <- function(x) Spectra::peaksData(x,
        columns = Spectra::peaksVariables(scalar))
    for (k in key(scalar)) {
        a <- scalar[key(scalar) == k]
        b <- typed[key(typed) == k]
        expect_equal(peaks(b)[[1L]], peaks(a)[[1L]], info = k)
        expect_identical(b$msPurityContributors[[1L]],
            a$msPurityContributors[[1L]], info = k)
        expect_identical(b$msPurityClusterContributors[[1L]],
            a$msPurityClusterContributors[[1L]], info = k)
    }
})

test_that("averages with and without param coexist at different levels", {
    pa <- readRDS(system.file("extdata", "tests", "purityA",
        "3_filterFragSpectra_pa.rds", package = "msPurity"))
    param <- AverageFragSpectraParam(combine = "collisionEnergy")
    peaks_of <- function(pa, level) {
        x <- averagedSpectra(pa, legacy = FALSE)
        x <- x[x$av_level == level]
        stats::setNames(as.list(Spectra::peaksData(x)),
            paste(x$grpid, x$sample))
    }

    scalar_all <- averageAllFragSpectra(pa)
    typed_intra <- averageIntraFragSpectra(scalar_all, param = param)
    expect_identical(peaks_of(typed_intra, "av_all"),
        peaks_of(scalar_all, "av_all"))

    scalar_levels <- averageInterFragSpectra(averageIntraFragSpectra(pa))
    typed_all <- averageAllFragSpectra(scalar_levels, param = param)
    for (level in c("av_intra", "av_inter"))
        expect_identical(peaks_of(typed_all, level),
            peaks_of(scalar_levels, level), info = level)
})
