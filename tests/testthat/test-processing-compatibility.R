test_that("metadata compatibility uses public Annotated accessors", {
    x <- S4Vectors::SimpleList(value = 1L)
    get_metadata <- msPurity:::.pa_spectra_metadata
    set_metadata <- get(".pa_spectra_metadata<-", asNamespace("msPurity"))
    expect_true(methods::hasMethod("metadata", class(x),
        where = asNamespace("S4Vectors")))
    expected <- list(source = "public-accessor", settings = list(ppm = 5))
    x <- set_metadata(x, expected)
    expect_identical(S4Vectors::metadata(x), expected)
    expect_identical(get_metadata(x), expected)
})

test_that("metadata compatibility round-trips the supported Spectra layout", {
    x <- Spectra::Spectra(Spectra::MsBackendMemory())
    expected <- list(source = "spectra", settings = list(ppm = 5))
    set_metadata <- get(".pa_spectra_metadata<-", asNamespace("msPurity"))
    x <- set_metadata(x, expected)
    expect_identical(msPurity:::.pa_spectra_metadata(x), expected)
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(x, path)
    expect_identical(msPurity:::.pa_spectra_metadata(readRDS(path)), expected)
    expect_error(set_metadata(x, numeric()), "list")
})
