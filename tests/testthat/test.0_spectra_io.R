context("Raw data read through Spectra")

raw_files <- c(
    list.files(system.file("extdata", "lcms", "mzML", package = "msPurityData"),
               full.names = TRUE),
    list.files(system.file("extdata", "dims", "mzML", package = "msPurityData"),
               full.names = TRUE))

test_that("the header adapter equals mzR::header()", {
    for (f in raw_files) {
        mr <- mzR::openMSfile(f)
        expect_identical(msPurity:::.msp_header(msPurity:::.msp_read(f)),
                         mzR::header(mr), label = basename(f))
    }
})

test_that("peaks equal mzR::peaks()", {
    for (f in raw_files) {
        mr <- mzR::openMSfile(f)
        expect_identical(msPurity:::.msp_peaks(msPurity:::.msp_read(f)),
                         mzR::peaks(mr), label = basename(f))
    }
})

test_that("getmrdf() and getscans() accept paths or Spectra", {
    f <- raw_files[grepl("LCMSMS", raw_files)]
    sp <- msPurity:::.msp_read(f)
    expect_identical(msPurity:::getmrdf(sp), msPurity:::getmrdf(f))
    expect_identical(msPurity:::getscans(sp), msPurity:::getscans(f))
    expect_identical(msPurity:::getscans(f[1]),
                     msPurity:::.msp_peaks(msPurity:::.msp_read(f[1])))
})

test_that("isolation windows agree between file-backed and memory spectra", {
    f <- raw_files[grepl("LCMSMS_1", raw_files)]
    sp <- Spectra::filterMsLevel(msPurity:::.msp_read(f), 2L)
    recorded <- msPurity:::.msp_isolation_settings(sp)
    expect_identical(unique(recorded$lowerOffset), 0.5)
    expect_identical(unique(recorded$upperOffset), 0.5)
    mem <- Spectra::setBackend(sp, Spectra::MsBackendMemory())
    expect_identical(msPurity:::.msp_isolation_settings(mem), recorded)
})

test_that("a non-default mzRback is reported as deprecated", {
    expect_message(msPurity:::.msp_deprecate_mzRback("ramp"), "deprecated")
    expect_silent(msPurity:::.msp_deprecate_mzRback("pwiz"))
})
