test_that("msfr peak lists are read through Spectra", {
    files <- list.files(system.file("extdata", "dims", "msfr-peaks",
                                    package = "msPurityData"),
                        full.names = TRUE)
    for (f in files) {
        sp <- msPurity:::.msp_read_msfr(f)
        csv <- utils::read.csv(f)
        expect_s4_class(sp, "Spectra")
        expect_length(sp, length(unique(csv$scanid)))
        expect_true(all(c("snr", "noise") %in% Spectra::peaksVariables(sp)))
        expect_identical(msPurity:::.msp_msfr_frame(sp), csv)
    }
    expect_error(msPurity:::.msp_read_msfr(
        system.file("extdata", "tests", "external_annotations", "beams.tsv",
                    package = "msPurity")), "msfr peak list")
})
