## ------------------------------------------------------------------------
library(msPurity)
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)

## ------------------------------------------------------------------------
pa  <- purityA(msmsPths, interpol = "linear")
pa <- frag4feature(pa, xset)

## ------------------------------------------------------------------------
result <- spectral_matching_lcmsms(pa, xset, out_dir = tempdir())

