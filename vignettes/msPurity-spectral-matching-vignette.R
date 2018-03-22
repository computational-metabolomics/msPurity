## --------------------------------------------------------------------------
library(msPurity)
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)

## --------------------------------------------------------------------------
pa  <- purityA(msmsPths, interpol = "linear")
pa <- frag4feature(pa, xset, create_db=TRUE)

## --------------------------------------------------------------------------
result <- spectral_matching( pa@db_path, out_dir = tempdir())

