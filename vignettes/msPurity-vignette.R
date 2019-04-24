## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
library(msPurity)
mzMLpths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")

## --------------------------------------------------------------------------
pa <- purityA(mzMLpths)
print(pa@puritydf[1:3,])

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
pa_norm <- purityA(msmsPths[1], iwNorm=TRUE, iwNormFun=iwNormGauss(sdlim=3, minOff=-0.5, maxOff=0.5))

## ----results='hide', message=FALSE, warning=FALSE,  echo = T---------------

suppressPackageStartupMessages(library(xcms))

xset <- xcms::xcmsSet(mzMLpths)
xset <- xcms::group(xset)

