## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
library(msPurity)
mzMLpths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)
xset <- xcms::xcmsSet(mzMLpths)
xset <- xcms::group(xset)

