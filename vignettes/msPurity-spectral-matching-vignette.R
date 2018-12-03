## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
library(msPurity)
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
pa  <- purityA(msmsPths)
pa <- frag4feature(pa, xset, create_db=TRUE)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
result <- spectral_matching(pa@db_path,  scan_ids = c(1120,  366, 1190, 601,  404,1281, 1323, 1289))

## ---- echo=FALSE-----------------------------------------------------------
htmltools::includeHTML("query_sqlite_schema.html")


## ---- echo=FALSE-----------------------------------------------------------
htmltools::includeHTML("library_sqlite_schema.html")


