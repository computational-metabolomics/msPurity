## ----include=FALSE---------------------------------------------------------
library_database_schema_path <- system.file("extdata", "schemas", "library_database_schema.html", package="msPurity")
query_database_schema_path <- system.file("extdata", "schemas", "query_database_schema.html", package="msPurity")

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
library(msPurity)
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
pa  <- purityA(msmsPths, interpol = "linear")
pa <- frag4feature(pa, xset, create_db=TRUE)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
#result <- spectral_matching( pa@db_path, out_dir = tempdir())

