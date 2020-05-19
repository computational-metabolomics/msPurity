## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
library(msPurity)
mzMLpths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)
xset <- xcms::xcmsSet(mzMLpths)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
pa  <- purityA(mzMLpths)
pa <- frag4feature(pa, xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
pa <- filterFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
pa <- averageIntraFragSpectra(pa) # use parameters specific to intra spectra 
pa <- averageInterFragSpectra(pa) # use parameters specific to inter spectra

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
pa <- averageAllFragSpectra(pa) 

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
td <- tempdir()
q_dbPth <- createDatabase(pa, xset, outDir = td, dbName = 'lcmsms-processing.sqlite')

## -----------------------------------------------------------------------------
result <- spectralMatching(q_dbPth, q_xcmsGroups = c(17, 41), l_accessions=c('CCMSLIB00000577898','CE000616'))
print(result)

