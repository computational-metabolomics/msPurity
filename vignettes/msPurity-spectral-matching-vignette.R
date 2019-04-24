## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
library(msPurity)
mzMLpths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)
xset <- xcms::xcmsSet(mzMLpths)
xset <- xcms::group(xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
pa  <- purityA(mzMLpths)
pa <- frag4feature(pa, xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
pa <- filterFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
pa <- averageIntraFragSpectra(pa) # use parameters specific to intra spectra 
pa <- averageInterFragSpectra(pa) # use parameters specific to inter spectra

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
pa <- averageAllFragSpectra(pa) 

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
q_dbPth <- createDatabase(pa, xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE------------
q_xcmsGroups = c(17, 41) # only search chosen xcms groups
q_spectraTypes = 'av_all' # only search spectra averaged using averageAllFragSpectra 
l_accessions=c('CCMSLIB00000577898', 'CE000616') # only search against these 2 accessions based from MoNA
result <- spectralMatching(q_dbPth, q_xcmsGroups = q_xcmsGroups, q_spectraTypes = q_spectraTypes,
                             cores = 1, q_pol = NA, l_accessions = l_accessions,
                             l_pol = 'positive', updateDb = FALSE)

