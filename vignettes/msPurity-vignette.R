## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
library(msPurity)
mzMLpths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")

## ------------------------------------------------------------------------
pa <- purityA(mzMLpths)
print(pa@puritydf[1:3,])

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
pa_norm <- purityA(msmsPths[1], iwNorm=TRUE, iwNormFun=iwNormGauss(sdlim=3, minOff=-0.5, maxOff=0.5))

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------

suppressPackageStartupMessages(library(xcms))

xset <- xcms::xcmsSet(mzMLpths)
xset <- xcms::group(xset)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------
pa <- frag4feature(pa, xset)

## ------------------------------------------------------------------------
print(head(pa@grped_df[1:3]))

## ------------------------------------------------------------------------
print(pa@grped_ms2[[1]])  # fragmentation associated with the first XCMS grouped feature (i.e. xset@groups[1,])

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------
pa <- filterFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------
pa <- averageAllFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------
pa <- averageIntraFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------
pa <- averageInterFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------
td <- tempdir()
createMSP(pa, msp_file_pth = file.path(td, 'out.msp'))

## ----results='hide', message=FALSE, warning=FALSE,  echo = T-------------
q_dbPth <- createDatabase(pa, xset, outDir = td, dbName = 'test-mspurity-vignette.sqlite')

## ------------------------------------------------------------------------
result <- spectralMatching(q_dbPth, q_xcmsGroups = c(12, 27), cores=1, l_accessions=c('CCMSLIB00000577898','CE000616'))

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")

#Run xcms
#xset <- xcmsSet(msPths)
#xset <- group(xset)

# Or load an XCMS xcmsSet object saved earlier
xset <- readRDS(system.file("extdata", "tests", "xcms", "ms_only_xset.rds", package="msPurity"))
# Make sure the file paths are correct
xset@filepaths[1] <- msPths[basename(msPths)=="LCMS_1.mzML"]
xset@filepaths[2] <- msPths[basename(msPths)=="LCMS_2.mzML"]

## ------------------------------------------------------------------------
px <- purityX(xset, cores = 1, xgroups = c(1, 2), ilim=0)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE)
ppDIMS <- purityD(inDF, mzML=TRUE)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
ppDIMS <- averageSpectra(ppDIMS, snMeth = "median", snthr = 5)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
ppDIMS <- filterp(ppDIMS, thr=5000, rsd = 10)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
ppDIMS <- subtract(ppDIMS)

## ------------------------------------------------------------------------
ppDIMS <- dimsPredictPurity(ppDIMS)

print(head(ppDIMS@avPeaks$processed$B02_Daph_TEST_pos))

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE----------
mzpth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
predicted <- dimsPredictPuritySingle(filepth = mzpth, mztargets = c(111.0436, 113.1069))
print(predicted)

