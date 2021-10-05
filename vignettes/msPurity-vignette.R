## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
library(msPurity)
msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)

## -----------------------------------------------------------------------------
pa <- purityA(msPths)

print(pa@puritydf[1:3,])

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
pa_norm <- purityA(msPths[3], iwNorm=TRUE, iwNormFun=iwNormGauss(sdlim=3, minOff=-0.5, maxOff=0.5))

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------

##for xcms version 3+
suppressPackageStartupMessages(library(xcms))
suppressPackageStartupMessages(library(MSnbase))
suppressPackageStartupMessages(library(magrittr))

#read in data and subset to use data between 30 and 90 seconds and 100 and 200 m/z
msdata = MSnbase::readMSData(msPths, mode = 'onDisk', msLevel. = 1)
rtr = c(30, 90)
mzr = c(100, 200)
msdata = msdata %>%  MSnbase::filterRt(rt = rtr) %>%  MSnbase::filterMz(mz = mzr)

#perform feature detection in individual files
cwp <- CentWaveParam(snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30))
xcmsObj <- xcms::findChromPeaks(msdata, param = cwp)
#update metadata
#for(i in 1:length(msPths)){
#  xcmsObj@processingData@files[i] <- msPths[i]
#}

xcmsObj@phenoData@data$class = c('blank', 'blank', 'sample', 'sample')
xcmsObj@phenoData@varMetadata = data.frame('labelDescription' = c('sampleNames', 'class'))
#group chromatographic peaks across samples (correspondence analysis)
pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = 0, bw = 5, binSize = 0.017)
xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------
pa <- frag4feature(pa = pa, xcmsObj = xcmsObj)

## -----------------------------------------------------------------------------
print(pa@grped_df[c(48,49),])

## -----------------------------------------------------------------------------
print(pa@grped_ms2[[18]])  # fragmentation associated with the first XCMS grouped feature (i.e. xcmsObj@groups[432,] for xcms versions < 3 and featureDefinitions(xcmsObj)[432,] for xcms v3+)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------
pa <- filterFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------
pa <- averageAllFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------
pa <- averageIntraFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------
pa <- averageInterFragSpectra(pa)

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------
td <- tempdir()
createMSP(pa, msp_file_pth = file.path(td, 'out.msp'))

## ----results='hide', message=FALSE, warning=FALSE,  echo = T------------------
q_dbPth <- createDatabase(pa = pa, xcmsObj = xcmsObj, outDir = td, dbName = 'test-mspurity-vignette.sqlite')

## -----------------------------------------------------------------------------
result <- spectralMatching(q_dbPth, q_xcmsGroups = c(432), cores=1, l_accessions=c('CCMSLIB00003740033'))

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")


## run xcms (version 3+)
# suppressPackageStartupMessages(library(xcms))
# suppressPackageStartupMessages(library(MSnbase))
# suppressPackageStartupMessages(library(magrittr))
#
# #read in data and subset to use data between 30 and 90 seconds and 100 and 200 m/z
# msdata = readMSData(msPths, mode = 'onDisk', msLevel. = 1)
# rtr = c(30, 90)
# mzr = c(100, 200)
# msdata = msdata %>%  MSnbase::filterRt(rt = rtr) %>%  MSnbase::filterMz(mz = mzr)
#
# #perform feature detection in individual files
# cwp <- CentWaveParam(snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30))
# xcmsObj <- xcms::findChromPeaks(msdata, param = cwp)
# #update metadata
# for(i in 1:length(msPths)){
#   xcmsObj@processingData@files[i] <- msPths[i]
# }
#
# xcmsObj@phenoData@data$class = c('sample', 'sample')
# xcmsObj@phenoData@varMetadata = data.frame('labelDescription' = c('sampleNames', 'class'))
# #group chromatographic peaks across samples (correspondence analysis)
# pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = 0, bw = 5, binSize = 0.017)
# xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)

## Or load an XCMS xcmsSet object saved earlier
xcmsObj <- readRDS(system.file("extdata", "tests", "xcms", "ms_only_xcmsnexp.rds", package="msPurity"))

## Make sure the file paths are correct
xcmsObj@processingData@files[1] = msPths[basename(msPths)=="LCMS_1.mzML"]
xcmsObj@processingData@files[2] = msPths[basename(msPths)=="LCMS_2.mzML"]

## -----------------------------------------------------------------------------
px <- purityX(xset =  as(xcmsObj, 'xcmsSet'), cores = 1, xgroups = c(1, 2), ilim=0)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE)
ppDIMS <- purityD(inDF, mzML=TRUE)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
ppDIMS <- averageSpectra(ppDIMS, snMeth = "median", snthr = 5)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
ppDIMS <- filterp(ppDIMS, thr=5000, rsd = 10)

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
detach("package:magrittr", unload=TRUE)
ppDIMS <- subtract(ppDIMS)

## -----------------------------------------------------------------------------
ppDIMS <- dimsPredictPurity(ppDIMS)

print(head(ppDIMS@avPeaks$processed$B02_Daph_TEST_pos))

## ----results='hide', message=FALSE, warning=FALSE,  echo = TRUE---------------
mzpth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
predicted <- dimsPredictPuritySingle(filepth = mzpth, mztargets = c(111.0436, 113.1069))
print(predicted)

