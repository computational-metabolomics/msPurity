pkgname <- "msPurity"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('msPurity')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Getfiles")
### * Getfiles

flush(stderr()); flush(stdout())

### Name: Getfiles
### Title: Get files for DI-MS processing
### Aliases: Getfiles

### ** Examples


datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)



cleanEx()
nameEx("assessPuritySingle")
### * assessPuritySingle

flush(stderr()); flush(stdout())

### Name: assessPuritySingle
### Title: Assess the purity of a single LC-MS/MS or DI-MS/MS file
### Aliases: assessPuritySingle

### ** Examples

filepth <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")

puritydf <- assessPuritySingle(filepth)



cleanEx()
nameEx("averageSpectra-purityD-method")
### * averageSpectra-purityD-method

flush(stderr()); flush(stdout())

### Name: averageSpectra,purityD-method
### Title: Using purityD object, calculates to average mz, intensity and
###   signal-to-noise of multiple scans from multiple MS datafiles (mzML or
###   .csv)
### Aliases: averageSpectra,purityD-method averageSpectra

### ** Examples


datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
ppDIMS <- averageSpectra(ppDIMS)



cleanEx()
nameEx("averageSpectraSingle")
### * averageSpectraSingle

flush(stderr()); flush(stdout())

### Name: averageSpectraSingle
### Title: Calculates to average mz, intensity and signal-to-noise of
###   multiple scans from 1 MS datafile (mzML or .csv)
### Aliases: averageSpectraSingle

### ** Examples

mzmlPth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
avP <- averageSpectraSingle(mzmlPth)



cleanEx()
nameEx("dimsPredictPurity-purityD-method")
### * dimsPredictPurity-purityD-method

flush(stderr()); flush(stdout())

### Name: dimsPredictPurity,purityD-method
### Title: Using purityD object, assess anticipated purity from a DI-MS run
### Aliases: dimsPredictPurity,purityD-method dimsPredictPurity

### ** Examples


datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
ppDIMS <- averageSpectra(ppDIMS)
ppDIMS <- filterp(ppDIMS)
ppDIMS <- subtract(ppDIMS)
ppDIMS <- dimsPredictPurity(ppDIMS)



cleanEx()
nameEx("dimsPredictPuritySingle")
### * dimsPredictPuritySingle

flush(stderr()); flush(stdout())

### Name: dimsPredictPuritySingle
### Title: Predict the precursor purity from a DI-MS dataset
### Aliases: dimsPredictPuritySingle

### ** Examples

mzmlPth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
predicted <- dimsPredictPuritySingle(c(173.0806, 216.1045), filepth=mzmlPth , minOffset=0.5, maxOffset=0.5, ppm=5, mzML=TRUE)



cleanEx()
nameEx("filterp-purityD-method")
### * filterp-purityD-method

flush(stderr()); flush(stdout())

### Name: filterp,purityD-method
### Title: Filter out peaks based on intensity and RSD criteria
### Aliases: filterp,purityD-method filterp

### ** Examples


datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)

ppDIMS <- purityD(inDF, cores=1)
ppDIMS <- averageSpectra(ppDIMS)
ppDIMS <- filterp(ppDIMS, thr = 5000)



cleanEx()
nameEx("frag4feature-purityA-method")
### * frag4feature-purityA-method

flush(stderr()); flush(stdout())

### Name: frag4feature,purityA-method
### Title: Assign precursor purity scored fragmentation spectra to XCMS
###   features
### Aliases: frag4feature,purityA-method frag4feature

### ** Examples


msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)

pa  <- purityA(msmsPths, interpol = "linear")
pa <- frag4feature(pa, xset)




cleanEx()
nameEx("getP-purityD-method")
### * getP-purityD-method

flush(stderr()); flush(stdout())

### Name: getP,purityD-method
### Title: Get peaklist for a purityD object
### Aliases: getP,purityD-method getP

### ** Examples

datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
peaks <- getP(ppDIMS)



cleanEx()
nameEx("get_additional_mzml_meta")
### * get_additional_mzml_meta

flush(stderr()); flush(stdout())

### Name: get_additional_mzml_meta
### Title: Get additional mzML meta
### Aliases: get_additional_mzml_meta

### ** Examples

mzml_pth <- system.file("extdata", "dims", "mzML", 'B02_Daph_TEST_pos.mzML', package="msPurityData")
meta_df <- get_additional_mzml_meta(mzml_pth)



cleanEx()
nameEx("groupPeaks-purityD-method")
### * groupPeaks-purityD-method

flush(stderr()); flush(stdout())

### Name: groupPeaks,purityD-method
### Title: Using purityD object, group multiple peaklists by similar mz
###   values (mzML or .csv)
### Aliases: groupPeaks,purityD-method groupPeaks

### ** Examples


datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
ppDIMS <- averageSpectra(ppDIMS)
grpedP <- groupPeaks(ppDIMS)



cleanEx()
nameEx("groupPeaksEx")
### * groupPeaksEx

flush(stderr()); flush(stdout())

### Name: groupPeaksEx
### Title: Group peaklists from a list of dataframes
### Aliases: groupPeaksEx

### ** Examples


datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
ppDIMS <- averageSpectra(ppDIMS)
grpedP <- groupPeaks(ppDIMS)



cleanEx()
nameEx("initialize-purityD-method")
### * initialize-purityD-method

flush(stderr()); flush(stdout())

### Name: initialize,purityD-method
### Title: Constructor for S4 class to represent a DI-MS purityD
### Aliases: initialize,purityD-method

### ** Examples

datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)



cleanEx()
nameEx("iwNormGauss")
### * iwNormGauss

flush(stderr()); flush(stdout())

### Name: iwNormGauss
### Title: Gaussian normalisation for isolation window efficiency
### Aliases: iwNormGauss

### ** Examples


iwNormFun <- iwNormGauss(minOff=-0.5, maxOff=0.5)
pm <- data.frame(mz=c(99.5, 99.9, 100, 100.1, 100.5),i=c(1000, 1000, 1000, 1000, 1000))
mzmax = 100.5
mzmin = 99.5
middle <- mzmax-(mzmax-mzmin)/2
adjustmz = pm$mz-middle

# normalise the intensities
pm$normi = pm$i*iwNormFun(adjustmz)





cleanEx()
nameEx("iwNormQE.5")
### * iwNormQE.5

flush(stderr()); flush(stdout())

### Name: iwNormQE.5
### Title: Q-Exactive +/- 0.5 range, normalisation for isolation window
###   efficiency
### Aliases: iwNormQE.5

### ** Examples

iwNormFun <- iwNormQE.5()
pm <- data.frame(mz=c(99.5, 99.9, 100, 100.1, 100.5),i=c(1000, 1000, 1000, 1000, 1000))
mzmax = 100.5
mzmin = 99.5
middle <- mzmax-(mzmax-mzmin)/2
adjustmz = pm$mz-middle

# normalise the intensities
pm$normi = pm$i*iwNormFun(adjustmz)




cleanEx()
nameEx("iwNormRcosine")
### * iwNormRcosine

flush(stderr()); flush(stdout())

### Name: iwNormRcosine
### Title: Raised cosine normalisation for isolation window efficiency
### Aliases: iwNormRcosine

### ** Examples

iwNormFun <- iwNormRcosine()
pm <- data.frame(mz=c(99.5, 99.9, 100, 100.1, 100.5),i=c(1000, 1000, 1000, 1000, 1000))
mzmax = 100.5
mzmin = 99.5
middle <- mzmax-(mzmax-mzmin)/2
adjustmz = pm$mz-middle

# normalise the intensities
pm$normi = pm$i*iwNormFun(adjustmz)



cleanEx()
nameEx("pcalc")
### * pcalc

flush(stderr()); flush(stdout())

### Name: pcalc
### Title: Perform purity calculation on a peak matrix
### Aliases: pcalc

### ** Examples

pm <- rbind(c(100, 1000),c(101.003, 10))
pcalc(pm, mzmin = 98, mzmax = 102, mztarget=100, ppm=5)
pcalc(pm, mzmin = 98, mzmax = 102, mztarget=100, ppm=5, isotopes = TRUE)




cleanEx()
nameEx("purityA")
### * purityA

flush(stderr()); flush(stdout())

### Name: purityA
### Title: Assess the purity of multiple LC-MS/MS or DI-MS/MS files
###   (constructor)
### Aliases: purityA

### ** Examples

filepths <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")
pa <- purityA(filepths)



cleanEx()
nameEx("purityD-class")
### * purityD-class

flush(stderr()); flush(stdout())

### Name: purityD-class
### Title: An S4 class to represent a DI-MS purityD
### Aliases: purityD-class purityD

### ** Examples

datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)



cleanEx()
nameEx("purityX")
### * purityX

flush(stderr()); flush(stdout())

### Name: purityX
### Title: Assessing anticipated purity of XCMS features from an LC-MS run
### Aliases: purityX

### ** Examples

msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
xset <- xcms::xcmsSet(msPths)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)
ppLCMS <- purityX(xset, cores = 1, xgroups = c(1, 2))




cleanEx()
nameEx("spectral_matching")
### * spectral_matching

flush(stderr()); flush(stdout())

### Name: spectral_matching
### Title: Spectral matching
### Aliases: spectral_matching

### ** Examples

msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)

pa  <- purityA(msmsPths)
pa <- frag4feature(pa, xset)
#NOTE that scan_ids here are refer the unique scan id calculated by purityA (pids).
#Only required if you want to limit the spectral matching to certain scans
result <- spectral_matching(pa@db_path, scan_ids = c(1120,  366, 1190, 601,  404,1281, 1323, 1289))



cleanEx()
nameEx("subtract-purityD-method")
### * subtract-purityD-method

flush(stderr()); flush(stdout())

### Name: subtract,purityD-method
### Title: Using Subtract MZ values based on ppm tolerance and noise ratio
### Aliases: subtract,purityD-method subtract

### ** Examples

datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)

ppDIMS <- purityD(inDF, cores=1)
ppDIMS <- averageSpectra(ppDIMS)
ppDIMS <- filterp(ppDIMS, thr = 5000)
ppDIMS <- subtract(ppDIMS)



cleanEx()
nameEx("subtractMZ")
### * subtractMZ

flush(stderr()); flush(stdout())

### Name: subtractMZ
### Title: Subtract MZ values based on ppm tolerance and noise ratio
### Aliases: subtractMZ

### ** Examples

mz1 <- c(100.001, 200.002, 300.302)
mz2 <- c(100.004, 200.003, 500.101)
i1 <- c(100, 100, 100)
i2 <- c(100, 10000, 100)

subtractMZ(mz1, mz2, i1, i2, ppm=5, s2bthres =10)



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
