##general workflow from https://ccms-ucsd.github.io/GNPSDocumentation/featurebasedmolecularnetworking-with-xcms3/
#Import data (readMSData)
#Peak picking (findChromPeaks)
#Retention time alignment (adjustRtime).
#Peak grouping (groupChromPeaks).
#Gap filling (fillChromPeaks).
#Run CAMERA for adduct annotation (xsAnnotate).

load('C:/Users/m_r_j/Desktop/Rdata_testing_msPurityupdates.Rdata')


library(msPurity)
library(xcms)

#read in files and data
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE,
                       pattern = "MSMS")
ms_data = readMSData(msmsPths, mode = 'onDisk', msLevel. = 1)

#perform feature detection in individual files
cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10,
                     peakwidth = c(3, 30))
xdata <- findChromPeaks(ms_data, param = cwp)

#perform retention time correction
xdata_adj <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))

#group features across samples
sg = rep(1, length(xdata_adj$sampleNames))
pdp <- PeakDensityParam(sampleGroups = sg, minFraction = 0, bw = 30)
xdata_adj <- groupChromPeaks(xdata_adj, param = pdp)
xdata_adj_store = xdata_adj

##################
# msPurity section
##################

#create purity A object
pa = purityA(msmsPths)

######
obj = xdata_adj
peaklist = featureDefinitions(xdata_adj)

xset <- as(filterMsLevel(xdata_adj, msLevel = 1L), "xcmsSet")
xsa <- xsAnnotate(xset, polarity = "positive")
xsaF <- groupFWHM(xsa, sigma = 6, perfwhm = 1)
xsaC <- groupCorr(xsaF, cor_eic_th = 0.6, pval = 0.05, graphMethod = "hcs",
                  calcCiS = TRUE, calcCaS = TRUE, calcIso = FALSE)
xsaFI <- findIsotopes(xsaC, maxcharge = 2, maxiso = 3, minfrac = 0.5,
                      ppm = 10, intval = "maxo")
xsaFA <- findAdducts(xsaFI, polarity = "positive",
                     max_peaks = 100, multiplier = 3, ppm = 10)

