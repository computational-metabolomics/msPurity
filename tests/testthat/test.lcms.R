context ("checking lcms based functions")

test_that("checking lcms based functions", {
  print ("\n")
  print("########################################################")
  print("## Checking LCMS based class and functions            ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  library(xcms)

  xset <- xcmsSet(msmsPths)
  xset <- group(xset)

  pa  <- purityA(msmsPths)
  pa <- frag4feature(pa, xset, create_db=FALSE)

  expect_equal(nrow(pa@puritydf), 1658)
  expect_equal(round(pa@puritydf$inPurity[[1]],4), 1)
  expect_equal(round(median(pa@puritydf$inPurity),2), 0.85)
  expect_equal(round(pa@puritydf$aMz[[1]],4),  391.2838)

  expect_equal(round(pa@grped_df$inPurity[1],4), 1)
  expect_equal(round(pa@grped_df$precurMtchPPM[1], 4), 1.0048)

  expect_equal(length(pa@grped_ms2), 77)
  expect_equal(nrow(pa@grped_ms2[[2]][[1]]), 4)
  expect_equal(round(pa@grped_ms2[[1]][[1]][1],4), 112.0509)

  pa <- averageIntraFragSpectra(pa)
  expect_equal(length(pa@av_spectra), 77)
  expect_equal(length(pa@av_spectra$`12`$av_intra), 2)
  expect_equal(round(pa@av_spectra$`12`$av_intra$`1`$mz, 4), c(107.2701, 116.0165, 116.0709, 116.1073))
  expect_equal(round(pa@av_spectra$`12`$av_intra$`2`$mz, 4), c(103.1290, 116.0168, 116.0709, 116.1073, 130.0276))
  expect_equal(round(pa@av_spectra$`12`$av_intra$`1`$frac, 4), c(0.3333, 0.6667, 1.0000, 0.6667))
  expect_equal(round(pa@av_spectra$`12`$av_intra$`2`$frac, 4), c(0.3333, 1.0000, 1.0000, 0.6667, 0.3333))
  expect_equal(round(pa@av_spectra$`12`$av_intra$`1`$i, 2), c(1726.61, 7324.86, 2114202.40, 13797.10))
  expect_equal(round(pa@av_spectra$`12`$av_intra$`2`$i, 2), c( 4419.71, 23772.87, 2029273.85, 12305.23, 4081.14))

  pa <- averageInterFragSpectra(pa)
  expect_equal(length(pa@av_spectra$`12`$av_inter), 15)
  expect_equal(round(pa@av_spectra$`12`$av_inter$mz, 4), c(116.0166,116.0709,116.1073))
  expect_equal(round(pa@av_spectra$`12`$av_inter$frac, 4), c(1,1,1))
  expect_equal(round(pa@av_spectra$`12`$av_inter$i, 2), c(31097.73, 4143476.25, 26102.33))

  pa <- averageAllFragSpectra(pa)
  expect_equal(length(pa@av_spectra$`12`$av_all), 15)
  expect_equal(round(pa@av_spectra$`12`$av_all$mz, 4), c(103.1290, 107.2701, 116.0166, 116.0709, 116.1073, 130.0276))
  expect_equal(round(pa@av_spectra$`12`$av_all$frac, 4), c(0.1667, 0.1667, 0.8333, 1.0000, 0.6667, 0.1667))
  expect_equal(round(pa@av_spectra$`12`$av_all$i, 2), c(4419.71, 1726.61, 31097.73, 4143476.25, 26102.33, 4081.14))
  expect_equal(round(pa@av_spectra$`12`$av_all$mz, 4), c(103.1290, 107.2701, 116.0166, 116.0709, 116.1073, 130.0276))
  expect_equal(round(pa@av_spectra$`12`$av_all$frac, 4), c(0.1667, 0.1667, 0.8333, 1.0000, 0.6667, 0.1667))
  expect_equal(round(pa@av_spectra$`12`$av_all$i, 2), c(4419.71, 1726.61, 31097.73, 4143476.25, 26102.33, 4081.14))


  msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
  xset2 <- xcmsSet(msPths)
  xset2 <- group(xset2)

  ppLCMS <- purityX(xset2, cores = 1, xgroups = c(1, 2), ilim=0)


  # msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
  # xset2 <- xcmsSet(msPths,method = "centWave")
  # xset2 <- group(xset2)
  # ppLCMS <- purityX(xset2)

  # rtraw <- xset2@peaks[,c('rt', 'rtmin', 'rtmax')]
  # colnames(rtraw) <- c('rt_raw','rtmin_raw','rtmax_raw')
  # xset2@peaks <- cbind(xset2@peaks, rtraw)
  #
  # xset2 <- group(xset2)
  #
  # xset2 <- retcor(xset2, 'obiwarp')
  # xset2 <- group(xset2)
  #
  #
  #
  # raw1 <- getXcmsRaw(xset2, sampleidx=1, rt='corrected')
  # raw2 <- getXcmsRaw(xset2, sampleidx=2, rt='corrected')
  # pl1 <- pl[pl[,'sample']==1,]
  # pl2 <- pl[pl[,'sample']==2,]


  expect_equal(round(median(ppLCMS@predictions$grpid),3), 1.5)
  expect_equal(round(median(ppLCMS@predictions$mean),3), 0.995)
  expect_equal(round(median(ppLCMS@predictions$median),3), 0.995)
  expect_equal(round(median(ppLCMS@predictions$sd, na.rm=TRUE),3), 0.001)
  expect_equal(round(median(ppLCMS@predictions$i),3), 43822022)
  expect_equal(round(median(ppLCMS@predictions$mz),3), 102.573)


  ppLCMS_sf <- purityX(xset2, singleFile = 1)

  #expect_equal(round(median(ppLCMS_sf@predictions$id),3), 514) Note that the new version of xcms
  #                                                             now gives a different number of peaks using
  #                                                             standard settings.. so we can't do this check anymore
  #                                                             as it will fail with newer builds
  expect_equal(round(median(ppLCMS_sf@predictions$purityFWHMmedian, na.rm = T),3), 0.935)
  expect_equal(round(median(ppLCMS_sf@predictions$purityFWmedian),3), 0.505)
  expect_equal(round(median(ppLCMS_sf@predictions$pknmFWHMmedian, na.rm=T ),3), 2)
  expect_equal(round(median(ppLCMS_sf@predictions$pknmFWmedian),3), 3)



})


  # For some reason the following setup does not work on the mac build
  ##############################################################
  # Perform XCMS
  ##############################################################
  # Make sure we get the file list in correct order
#   datapth <- system.file("extdata", "lcms", package="msPurityData")
#
#   files <- c("LCMS_1.mzML","LCMS_2.mzML", "LCMSMS_1.mzML", "LCMSMS_2.mzML")
#   mzdatafiles <- vector()
#   for (i in 1:length(files)){
#     mzdatafiles[i] <- file.path(datapth, "mzML", files[i])
#   }
#
#   xset <- xcms::xcmsSet(mzdatafiles,
#                    method = "centWave",
#                    nSlaves=1,
#                    ppm=5,
#                    peakwidth=c(5,20),
#                    snthresh=10000000, # just so the processing is faster
#                    prefilter=c(3,100),
#                    mzCenterFun="wMean",
#                    integrate=1,
#                    mzdiff=0.001
#   )
#   xset <- xcms::group(xset, bw=5, mzwid = 0.01, minfrac = 1)
#
#   pa2 <- purityA(xset@filepaths, cores=1,
#                  interpol = "linear", iwNorm = TRUE, ilim = 0.05)
#
#   expect_equal(nrow(pa2@puritydf), 1658)
#   expect_equal(round(pa2@puritydf$inPurity[[1]],4), 1)
#   expect_equal(round(median(pa2@puritydf$inPurity),2), 0.85)
#   expect_equal(round(pa2@puritydf$aMz[[1]],4),  391.2838)



#   # Get associated spectra for each feature from xcms
#   pa3 <- frag4feature(pa2, xset)
#
#   expect_equal(round(pa3@grped_df$inPurity[1],4), 0.7166)
#   expect_equal(round(pa3@grped_df$precurMtchPPM[1], 4), 0.1651)
#
#   expect_equal(length(pa3@grped_ms2), 8)
#   expect_equal(nrow(pa3@grped_ms2[[1]][[1]]), 14)
#   expect_equal(round(pa3@grped_ms2[[1]][[1]][1],4), 102.0916)
#
#
#   ##########################################################
#   # Precursor purity predictions LC-MS
#   ##########################################################
#
#   expect_equal(round(median(ppLCMS@predictions$grpid),3), 5.5)
#   expect_equal(round(median(ppLCMS@predictions$mean),3), 0.995)
#   expect_equal(round(median(ppLCMS@predictions$median),3), 0.995)
#   expect_equal(round(median(ppLCMS@predictions$sd, na.rm=TRUE),3), 0.001)
#   expect_equal(round(median(ppLCMS@predictions$i),3), 68391965)
#   expect_equal(round(median(ppLCMS@predictions$mz),3), 303.66)
#

