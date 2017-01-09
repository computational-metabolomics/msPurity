test_that("checking lcms based functions", {
  print("########################################################")
  print("## Checking LCMS based class and functions            ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  library(xcms)

  xset <- xcmsSet(msmsPths)
  xset <- group(xset)

  pa  <- purityA(msmsPths, interpol = "linear")
  pa <- frag4feature(pa, xset)

  expect_equal(nrow(pa@puritydf), 1658)
  expect_equal(round(pa@puritydf$inPurity[[1]],4), 1)
  expect_equal(round(median(pa@puritydf$inPurity),2), 0.85)
  expect_equal(round(pa@puritydf$aMz[[1]],4),  391.2838)

  expect_equal(round(pa@grped_df$inPurity[1],4), 1)
  expect_equal(round(pa@grped_df$precurMtchPPM[1], 4), 0.5516)

  expect_equal(length(pa@grped_ms2), 77)
  expect_equal(nrow(pa@grped_ms2[[2]][[1]]), 4)
  expect_equal(round(pa@grped_ms2[[1]][[1]][1],4), 112.0509)


  msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
  xset2 <- xcmsSet(msPths)
  xset2 <- group(xset2)
  ppLCMS <- purityX(xset2, cores = 1, xgroups = c(1, 2))

  expect_equal(round(median(ppLCMS@predictions$grpid),3), 1.5)
  expect_equal(round(median(ppLCMS@predictions$mean),3), 0.995)
  expect_equal(round(median(ppLCMS@predictions$median),3), 0.995)
  expect_equal(round(median(ppLCMS@predictions$sd, na.rm=TRUE),3), 0.001)
  expect_equal(round(median(ppLCMS@predictions$i),3), 43822022)
  expect_equal(round(median(ppLCMS@predictions$mz),3), 102.573)

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

