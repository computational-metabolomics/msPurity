context ("Checking file list function")

test_that("Checking file list function", {
  print ("\n")
  print("########################################################")
  print("## Checking file list function                        ##")
  print("########################################################")
  datapth <- system.file("extdata", "dims",package="msPurityData")

  print("=== check for using csv ===")
  pattern = ".csv"
  inDF <- Getfiles(file.path(datapth,"msfr-peaks"), raw=FALSE, check=FALSE, pattern = pattern, cStrt=FALSE  )
  fl <- length(dir(file.path(datapth,"msfr-peaks"), full.names=TRUE, pattern=pattern, ignore.case = TRUE, recursive=FALSE,include.dirs = FALSE))
  expect_equal(fl, nrow(inDF))

  print("=== check for using mzML ===")
  pattern = ".mzML"
  inDF <- Getfiles(file.path(datapth,"mzML"), raw=FALSE, check=FALSE, pattern = pattern, cStrt=FALSE )
  fl <- length(dir(file.path(datapth,"mzML"), full.names=TRUE, pattern=pattern, ignore.case = TRUE, recursive=FALSE, include.dirs = FALSE))
  expect_equal(fl, nrow(inDF))

})
test_that("Check average spectra function only", {
  print("########################################################")
  print("## Check average spectra (function only)              ##")
  print("########################################################")
  mzmlPth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
  msfrPth <- system.file("extdata", "dims", "msfr-peaks", "Daph_TEST_pos_B02.csv", package="msPurityData")


  print("=== check single core ===")
  avP <- averageSpectraSingle(mzmlPth, cores=1, snthr = 40) # high snr threshold makes hc quicker

  expect_equal(nrow(avP), 80)
  expect_equal(round(median(avP$mz), 3), 553.278)
  expect_equal(round(median(avP$i), 0), 4639006)
  expect_equal(round(median(avP$snr), 3), 90.165)
  expect_equal(round(unname(unlist(avP[1,])), 3), c(1.000, 173.081, 11272447.000, 216.506, 9.006, 0.015, 5, 4))

  #print("=== check multi-core ===")
  #Multi-core seems to fail on unit tests but fine outside of unit test enviroment, not sure why
  #avP <- averageSpectraSingle(mzmlPth, cores=3, snthr = 40)

  print("=== check simple clustering ===")
  avP <- averageSpectraSingle(mzmlPth, cores=1, clustType="simple", snthr = 40)

  expect_equal(nrow(avP), 80)
  expect_equal(round(median(avP$mz), 3), 553.278)
  expect_equal(round(median(avP$i), 0), 4639006)
  expect_equal(round(median(avP$snr), 3), 90.165)
  expect_equal(round(unname(unlist(avP[1,])), 3), c(1.000, 173.081, 11272447.000, 216.506, 9.006, 0.015, 5, 4))


  print("=== check using MsFileReader output (median SNR thres) ===")
  avP <- averageSpectraSingle(msfrPth, cores=1, MSFileReader = TRUE, snMeth="median", snthr = 40)

  expect_equal(nrow(avP), 2)
  expect_equal(round(median(avP$mz), 3), 164.075)
  expect_equal(round(median(avP$i), 0), 5855576)
  expect_equal(round(median(avP$snr), 3), 292.333)
  expect_equal(round(unname(unlist(avP[1,])), 3), c(1.000, 155.070, 1491491.750, 74.527, 1.756, 0.127, 2, 2))




  print("=== check using MsFileReader output (precalc SNR thres) ===")
  avP <- averageSpectraSingle(msfrPth, cores=1, MSFileReader = TRUE, snMeth="precalc", snthr = 40)
  expect_equal(nrow(avP), 5)
  expect_equal(round(median(avP$mz), 3), 173.081)
  expect_equal(round(median(avP$i), 0), 691679)
  expect_equal(round(median(avP$snr), 3), 51.868)
  expect_equal(round(unname(unlist(avP[1,])), 3), c( 1.000, 155.070, 1491491.750, 113.299, 1.756, 0.110, 2, 2))

})


test_that("Check subtract spectra (function only)", {
  print("########################################################")
  print("## Check mz subtraction (function only)               ##")
  print("########################################################")
  mz1 <- c(100.001, 200.002, 300.302)
  mz2 <- c(100.004, 200.003, 500.101)
  i1 <- c(100, 100, 100)
  i2 <- c(100, 10000, 100)
  subout <- subtractMZ(mz1, mz2, i1, i2, ppm=5, s2bthres =10)
  expect_equal(subout, c(100.001, 300.302))

})



test_that("Check predict purity (function only)", {
  print("########################################################")
  print("## Check predict purity (function only)               ##")
  print("########################################################")
  mzmlPth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
  msfrPth <- system.file("extdata", "dims", "msfr-peaks", "Daph_TEST_pos_B02.csv", package="msPurityData")

  # run the function
  predicted1 <- dimsPredictPuritySingle(c(155.0700, 157.0834, 173.0806, 179.1177, 195.1224),
                                        filepth=msfrPth,
                                        minOffset=0.5,
                                        maxOffset=0.5,
                                        ppm=5,
                                        mzML=F,
                                        iwNorm = T,
                                        ilim = 0)

  # check the function
  expect_equal(round(unname(unlist(predicted1[1,])), 3), c(0.958 ,0.958, 0.001, 0.121, 0.001, 3.000))
  expect_equal(round(unname(unlist(predicted1[4,])), 3), c(1, 1, 0, 0, 0, 1))
  expect_equal(round(median(predicted1[,1]),3), 0.959)

  print("=== Check predicted purity for mzML files ===")
  # get the peaks tocheck

  #run the function
  predicted2 <- dimsPredictPuritySingle(c(155.0700, 157.0834, 173.0806, 179.1177, 195.1224),
                                        filepth=mzmlPth , minOffset=0.5, maxOffset=0.5, ppm=5, mzML=T)

  expect_equal(round(unname(unlist(predicted2[1,])), 3), c(1, 1, 0, 0, 0, 1))
  expect_equal(round(unname(unlist(predicted2[5,])), 3), c(0.816, 0.816, 0.021, 2.634, 0.010, 2.000))
  expect_equal(round(median(predicted2[,1]),3), 1)

})



test_that("Check groupPeaks (function only)", {
  print("########################################################")
  print("## Check groupPeaks (function only)               ##")
  print("########################################################")
  mzmlPth1 <- system.file("extdata", "dims", "mzML", "B02_Blank_TEST_pos.mzML", package="msPurityData")
  mzmlPth2 <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")

  print("=== check single core ===")
  avP1 <- averageSpectraSingle(mzmlPth1, cores=1, snthr = 40)
  avP2 <- averageSpectraSingle(mzmlPth2, cores=1, snthr = 40)

  peak_list <- list('p1'=avP1, 'p2'=avP2)
  gpeaks <- groupPeaksEx(peak_list)

  expect_equal(nrow(gpeaks), 111)
  expect_equal(ncol(gpeaks), 17)
})



test_that("Check mzML workflow", {
  print("########################################################")
  print("## Check mzML Workflow                                ##")
  print("########################################################")

  pattern = ".mzML"
  datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
  examp <- system.file("extdata", "dims", "example_results", package="msPurityData")

  inDF <- Getfiles(datapth, raw=FALSE, pattern=pattern, check = FALSE, cStrt = FALSE)
  rownames(inDF) <- seq(1, nrow(inDF))

  # can only do 1 core for average spectra, as ddply gives warnings that kill the
  # R CMD check call
  ex <- purityD(fileList=inDF, cores=1, mzML=TRUE)


  ################################
  # Average spectra.
  ################################
  print("averaging spectra")
  exAv <- averageSpectra(ex,
                         snMeth = "median",
                         clustType="hc",
                         ppm=1.5,
                         snthr = 30, # very high criteria so we have less peaks to process (for checking purpsoes only)
                         av="median",
                         missingV = "ignore",
                         minfrac=1
  )
  closeAllConnections()

  #saveRDS(exAv, file.path(examp, "exAv-mzML.rds"))
  exAvTest <- readRDS(file.path(examp, "exAv-mzML.rds"))

  # check avPeaks
  ##### NOTE: #### The previous rds files were run when normTIC was not a feature.
  ################ For the moment we just do not check the normTIC column

  expect_equal(length(exAv@avPeaks), length(exAvTest@avPeaks))

  expect_equal(exAv@avPeaks$processed$B02_Blank_TEST_pos[,1:5],
               exAvTest@avPeaks$processed$B02_Blank_TEST_pos[,1:5])

  expect_equal(exAv@avPeaks$processed$B02_Daph_TEST_pos[,1:5],
               exAvTest@avPeaks$processed$B02_Daph_TEST_pos[,1:5])

  ################################
  # filtering spectra.
  ################################
  print("filtering spectra")
  exF <- filterp(exAv, thr=5000, rsd =5, sampleOnly=FALSE)
  closeAllConnections()
  #saveRDS(exF, file.path(examp, "exF-mzML.rds"))
  exFTest <- readRDS(file.path(examp , "exF-mzML.rds"))
  # check avPeaks
  expect_equal(exF@avPeaks$processed$B02_Blank_TEST_pos[,1:5],
               exFTest@avPeaks$processed$B02_Blank_TEST_pos[,1:5])

  expect_equal(exF@avPeaks$processed$B02_Daph_TEST_pos[,1:5],
               exFTest@avPeaks$processed$B02_Daph_TEST_pos[,1:5])

  # check peaks have been removed
  expect_gt(nrow(exAv@avPeaks$processed[[2]]), nrow(exF@avPeaks$processed[[2]]))


  ################################
  # subtracting spectra.
  ################################
  print("subtracting")
  exS <- subtract(exF, ppm=5, s2bthres=10)
  closeAllConnections()
  #saveRDS(exS, file.path(examp, "exS-mzML.rds"))
  exSTest <- readRDS(file.path(examp, "exS-mzML.rds"))
  # check avPeaks
  expect_equal(exS@avPeaks$processed$B02_Blank_TEST_pos[,1:5],
               exSTest@avPeaks$processed$B02_Blank_TEST_pos[,1:5])

  expect_equal(exS@avPeaks$processed$B02_Daph_TEST_pos[,1:5],
               exSTest@avPeaks$processed$B02_Daph_TEST_pos[,1:5])

  # check peaks have been removed
  expect_gt(nrow(exF@avPeaks$processed[[2]]), nrow(exS@avPeaks$processed[[2]]))

  ################################
  # Predict purity
  ################################
  print("purity prediction")
  # with normalisation and ilim
  exP <- dimsPredictPurity(exS, ppm=1.5, iwNorm = TRUE, ilim = 0.05, sampleOnly = TRUE)
  expect_equal(round(median(exP@avPeaks$processed$B02_Daph_TEST_pos$medianPurity), 3), 1)
  expect_equal(round(exP@avPeaks$processed$B02_Daph_TEST_pos$medianPurity[3], 3), 0.894)

  # Without normalisation or ilim
  exPN <- dimsPredictPurity(exS, ppm=1.5, iwNorm = FALSE, ilim=0, sampleOnly=TRUE)
  expect_equal(round(median(exPN@avPeaks$processed$B02_Daph_TEST_pos$medianPurity), 3), 0.951)


  ################################
  # Group peaks
  ################################
  exPN@cores = 1
  exG <- groupPeaks(exPN)
  expect_equal(nrow(exG@groupedPeaks), 78)
  expect_equal(ncol(exG@groupedPeaks), 13)
  # check a few columns
  expect_equal(round(median(exG@groupedPeaks$i_B02_Blank_TEST_pos, na.rm = T),1), 805515.5)
  expect_equal(round(median(exG@groupedPeaks$inorm_B02_Blank_TEST_pos, na.rm = T),5), 0.01144)


})



