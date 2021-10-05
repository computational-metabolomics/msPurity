context ("checking purityA")

test_that("checking purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking purityA                                   ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")

  pa  <- purityA(msmsPths)

  #saveRDS(pa, file.path("inst", "extdata", "tests", "purityA", "1_purityA_pa.rds"))

  expect_equal(nrow(pa@puritydf), 1658)
  expect_equal(round(pa@puritydf$inPurity[[1]],4), 1)
  expect_equal(round(median(pa@puritydf$inPurity),2), 0.85)
  expect_equal(round(pa@puritydf$aMz[[1]],4),  391.2838)

  pa_saved <- readRDS(system.file("extdata", "tests", "purityA", "1_purityA_pa.rds", package="msPurity"))
  expect_equal(pa@puritydf, pa_saved@puritydf)

})


test_that("checking frag4feature", {
  print("\n")
  print("########################################################")
  print("## Checking frag4feature                              ##")
  print("########################################################")

  #read in files and data
  #msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  #ms_data = readMSData(msmsPths, mode = 'onDisk', msLevel. = 1)

  #rtr = c(30, 90)
  #mzr = c(100, 200)

  #ms_data = ms_data %>%
  #  filterRt(rt = rtr) %>%
  #  filterMz(mz = mzr)

  #perform feature detection in individual files
  #cwp <- CentWaveParam(snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30))
  #xcmsObj <- xcms::findChromPeaks(ms_data, param = cwp)

  #perform retention time correction
  #xcmsObj <- adjustRtime(xcmsObj, param = ObiwarpParam(binSize = 0.6))

  #group features across samples
  #sg = rep(1, length(xcmsObj$sampleNames))
  #pdp <- PeakDensityParam(sampleGroups = sg, minFraction = 0, bw = 5, binSize = 0.017)
  #xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)

  #create purity A object
  #pa = purityA(msmsPths)

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")

  fns = c(system.file("extdata","tests", "xcms", "msms_only_xcmsnexp.rds", package="msPurity"),
          system.file("extdata","tests", "xcms", "msms_only_xset.rds", package="msPurity"))

  for(fn in fns){
    pa <- readRDS(system.file("extdata", "tests", "purityA", "1_purityA_pa.rds", package="msPurity"))
    pa@fileList[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
    pa@fileList[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
    xcmsObj <- readRDS(fn)
    if('xcmsSet' == class(xcmsObj)){
      xcmsObj@filepaths[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
      xcmsObj@filepaths[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
    }else{
      xcmsObj@processingData@files[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
      xcmsObj@processingData@files[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
    }

    pa <- frag4feature(pa = pa, xcmsObj = xcmsObj, create_db=FALSE)

    #saveRDS(pa, file.path("inst", "extdata", "test_data", "purityA", "2_frag4feature_pa.rds"))

    # LCMSMS_1, precursor = 150.0583 m/z @ ~105 seconds [top 3 MS2 fragments (high to low): 104.0532, 133.0318 and 102.0554 m/z]
    expect_equal(round(pa@grped_df[pa@grped_df$grpid==187,]$inPurity[1],4), 1)
    expect_equal(round(pa@grped_df[pa@grped_df$grpid==187,]$precurMtchPPM[1], 4), 0.1294)
    expect_equal(length(pa@grped_ms2), 32)
    expect_equal(nrow(pa@grped_ms2[[18]][[1]]), 7) #first MS2 spectrum for grpid==140, file # 1
    expect_equal(round(pa@grped_ms2[[18]][[1]][2],4), 104.0532) #base peak in second MS2 spectrum for grpid==84, file #1

    pa_saved <- readRDS(system.file("extdata", "tests", "purityA", "2_frag4feature_pa.rds", package="msPurity"))
    expect_equal(pa@grped_ms2, pa_saved@grped_ms2)
    expect_equal(pa@grped_df, pa_saved@grped_df)

  }

})


test_that("checking frag4feature (fillpeaks)", {
  print("\n")
  print("########################################################")
  print("## Checking frag4feature (fillpeaks)                  ##")
  print("########################################################")
  library(xcms)
  # library(msPurity)
  # library(xcms)
  # library(magrittr)
  #
  # #### read in files and data
  # msPths <- dirname(list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = 'MSMS'))
  # msPths[1] <- file.path(msPths[1], 'LCMSMS_1.mzML')
  # msPths[2] <- file.path(msPths[2], 'LCMSMS_2.mzML')
  # ms_data = readMSData(msPths, mode = 'onDisk', msLevel. = 1)
  #
  # ##### subset the data to focus on retention times 30-90 seconds and m/z values between 100 and 200 m/z.
  # rtr = c(30, 90)
  # mzr = c(100, 200)
  # ms_data = ms_data %>%  filterRt(rt = rtr) %>%  filterMz(mz = mzr)
  #
  # ##### perform feature detection in individual files
  # cwp <- CentWaveParam(snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30))
  # xcmsObj <- xcms::findChromPeaks(ms_data, param = cwp)
  # xcmsObj@phenoData@data$class = c('sample', 'sample')
  # xcmsObj@phenoData@varMetadata = data.frame('labelDescription' = 'sampleNames', 'class')
  # pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = 0, bw = 5, binSize = 0.017)
  # xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)
  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")

  fns = c(system.file("extdata","tests", "xcms", "msms_only_xcmsnexp.rds", package="msPurity"),
          system.file("extdata","tests", "xcms", "msms_only_xset.rds", package="msPurity"))

  for(fn in fns){
    print(fn)
    pa <- readRDS(system.file("extdata", "tests", "purityA", "1_purityA_pa.rds", package="msPurity"))
    pa@fileList[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
    pa@fileList[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
    xcmsObj <- readRDS(fn)

    if('XCMSnExp' == class(xcmsObj)[1]){
      #xcmsObj <- xcms::fillChromPeaks(xcmsObj, param = ChromPeakAreaParam())
      xcmsObj@processingData@files[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
      xcmsObj@processingData@files[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
      xcmsObj <- xcms::fillChromPeaks(xcmsObj, param = FillChromPeaksParam(expandMz = 0, expandRt = 0, ppm = 0))
    }else if('xcmsSet' == class(xcmsObj)[1]){
      xcmsObj@filepaths[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
      xcmsObj@filepaths[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
      xcmsObj <- xcms::fillPeaks.chrom(xcmsObj, expand.mz = 0, expand.rt = 0)
    }else{
      stop('obj is not of class XCMSnExp or xcmsSet')
    }

    pa <- frag4feature(pa = pa, xcmsObj = xcmsObj, createDb = FALSE)

    #saveRDS(pa, file.path("inst", "extdata", "test_data", "purityA", "2_frag4feature_pa.rds"))

    # LCMSMS_1, precursor = 150.0583 m/z @ ~105 seconds
    expect_equal(round(pa@grped_df[pa@grped_df$grpid==187,]$inPurity[2],4), 1)
    expect_equal(round(pa@grped_df[pa@grped_df$grpid==187,]$precurMtchPPM[1], 4), 0.1294)
    #expect_equal(round(pa@grped_df$precurMtchPPM[1], 4), 1.0048)
    expect_equal(length(pa@grped_ms2), 32)
    #expect_equal(length(pa@grped_ms2), 77)
    expect_equal(nrow(pa@grped_ms2[[18]][[1]]), 7) #first MS2 spectrum for grpid==140, file # 1
    #expect_equal(nrow(pa@grped_ms2[[2]][[1]]), 4)
    expect_equal(round(pa@grped_ms2[[18]][[1]][2],4), 104.0532) #base peak in first MS2 spectrum, for grpid==140, file #1
    #expect_equal(round(pa@grped_ms2[[1]][[1]][1],4), 112.0509)

    pa_saved <- readRDS(system.file("extdata", "tests", "purityA", "2_frag4feature_pa.rds", package="msPurity"))
    expect_equal(pa@grped_ms2, pa_saved@grped_ms2)
    expect_equal(pa@grped_df, pa_saved@grped_df)

  }
})



test_that("checking filterFragSpectra purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking filterFragSpectra                         ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "2_frag4feature_pa.rds", package="msPurity"))

  pa <- filterFragSpectra(pa, plim = 0.7, snr = 3)

  #saveRDS(pa, file.path("inst", "extdata", "purityA_tests", "3_filterFragSpectra_pa.rds"))
  expect_equal(colnames(pa@grped_ms2[[18]][[1]]), c("mz", "i","snr", "ra" ,
                                                   "purity_pass_flag",    "intensity_pass_flag",
                                                   "ra_pass_flag","snr_pass_flag",   "pass_flag"  ))

  ## grpid == 84, precursor mz = 150.0581, retention time ~ 63 s
  expect_equal(round(pa@grped_ms2[[18]][[1]][,'mz'],4), c(102.0554, 104.0532, 105.0009, 105.0372, 133.0318, 137.4163, 150.0583))
  expect_equal(round(pa@grped_ms2[[18]][[1]][,'i'],0), c(4631614, 16147574, 196574, 273244, 10814390, 30550, 1973326))
  expect_equal(round(pa@grped_ms2[[18]][[1]][,'snr'], 2), c(2.35, 8.18, 0.10, 0.14, 5.48, 0.02, 1.00))
  expect_equal(round(pa@grped_ms2[[18]][[1]][,'ra'], 2), c(28.68, 100.00, 1.22, 1.69, 66.97, 0.19, 12.22))
  expect_equal(pa@grped_ms2[[18]][[1]][,'purity_pass_flag'], c(1,1,1,1,1,1,1))
  expect_equal(pa@grped_ms2[[18]][[1]][,'intensity_pass_flag'], c(1,1,1,1,1,1,1))
  expect_equal(pa@grped_ms2[[18]][[1]][,'ra_pass_flag'], c(1,1,1,1,1,1,1))
  expect_equal(pa@grped_ms2[[18]][[1]][,'snr_pass_flag'], c(0, 1, 0, 0, 1, 0, 0))
  expect_equal(pa@grped_ms2[[18]][[1]][,'pass_flag'], c(0, 1, 0, 0, 1, 0, 0))

  pa_saved <- readRDS(system.file("extdata", "tests", "purityA", "3_filterFragSpectra_pa.rds", package="msPurity"))
  pa_saved@fileList <- basename(pa_saved@fileList)
  pa@fileList <- basename(pa@fileList)
  expect_equal(pa, pa_saved)

})





test_that("checking averageIntraFragSpectra (no filter) purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking averageIntraFragSpectra                   ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "2_frag4feature_pa.rds", package="msPurity"))

  pa <- averageIntraFragSpectra(pa)

  #saveRDS(pa, file.path("inst", "extdata", "test_data", "purityA", "4_averageIntraFragSpectra_no_filter_pa.rds"))
  expect_equal(length(pa@av_spectra), 32)
  #expect_equal(length(pa@av_spectra), 77)
  expect_equal(length(pa@av_spectra[[18]]$av_intra), 2)
  #expect_equal(length(pa@av_spectra$`12`$av_intra), 2)
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[1]]$mz, 4), c(102.0554, 104.0532, 105.0009, 105.0372, 133.0318, 137.4163, 150.0583))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[2]]$mz, 4), c(101.5791, 102.0554, 104.0533, 105.0010, 105.0372, 133.0319,
                                                                 150.0387, 150.0583, 155.0976))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[1]]$frac, 4), c(1, 1, 1, 1, 1, 1, 1))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[2]]$frac, 4), c(1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[1]]$i, 2), c(4631614.50, 16147574.00, 196574.16, 273244.06, 10814390.00, 30549.83, 1973325.50))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[2]]$i, 2), c(20766.21, 4154222.25, 13982832.00, 131116.56, 266003.91,
                                                                9440187.00, 92384.07, 1694299.12, 22783.27))

  pa_saved <- readRDS(system.file("extdata","tests", "purityA", "4_averageIntraFragSpectra_no_filter_pa.rds", package="msPurity"))
  pa_saved@fileList <- basename(pa_saved@fileList)
  pa@fileList <- basename(pa@fileList)
  expect_equal(pa, pa_saved)


})

test_that("checking averageInterFragSpectra (no filter) purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking averageInterFragSpectra                   ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "4_averageIntraFragSpectra_no_filter_pa.rds", package="msPurity"))

  pa <- averageInterFragSpectra(pa)

  #saveRDS(pa, file.path("inst", "extdata", "purityA_tests", "5_averageInterFragSpectra_no_filter_pa.rds"))
  expect_equal(length(pa@av_spectra[[18]]$av_inter), 15)
  expect_equal(round(pa@av_spectra[[18]]$av_inter$mz, 4), c(101.5791, 102.0554, 104.0533, 105.0010, 105.0372,
                                                            133.0319, 137.4163, 150.0387, 150.0583, 155.0976))
  expect_equal(round(pa@av_spectra[[18]]$av_inter$frac, 4), c(0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 1.0, 0.5))
  expect_equal(round(pa@av_spectra[[18]]$av_inter$i, 2), c(20766.21, 8785836.75, 30130406.00, 327690.72, 539247.97, 20254577.00,
                                                           30549.83, 92384.07, 3667624.62, 22783.27))

  pa_saved <- readRDS(system.file("extdata","tests", "purityA", "5_averageInterFragSpectra_no_filter_pa.rds", package="msPurity"))

  expect_equal(pa, pa_saved)

})


test_that("checking averageAllFragSpectra (no filter) purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking averageAllFragSpectra                   ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "5_averageInterFragSpectra_no_filter_pa.rds", package="msPurity"))

  pa <- averageAllFragSpectra(pa)
  #saveRDS(pa, file.path("inst", "extdata", "tests", "purityA", "6_averageAllFragSpectra_no_filter_pa.rds"))

  expect_equal(length(pa@av_spectra[[18]]$av_all), 15)
  expect_equal(round(pa@av_spectra[[18]]$av_all$mz, 4), c(101.5791, 102.0554, 104.0533, 105.0010, 105.0372,
                                                          133.0319, 137.4163, 150.0387, 150.0583, 155.0976))
  expect_equal(round(pa@av_spectra[[18]]$av_all$frac, 4), c(0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 1.0, 0.5))
  expect_equal(round(pa@av_spectra[[18]]$av_all$i, 2), c(20766.21, 8785836.75, 30130406.00, 327690.72, 539247.97, 20254577.00,
                                                         30549.83, 92384.07, 3667624.62, 22783.27))
  expect_equal(round(pa@av_spectra[[18]]$av_all$mz, 4), c(101.5791, 102.0554,104.0533, 105.0010, 105.0372, 133.0319, 137.4163,
                                                          150.0387, 150.0583, 155.0976))
  expect_equal(round(pa@av_spectra[[18]]$av_all$frac, 4), c(0.5, 1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.5, 1.0, 0.5))
  expect_equal(round(pa@av_spectra[[18]]$av_all$i, 2), c(20766.21, 8785836.75, 30130406.00, 327690.72, 539247.97, 20254577.00,
                                                         30549.83, 92384.07, 3667624.62, 22783.27))

  pa_saved <- readRDS(system.file("extdata", "tests", "purityA", "6_averageAllFragSpectra_no_filter_pa.rds", package="msPurity"))
  expect_equal(pa, pa_saved)

})






test_that("checking averageIntraFragSpectra (with filter) purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking averageIntraFragSpectra (with filter)     ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "3_filterFragSpectra_pa.rds", package="msPurity"))

  pa <- averageIntraFragSpectra(pa)

  expect_equal(length(pa@av_spectra), 32)
  expect_equal(length(pa@av_spectra[[18]]$av_intra), 2)
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[1]]$mz, 4), c(104.0532, 133.0318))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[2]]$mz, 4), c(102.0554, 104.0533, 133.0319, 150.0583))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[1]]$frac, 4), c(1, 1))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[2]]$frac, 4), c(1, 1, 1, 1))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[1]]$i, 0), c(16147574, 10814390))
  expect_equal(round(pa@av_spectra[[18]]$av_intra[[2]]$i, 0), c(4154222, 13982832, 9440187, 1694299))


  #saveRDS(pa, file.path("inst", "extdata", "tests", "purityA", "7_averageIntraFragSpectra_with_filter_pa.rds"))

  pa_saved <- readRDS(system.file("extdata","tests", "purityA", "7_averageIntraFragSpectra_with_filter_pa.rds", package="msPurity"))
  expect_equal(pa, pa_saved)


})

test_that("checking averageInterFragSpectra (with filter) purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking averageInterFragSpectra (with filter)     ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "7_averageIntraFragSpectra_with_filter_pa.rds", package="msPurity"))

  pa <- averageInterFragSpectra(pa)

  expect_equal(length(pa@av_spectra[[18]]$av_inter), 15)
  expect_equal(round(pa@av_spectra[[18]]$av_inter$mz, 4), c(102.0554, 104.0533, 133.0319, 150.0583))
  expect_equal(round(pa@av_spectra[[18]]$av_inter$frac, 4), c(0.5, 1.0, 1.0, 0.5))
  expect_equal(round(pa@av_spectra[[18]]$av_inter$i, 0), c(4154222, 30130406, 20254577, 1694299))

  #saveRDS(pa, file.path("inst", "extdata", "tests",  "purityA", "8_averageInterFragSpectra_with_filter_pa.rds"))

  pa_saved <- readRDS(system.file("extdata","tests", "purityA", "8_averageInterFragSpectra_with_filter_pa.rds", package="msPurity"))
  expect_equal(pa, pa_saved)

})


test_that("checking averageAllFragSpectra (with filter) purityA", {
  print ("\n")
  print("########################################################")
  print("## Checking averageAllFragSpectra  (with filter)    ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "8_averageInterFragSpectra_with_filter_pa.rds", package="msPurity"))

  pa <- averageAllFragSpectra(pa)

  expect_equal(length(pa@av_spectra[[18]]$av_all), 15)
  expect_equal(round(pa@av_spectra[[18]]$av_all$mz, 4), c(102.0554, 104.0533, 133.0319, 150.0583))
  expect_equal(round(pa@av_spectra[[18]]$av_all$frac, 4), c(0.5, 1.0, 1.0, 0.5))
  expect_equal(round(pa@av_spectra[[18]]$av_all$i, 0), c(4154222, 30130406, 20254577, 1694299))
  expect_equal(round(pa@av_spectra[[18]]$av_all$mz, 4), c(102.0554, 104.0533, 133.0319, 150.0583))
  expect_equal(round(pa@av_spectra[[18]]$av_all$frac, 4), c(0.5, 1.0, 1.0, 0.5))
  expect_equal(round(pa@av_spectra[[18]]$av_all$i, 0), c(4154222, 30130406, 20254577, 1694299))

  #saveRDS(pa, file.path("inst", "extdata", "tests", "purityA", "9_averageAllFragSpectra_with_filter_pa.rds"))

  pa_saved <- readRDS(system.file("extdata", "tests", "purityA", "9_averageAllFragSpectra_with_filter_pa.rds", package="msPurity"))
  expect_equal(pa, pa_saved)

})



test_that("checking createMSP based functions", {
  print ("\n")
  print("########################################################")
  print("## Checking createMSP functions                       ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "9_averageAllFragSpectra_with_filter_pa.rds", package="msPurity"))

  get_msp_str <- function(msp_pth){
    msp_str <- readChar(msp_pth, file.info(msp_pth)$size)

    msp_str <- gsub('\n','',msp_str)
    msp_str <- gsub('\r','',msp_str)
    msp_str <- gsub('msPurity version:\\d+\\.\\d+\\.\\d+','', msp_str)
    return(msp_str)
  }


  metadata <- data.frame('grpid'=c(187), 'MS$FOCUSED_ION: PRECURSOR_TYPE'=c('[M+H]+'),
                         'AC$MASS_SPECTROMETRY: ION_MODE'=c("POSITIVE"), 'CH$NAME:'=c('Methionine'),
                         check.names = FALSE, stringsAsFactors = FALSE)

  tmp_dir <- tempdir()


  ################################
  # Check all method
  ################################
  all_msp_new_pth <- file.path(tmp_dir,'all.msp')
  createMSP(pa, msp_file = all_msp_new_pth, metadata = metadata,
            method = "all", xcms_groupids = c(187))

  all_msp_new <- get_msp_str(all_msp_new_pth)
  all_msp_old <- get_msp_str(system.file("extdata", "tests", "msp", "all.msp", package="msPurity"))
  expect_equal(all_msp_new, all_msp_old)


  ################################
  # Check max method
  ################################
  max_msp_new_pth <- file.path(tmp_dir,'max.msp')
  createMSP(pa, msp_file = max_msp_new_pth, metadata = metadata, method = "max",
            xcms_groupids = c(187), filter=FALSE)

  max_msp_new <- get_msp_str(max_msp_new_pth)
  max_msp_old <- get_msp_str(system.file("extdata", "tests","msp", "max.msp", package="msPurity"))
  expect_equal(max_msp_new, max_msp_old)


  ################################
  # Check av_inter method
  ################################
  av_inter_msp_new_pth <- file.path(tmp_dir,'av_inter.msp')
  createMSP(pa, msp_file = av_inter_msp_new_pth, metadata = metadata, method = "av_inter", xcms_groupids = c(187))

  av_inter_msp_new <- get_msp_str(av_inter_msp_new_pth)
  av_inter_msp_old <- get_msp_str(system.file("extdata","tests","msp", "av_inter.msp", package="msPurity"))
  expect_equal(av_inter_msp_new, av_inter_msp_old)


  ################################
  # Check av_intra method
  ################################
  av_intra_msp_new_pth <- file.path(tmp_dir,'av_intra.msp')
  createMSP(pa, msp_file = av_intra_msp_new_pth, metadata = metadata, method = "av_intra", xcms_groupids = c(187))

  av_intra_msp_new <- get_msp_str(av_intra_msp_new_pth)
  av_intra_msp_old <- get_msp_str(system.file("extdata","tests","msp", "av_intra.msp", package="msPurity"))
  expect_equal(av_intra_msp_new, av_intra_msp_old)

  ################################
  # Check av_all method
  ################################
  av_all_msp_new_pth <- file.path(tmp_dir,'av_all.msp')
  createMSP(pa, msp_file = av_all_msp_new_pth, metadata = metadata, method = "av_all", xcms_groupids = c(187))

  av_all_msp_new <- get_msp_str(av_all_msp_new_pth)
  av_all_msp_old <- get_msp_str(system.file("extdata", "tests","msp", "av_all.msp", package="msPurity"))
  expect_equal(av_all_msp_new, av_all_msp_old)


  ################################
  # When two metadata details for single XCMS feature
  ################################
  metadatad <- metadata
  metadatad[2, ] <- metadatad[1,]
  metadatad[2,4] <- 'possible other compound'

  av_all_msp_new_dupmeta_pth <- file.path(tmp_dir,'av_all_dupmeta.msp')
  createMSP(pa, msp_file = av_all_msp_new_dupmeta_pth, metadata = metadatad, method = "av_all", xcms_groupids = c(187))

  av_all_msp_dupmeta_new <- get_msp_str(av_all_msp_new_dupmeta_pth)
  av_all_msp_dupmeta_old <- get_msp_str(system.file("extdata", "tests","msp", "av_all_dupmeta.msp", package="msPurity"))
  expect_equal(av_all_msp_dupmeta_new, av_all_msp_dupmeta_old)

  ################################
  # When the RECORD_TITLE: is defined by user
  ################################
  metadata <- data.frame('grpid'=c(187), 'MS$FOCUSED_ION: PRECURSOR_TYPE'=c('[M+H]+'),
                         'AC$MASS_SPECTROMETRY: ION_MODE'=c("POSITIVE"), 'RECORD_TITLE:'=c('Methionine'),
                         check.names = FALSE, stringsAsFactors = FALSE)
  recrdt_msp_new_pth <- file.path(tmp_dir,'recrdt.msp')
  createMSP(pa, msp_file = recrdt_msp_new_pth, metadata = metadata,
            method = "all", xcms_groupids = c(187))

  recrdt_msp_new <- get_msp_str(recrdt_msp_new_pth)
  recrdt_msp_old <- get_msp_str(system.file("extdata", "tests", "msp", "recrdt.msp", package="msPurity"))
  expect_equal(recrdt_msp_new, recrdt_msp_old)
})






