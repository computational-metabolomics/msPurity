context ("checking flag and remove peaks")

test_that("checking flag and remove peaks", {
  print ("\n")
  print("########################################################")
  print("## checking flag and remove (lc-ms)                   ##")
  print("########################################################")
  # library(msPurity)
  # library(xcms)
  # library(magrittr)
  #read in files and data
  # msPths <- dirname(list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE))
  # msPths[1] <- file.path(msPths[1], 'LCMS_1.mzML')
  # msPths[2] <- file.path(msPths[2], 'LCMS_2.mzML')
  # msPths[3] <- file.path(msPths[3], 'LCMSMS_1.mzML')
  # msPths[4] <- file.path(msPths[4], 'LCMSMS_2.mzML')
  # ms_data = readMSData(msPths, mode = 'onDisk', msLevel. = 1)
  #
  # #subset the data to focus on retention times 30-90 seconds and m/z values between 100 and 200 m/z.
  # rtr = c(30, 90)
  # mzr = c(100, 200)
  # ms_data = ms_data %>%  filterRt(rt = rtr) %>%  filterMz(mz = mzr)
  #
  # ##### perform feature detection in individual files
  # cwp <- CentWaveParam(snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30))
  # xcmsObj <- xcms::findChromPeaks(ms_data, param = cwp)
  # xcmsObj@phenoData@data$class = c('blank', 'blank', 'sample', 'sample')
  # xcmsObj@phenoData@varMetadata = data.frame('labelDescription' = 'sampleNames', 'class')
  # pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = 0, bw = 5, binSize = 0.017)
  # xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)
  # saveRDS(xcmsObj, system.file("extdata", "tests", "purityA", "10_filterflagremove.rds", package="msPurity"))

  xcmsObj = readRDS(system.file("extdata", "tests", "purityA", "10_input_filterflagremove.rds", package="msPurity"))

  #check that this works for both xcmsObj of class 'XCMSnExp' and 'xcmsSet'
  xcmsObjs = list(xcmsObj, as(xcmsObj,'xcmsSet'))

  for(xcmsObj in xcmsObjs){
    if(class(xcmsObj) == 'xcmsSet'){
      xcmsObj@phenoData[,1] <- c('blank', 'blank', 'sample', 'sample')
    }
    fr = msPurity::flag_remove(xcmsObj, pol=NA, rsd_i_blank=NA, minfrac_blank=0.5,
                               rsd_rt_blank=NA, ithres_blank=NA, s2b=10, ref.class='blank',
                               egauss_thr=NA, rsd_i_sample=NA, minfrac_sample=0.7,
                               rsd_rt_sample=NA, ithres_sample=NA, minfrac_xcms=0.7,
                               mzwid=0.017, bw=5, out_dir='.', temp_save=FALSE, remove_spectra=TRUE, grp_rm_ids=NA, xset=NA)
    expect_equal(nrow(fr$grp_peaklist), 27) # passing on testthat but failing on check (not sure why)
  }
})
