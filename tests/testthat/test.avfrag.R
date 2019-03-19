context ("checking averaging functionality")

test_that("checking averaging functionality", {
  print ("\n")
  print("########################################################")
  print("## Checking averaging                                ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")


  xset <- xcms::xcmsSet(msmsPths)
  xset <- xcms::group(xset)


  pa  <- purityA(msmsPths)


  pa <- frag4feature(pa, xset)
  paIntra <- averageIntraFragSpectra(pa)
  paInter <- averageInterFragSpectra(paIntra)
  paAll <- averageAllFragSpectra(paInter)

  expect_equal(ncol(paAll@av_spectra$`12`$av_intra$`1`), 15)
  expect_equal(nrow(paAll@av_spectra$`12`$av_intra$`1`), 4)
  expect_equal(round(paAll@av_spectra$`12`$av_intra$`1`$mz, 3), c(107.270, 116.017, 116.071, 116.107))
  expect_equal(round(paAll@av_spectra$`12`$av_intra$`1`$i, 3), c(1726.613, 7324.864, 2114202.398,   13797.101))
  expect_equal(ncol(paAll@av_spectra$`12`$av_inter), 15)
  expect_equal(nrow(paAll@av_spectra$`12`$av_inter), 3)

  expect_equal(round(paAll@av_spectra$`12`$av_inter$mz, 3), c(116.017, 116.071, 116.107))
  expect_equal(round(paAll@av_spectra$`12`$av_inter$i, 1), c(  31097.7, 4143476.2,   26102.3))



  pa_no_peaks <- averageIntraFragSpectra(pa, minfrac=0.1, minnum = 1, snr=10000, ppm=5, rmp = T)
  expect_equal(unlist(pa_no_peaks@av_spectra), numeric(0))

  paf <- filterFragSpectra(pa, snr=100)
  pafIntra <- averageIntraFragSpectra(paf)
  pafInter <- averageInterFragSpectra(pafIntra)
  pafAll <- averageAllFragSpectra(pafInter)

  # only 1 peak passes
  expect_equal(paf@grped_ms2$`375`[[2]][,'pass_flag'], c(0,0,1,0,0,0,0,0))
  expect_equal(paf@grped_ms2$`375`[[1]][,'pass_flag'], c(0,0,0,0))

  # only 1 should be available
  expect_equal(nrow(pafIntra@av_spectra$`375`$av_intra$`2`), 1)
  expect_equal(nrow(pafInter@av_spectra$`375`$av_inter), 1)
  expect_equal(nrow(pafAll@av_spectra$`375`$av_all), 1)


})
