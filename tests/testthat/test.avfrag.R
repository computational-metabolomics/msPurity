test_that("checking averaging functionality", {
  print("########################################################")
  print("## Checking LCMS based class and functions            ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")


  xset <- xcms::xcmsSet(msmsPths)
  xset <- xcms::group(xset)


  pa  <- purityA(msmsPths)


  pa <- frag4feature(pa, xset)
  paIntra <- averageIntraFragSpectra(pa)
  paInter <- averageInterFragSpectra(paIntra)
  paAll <- averageAllFragSpectra(paInter)

  expect_equal(ncol(paAll@av_spectra$`12`$av_intra$`1`), 14)
  expect_equal(nrow(paAll@av_spectra$`12`$av_intra$`1`), 4)
  expect_equal(round(paAll@av_spectra$`12`$av_intra$`1`$mz, 3), c(107.270, 116.017, 116.071, 116.107))
  expect_equal(round(paAll@av_spectra$`12`$av_intra$`1`$i, 3), c(1726.613, 7324.864, 2114202.398,   13797.101))
  expect_equal(ncol(paAll@av_spectra$`12`$av_inter), 14)
  expect_equal(nrow(paAll@av_spectra$`12`$av_inter), 3)

  expect_equal(round(paAll@av_spectra$`12`$av_inter$mz, 3), c(116.017, 116.071, 116.107))
  expect_equal(round(paAll@av_spectra$`12`$av_inter$i, 1), c(  31097.7, 4143476.2,   26102.3))



  pa_no_peaks <- averageIntraFragSpectra(pa, minfrac=0.5, snr_pre = 100000, minnum = 1, ppm=5, plim = 0.8, remove_peaks = T)
  expect_equal(unlist(pa_no_peaks@av_spectra), NULL)
})
