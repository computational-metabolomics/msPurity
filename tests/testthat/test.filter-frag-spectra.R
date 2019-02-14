test_that("checking filterFragSpectra purityA", {
  print("########################################################")
  print("## Checking functions for filtering frag spectra      ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")


  xset <- xcms::xcmsSet(msmsPths)
  xset <- xcms::group(xset)

  pa  <- purityA(msmsPths)


  pa <- frag4feature(pa, xset)
  pa <- filterFragSpectra(pa, p = 0.7, snr = 3)

  expect_equal(colnames(pa@grped_ms2[[1]][[1]]), c("mz", "i","snr", "ra" ,
                                                   "purity_pass_flag",    "intensity_pass_flag",
                                                   "ra_pass_flag","snr_pass_flag",   "pass_flag"  ))


  expect_equal(round(pa@grped_ms2[[1]][[1]][,'mz'],4), c(112.0509, 126.5377))
  expect_equal(round(pa@grped_ms2[[1]][[1]][,'i'],0),  c(565117,   2499))
  expect_equal(pa@grped_ms2[[1]][[1]][,'snr'], c(1.991193651, 0.008806349))
  expect_equal(pa@grped_ms2[[1]][[1]][,'ra'], c(100.0000000,   0.4422648))
  expect_equal(pa@grped_ms2[[1]][[1]][,'purity_pass_flag'], c(1,   1))
  expect_equal(pa@grped_ms2[[1]][[1]][,'intensity_pass_flag'], c(1,   1))
  expect_equal(pa@grped_ms2[[1]][[1]][,'ra_pass_flag'], c(1,   1))
  expect_equal(pa@grped_ms2[[1]][[1]][,'snr_pass_flag'], c(0,   0))
  expect_equal(pa@grped_ms2[[1]][[1]][,'pass_flag'], c(0,   0))


  expect_equal(round(pa@grped_ms2[[10]][[2]][,'mz'],4), c(102.0013, 102.0376, 102.0553, 102.9030,
                                                          103.0827, 120.0116, 122.0163, 130.0501,
                                                          130.0860, 131.0527, 131.0899, 133.0192,
                                                          134.0267, 148.0426, 160.7018))
  expect_equal(round(pa@grped_ms2[[10]][[2]][,'i'],0),  c( 418361, 26940, 31788,  2297,  2656, 1534089,
                                                           2559,   17183,    4825,    3329,    3828,
                                                           31002,    3603, 1333086,    2493))
  expect_equal(pa@grped_ms2[[10]][[2]][,'snr'], c(86.7093726,   5.5836066,   6.5882984,   0.4761153,
                                                  0.5504728, 317.9551029,   0.5304692,   3.5613499,
                                                  1.0000000,   0.6899905,   0.7933997,   6.4254739,
                                                  0.7467862,276.2951708,   0.5167274))
  expect_equal(pa@grped_ms2[[10]][[2]][,'ra'], c( 27.2709486,   1.7560991,   2.0720845,   0.1497429,   0.1731291,
                                                  100.0000000,   0.1668378,   1.1200795,   0.3145098,   0.2170088,
                                                  0.2495320,   2.0208746,  0.2348716,  86.8975425,   0.1625158))
  expect_equal(pa@grped_ms2[[10]][[2]][,'purity_pass_flag'], c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(pa@grped_ms2[[10]][[2]][,'intensity_pass_flag'], c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(pa@grped_ms2[[10]][[2]][,'ra_pass_flag'], c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
  expect_equal(pa@grped_ms2[[10]][[2]][,'snr_pass_flag'], c(1, 1, 1, 0, 0, 1, 0, 1, 0, 0 ,0, 1, 0, 1, 0))
  expect_equal(pa@grped_ms2[[10]][[2]][,'pass_flag'], c(1, 1, 1, 0, 0, 1, 0, 1, 0, 0 ,0, 1, 0, 1, 0))



})
