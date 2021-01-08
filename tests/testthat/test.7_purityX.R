context ("checking purityX")

test_that("checking purityX (grouped)", {
  print ("\n")
  print("########################################################")
  print("## Checking purityX (grouped)                         ##")
  print("########################################################")

  msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
  #xset <- xcmsSet(msPths)
  #xset <- group(xset)
  #saveRDS(xset, file.path("inst", "extdata", "tests", "xcms", "ms_only_xset.rds"))

  xset <- readRDS(system.file("extdata", "tests", "xcms", "ms_only_xset.rds", package="msPurity"))

  xset@filepaths[1] <- msPths[basename(msPths)=="LCMS_1.mzML"]
  xset@filepaths[2] <- msPths[basename(msPths)=="LCMS_2.mzML"]


  px <- purityX(xset, cores = 1, xgroups = c(1, 2), ilim=0)




  #saveRDS(px, file.path("inst", "extdata", "tests", "purityX", "1_purityX_grouped_px.rds"))

  expect_equal(round(median(px@predictions$grpid),3), 1.5)
  expect_equal(round(median(px@predictions$mean),3), 0.995)
  expect_equal(round(median(px@predictions$median),3), 0.995)
  expect_equal(round(median(px@predictions$sd, na.rm=TRUE),3), 0.001)
  expect_equal(round(median(px@predictions$i),3), 43822022)
  expect_equal(round(median(px@predictions$mz),3), 102.573)


  px_saved <- readRDS(system.file("extdata", "tests", "purityX", "1_purityX_grouped_px.rds", package="msPurity"))
  expect_equal(px, px_saved)

})

test_that("checking purityX (single file)", {
  print ("\n")
  print("########################################################")
  print("## Checking purityX (single file)                     ##")
  print("########################################################")

  #msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
  #xset <- xcmsSet(msPths)
  #xset <- group(xset)
  #saveRDS(xset, file.path("inst", "extdata", "tests", "xcms", "ms_only_xset.rds"))
  msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
  xset <- readRDS(system.file("extdata", "tests", "xcms", "ms_only_xset.rds", package="msPurity"))
  xset@filepaths[1] <- msPths[basename(msPths)=="LCMS_1.mzML"]
  xset@filepaths[2] <- msPths[basename(msPths)=="LCMS_2.mzML"]
  px <- purityX(xset, singleFile = 1)


  #saveRDS(px, file.path("inst", "extdata", "tests", "purityX", "1_purityX_single_px.rds"))
  expect_equal(round(median(px@predictions$purityFWHMmedian, na.rm = TRUE),3), 0.935)
  expect_equal(round(median(px@predictions$purityFWmedian),3), 0.505)
  expect_equal(round(median(px@predictions$pknmFWHMmedian, na.rm=TRUE ),3), 2)
  expect_equal(round(median(px@predictions$pknmFWmedian),3), 3)

  px_saved <- readRDS(system.file("extdata", "tests", "purityX", "1_purityX_single_px.rds", package="msPurity"))
  expect_equal(px, px_saved)
})

