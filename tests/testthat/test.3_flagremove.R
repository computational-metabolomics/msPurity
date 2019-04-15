context ("checking flag and remove peaks")

test_that("checking flag and remove peaks", {
  print ("\n")
  print("########################################################")
  print("## checking flag and remove (lc-ms)                   ##")
  print("########################################################")
  library(msPurity)
  library(xcms)

  msPths <- dirname(list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE))

  msPths[1] <- file.path(msPths[1], 'LCMS_1.mzML')
  msPths[2] <- file.path(msPths[2], 'LCMS_2.mzML')
  msPths[3] <- file.path(msPths[3], 'LCMSMS_1.mzML')
  msPths[4] <- file.path(msPths[4], 'LCMSMS_2.mzML')

  xset <- xcms::xcmsSet(msPths)

  xset@phenoData[,1] <- c('blank', 'blank', 'sample', 'sample')

  #print(xset@phenoData)
  xset <- xcms::group(xset)
  fr = msPurity::flag_remove(xset)

  #print(head(fr$grp_peaklist))
  expect_equal(nrow(fr$grp_peaklist), 14) # passing on testthat but failng on check (not sure why)

  # testthat::expect_equal(nrow(fr[[2]]), 14)


})

