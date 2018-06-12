test_that("checking flag and remove peaks", {
  print("########################################################")
  print("## checking flag and remove (lc-ms)                  ##")
  print("########################################################")
  library(msPurity)

  msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)
  xset <- xcms::xcmsSet(msPths)
  xset@phenoData[,1] <- c('blank', 'blank', 'sample', 'sample')
  xset <- xcms::group(xset)
  fr = flag_remove(xset)
  expect_equal(nrow(fr[[1]]@peaks), 31)


})

