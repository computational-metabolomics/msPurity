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

  xset <- readRDS(system.file("extdata", "tests", "xcms", "ms_only_xset_OLD.rds", package="msPurity"))

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
  for (name in setdiff(methods::slotNames(px), "params")) {
        expect_equal(methods::slot(px, name), methods::slot(px_saved, name),
            info = name)
    }

})

test_that("purityX() plots the EICs of each feature", {
  skip_if_no_mspuritydata()
  td <- tempfile("purityX-plots-")
  dir.create(td)
  wd <- setwd(td)
  on.exit(setwd(wd))
  px <- suppressMessages(purityX(.golden_xcms("msms_only_xset_OLD.rds"),
                                 plotP = TRUE, xgroups = c(1, 2)))
  expect_setequal(unique(px@predictions$grpid), c(1, 2))
  plots <- list.files(file.path(td, "purityXplots"), pattern = "[.]png$")
  expect_gt(length(plots), 0)
})
