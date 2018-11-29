test_that("checking lcms based functions", {
  print("########################################################")
  print("## Checking LCMS based class and functions            ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")


  xset <- xcms::xcmsSet(msmsPths)
  xset <- xcms::group(xset)


  pa  <- purityA(msmsPths)


  pa <- frag4feature(pa, xset, create_db=FALSE)
  pa <- averageFragmentation(pa, minfrac_intra=0, minfrac_inter = 0, minfrac_all = 0, plim=0)

  get_msp_str <- function(msp_pth){
    msp_str <- readChar(msp_pth, file.info(msp_pth)$size)

    msp_str <- gsub('\n','',msp_str)
    msp_str <- gsub('\r','',msp_str)
    msp_str <- gsub('msPurity version:\\d+\\.\\d+\\.\\d+','', msp_str)
    return(msp_str)
  }



  metadata <- data.frame('grpid'=c(8, 12), 'MS$FOCUSED_ION: PRECURSOR_TYPE'=c('[M+H]+', '[M+H]+ 88.0158 [M+H+NH3]+ 105.042'),
                         'AC$MASS_SPECTROMETRY: ION_MODE'=c("POSITIVE","POSITIVE"), 'CH$NAME'=c('Sulfamethizole', 'Unknown'),
                         check.names = FALSE, stringsAsFactors = FALSE)

  tmp_dir <- tempdir()

  ################################
  # Check all method
  ################################
  all_msp_new_pth <- file.path(tmp_dir,'all.msp')
  createMSP(pa, msp_file = all_msp_new_pth, metadata = metadata, method = "all", xcms_groupids = c(8, 12))

  all_msp_new <- get_msp_str(all_msp_new_pth)
  all_msp_old <- get_msp_str(system.file("extdata", "msp_test", "all.msp", package="msPurity"))
  expect_equal(all_msp_new, all_msp_old)


  ################################
  # Check max method
  ################################
  max_msp_new_pth <- file.path(tmp_dir,'max.msp')
  createMSP(pa, msp_file = max_msp_new_pth, metadata = metadata, method = "max", xcms_groupids = c(8, 12))

  max_msp_new <- get_msp_str(max_msp_new_pth)
  max_msp_old <- get_msp_str(system.file("extdata", "msp_test", "max.msp", package="msPurity"))
  expect_equal(max_msp_new, max_msp_old)


  ################################
  # Check av_inter method
  ################################
  av_inter_msp_new_pth <- file.path(tmp_dir,'av_inter.msp')
  createMSP(pa, msp_file = av_inter_msp_new_pth, metadata = metadata, method = "av_inter", xcms_groupids = c(8, 12))

  av_inter_msp_new <- get_msp_str(av_inter_msp_new_pth)
  av_inter_msp_old <- get_msp_str(system.file("extdata","msp_test", "av_inter.msp", package="msPurity"))
  expect_equal(av_inter_msp_new, av_inter_msp_old)


  ################################
  # Check av_intra method
  ################################
  av_intra_msp_new_pth <- file.path(tmp_dir,'av_intra.msp')
  createMSP(pa, msp_file = av_intra_msp_new_pth, metadata = metadata, method = "av_intra", xcms_groupids = c(8, 12))

  av_intra_msp_new <- get_msp_str(av_intra_msp_new_pth)
  av_intra_msp_old <- get_msp_str(system.file("extdata","msp_test", "av_intra.msp", package="msPurity"))
  expect_equal(av_intra_msp_new, av_intra_msp_old)

  ################################
  # Check av_all method
  ################################
  av_all_msp_new_pth <- file.path(tmp_dir,'av_all.msp')
  createMSP(pa, msp_file = av_all_msp_new_pth, metadata = metadata, method = "av_all", xcms_groupids = c(8, 12))

  av_all_msp_new <- get_msp_str(av_all_msp_new_pth)
  av_all_msp_old <- get_msp_str(system.file("extdata", "msp_test", "av_all.msp", package="msPurity"))
  expect_equal(av_all_msp_new, av_all_msp_old)


  #createMSP(pa, msp_file = 'all.msp', metadata = metadata, method = "all", xcms_groupids = c(8, 12))
  #createMSP(pa, msp_file = 'max.msp', metadata = metadata, method = "max", xcms_groupids = c(8, 12))
  #createMSP(pa, msp_file = 'av_inter.msp', metadata = metadata, method = "av_inter", xcms_groupids = c(8, 12))
  #createMSP(pa, msp_file = 'av_intra.msp', metadata = metadata, method = "av_intra", xcms_groupids = c(8, 12))
  #createMSP(pa, msp_file = 'av_all.msp', metadata = metadata, method = "av_all", xcms_groupids = c(8, 12))



})



