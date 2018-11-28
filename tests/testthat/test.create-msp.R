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

  metadata <- data.frame('grpid'=8, 'isotope'=NA, 'MS$FOCUSED_ION: PRECURSOR_TYPE'='[M+H]+',
                         'AC$MASS_SPECTROMETRY: ION_MODE'="POSITIVE", 'CH$NAME'='Sulfamethizole',
                         check.names = FALSE, stringsAsFactors = FALSE)

  createMSP(pa, msp_file = 'test68.msp', metadata = metadata, method = "all")

  createMSP(pa, msp_file = 'test69.msp', metadata = metadata, method = "max")

  createMSP(pa, msp_file = 'test70.msp', metadata = metadata, method = "av_inter")

  createMSP(pa, msp_file = 'test71.msp', metadata = metadata, method = "av_intra")

  createMSP(pa, msp_file = 'test72.msp', metadata = metadata, method = "av_all")



})


