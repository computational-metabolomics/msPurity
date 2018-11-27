test_that("checking lcms based functions", {
  print("########################################################")
  print("## Checking LCMS based class and functions            ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")


  xset <- xcms::xcmsSet(msmsPths)
  xset <- xcms::group(xset)


  pa  <- purityA(msmsPths)
  pa <- frag4feature(pa, xset, create_db=FALSE)
  pa <- averageFragmentation(pa)

  metadata <- data.frame('grpid'=8, 'isotope'=NA, 'MS$FOCUSED_ION: PRECURSOR_TYPE'='[M+H]+',
                         'AC$MASS_SPECTROMETRY: ION_MODE'="POSITIVE", 'CH$NAME'='Sulfamethizole',
                         check.names = F, stringsAsFactors = F)

  createMSP(pa, msp_file = 'test7.msp', metadata = metadata, method = "all")

  createMSP(pa, msp_file = 'test9.msp', metadata = metadata, method = "max")

})


