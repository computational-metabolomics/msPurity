context ("checking spectralMatching (updated) functions")

test_that("checking spectralMatching (updated) functions", {
  print ("\n")
  print("########################################################")
  print("## Spectral matching functions                        ##")
  print("########################################################")
  library(msPurity)
  # msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  # xset <- xcms::xcmsSet(msmsPths)
  # xset <- xcms::group(xset)
  # xset <- xcms::retcor(xset)
  # xset <- xcms::group(xset)
  #
  # pa  <- purityA(msmsPths)
  # pa <- frag4feature(pa, xset)
  # pa <- filterFragSpectra(pa, allfrag=TRUE, plim=0, ilim = 0, snr = 0)
  # pa <- averageIntraFragSpectra(pa)
  # pa <- averageInterFragSpectra(pa)
  # pa <- averageAllFragSpectra(pa)
  # q_dbPth <- createDatabase(pa, xset, metadata=list('polarity'='positive', 'instrument'='Q-Exactive'))

  q_dbPth <- system.file("extdata", "createDatabase_example.sqlite", package="msPurity")

  result <- spectralMatching(q_dbPth, q_xcmsGroups = c(12, 27), cores=1, l_accessions=c('CCMSLIB00000577898',
                                                                                        'CE000616'))

  expect_equal(result$xcmsMatchedResults$grpid, c(12,27))
  expect_equal(result$xcmsMatchedResults$accession, c("CCMSLIB00000577898", "CE000616"))
  expect_equal(result$xcmsMatchedResults$inchikey, c("ONIBWKKTOPOVIA-UHFFFAOYSA-N", "AGPKZVBTJJNPAG-UHFFFAOYSA-N"))

  expect_equal(round(as.numeric(result$xcmsMatchedResults$dpc),3), c( 0.879, 0.941))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$rdpc),3), c( 0.995, 1.000))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$cdpc),3), c( 0.857, 0.896))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$mcount),3), c( 2, 1))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$allcount),3), c( 23, 20))
#
#   l_dbPth=NA
#   q_purity=NA
#   q_ppmProd=10
#   q_ppmPrec=5
#   q_raThres=0
#   q_pol='positive'
#   q_instrumentTypes=NA
#   q_instrument=NA
#   q_sources=NA
#   q_spectraType="av_all"
#   q_sids=NA
#   q_rtrange=c(NA, NA)
#   q_spectraFilter=TRUE
#   q_xcmsGroups=NA
#   q_pids = c(1659:1665)
#
#   l_purity=NA
#   l_ppmProd=10
#   l_ppmPrec=5
#   l_raThres=0
#   l_pol='positive'
#   l_instrumentTypes=NA
#   l_instrument=NA
#   l_sources=NA
#   l_spectraType=NA
#   l_sids=NA
#   l_rtrange=c(NA, NA)
#   l_spectraFilter=NA
#   l_xcmsGroups=NA
#   l_pids=NA
#
#
#
#
#
#   score_thres=0.6
#   raW=0.5
#   mzW=2
#   rttol=NA
#   match_alg='dpc'
#   cores=1
#   out_dir='.'
#   copy=FALSE
})
