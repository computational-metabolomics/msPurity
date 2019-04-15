context ("checking spectral matching functions ")

test_that("checking spectral matching functions (spectral_matching)", {
  print ("\n")
  print("########################################################")
  print("## Spectral matching functions (spectral_matching)    ##")
  print("########################################################")
  library(msPurity)
  td <- tempdir()

  # msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  # xset <- xcms::xcmsSet(msmsPths)
  # xset <- xcms::group(xset)
  # xset <- xcms::retcor(xset)
  # xset <- xcms::group(xset)
  #
  # pa  <- purityA(msmsPths)
  # pa <- frag4feature(pa, xset)
  # pa  <- averageAllFragSpectra(pa)

  # db_path <- create_database(pa, xset = xset, out_dir = td)
  q_dbPth <- system.file("extdata", "tests", "db", "create_database_example.sqlite", package="msPurity")

  q_copyDbPth <- file.path(td, paste0(format(Sys.time(), "%Y-%m-%d-%I%M%S"), 'copy.sqlite'))
  file.copy(from = q_dbPth, to=q_copyDbPth)

  result <- spectral_matching(q_copyDbPth, spectra_type_q = 'av_all', match_alg = 'dpc')

  expect_equal(result$xcms_summary_df$grpid, c(12, 27, 33, 42, 49, 65, 76, 213, 351))
  expect_equal(result$xcms_summary_df$best_name, c("proline", "Isoleucine", "N-methylnicotinate",
                                                   "Acetylcholine", "Oxypurinol", "alpha-Methyl-DL-histidine",
                                                   "Tyrosine", "Aspartame", "Doxazosin"  ))
  expect_equal(round(as.numeric(result$xcms_summary_df$best_median_score), 3),
               c(1.000, 0.929, 0.988, 0.949, 0.800, 0.922, 0.707, 0.718, 0.976))


})


test_that("checking spectral matching functions (spectralMatching", {
  print ("\n")
  print("########################################################")
  print("## Spectral matching functions                        ##")
  print("########################################################")
  td <- tempdir()
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

  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="msPurity")

  sm_out_pth <- file.path(td, 'sm_out.sqlite')

  result <- spectralMatching(q_dbPth,
                             q_xcmsGroups = c(12, 27),
                             cores=1,
                             q_pol = NA,
                             l_accessions=c('CCMSLIB00000577898', 'CE000616'),
                             updateDb = TRUE,
                             copyDb = TRUE,
                             outPth = sm_out_pth)

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
