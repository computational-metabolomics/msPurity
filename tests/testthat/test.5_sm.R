context ("checking spectral matching functions ")

test_that("checking spectral matching functions (spectral_matching)", {
  print ("\n")
  print("########################################################")
  print("## Testing spectral matching (spectral_matching)      ##")
  print("########################################################")
  td <- tempdir()

  #msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  #xset <- xcms::xcmsSet(msmsPths)
  #xset <- xcms::group(xset)
  #xset <- xcms::retcor(xset)
  #xset <- xcms::group(xset)

  #pa  <- purityA(msmsPths)
  #pa <- frag4feature(pa, xset)
  #pa  <- averageAllFragSpectra(pa)

  #q_copyDbPth <- create_database(pa, xset = xset, out_dir = td)
  q_dbPth <- system.file("extdata", "tests", "db", "create_database_example.sqlite", package="msPurity")
  q_copyDbPth <- file.path(td, paste0(format(Sys.time(), "%Y-%m-%d-%I%M%S"), 'copy_spectral_matching.sqlite'))
  file.copy(from = q_dbPth, to=q_copyDbPth)

  result <- spectral_matching(q_copyDbPth, spectra_type_q = 'av_all', match_alg = 'dpc')

  expect_equal(result$xcms_summary_df$grpid, c(12, 27, 33, 42, 49, 65, 76, 213, 351))
  expect_equal(result$xcms_summary_df$best_name, c("proline", "Isoleucine", "N-methylnicotinate",
                                                   "Acetylcholine", "Oxypurinol", "alpha-Methyl-DL-histidine",
                                                   "Tyrosine", "Aspartame", "Doxazosin"  ))
  expect_equal(round(as.numeric(result$xcms_summary_df$best_median_score), 3),
               c(1.000, 0.929, 0.988, 0.949, 0.800, 0.922, 0.707, 0.718, 0.976))


})


test_that("checking spectral matching functions (spectralMatching) query vs library", {
  print ("\n")
  print("########################################################")
  print("## Testing spectral matching (spectralMatching)  qvl  ##")
  print("########################################################")

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
  #

  td <- tempdir()
  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="msPurity")

  rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
  sm_out_pth <- file.path(td, rid)

  result <- spectralMatching(q_dbPth,
                             q_xcmsGroups = c(12, 27),
                             cores = 1,
                             q_pol = NA,
                             l_accessions=c('CCMSLIB00000577898', 'CE000616'),
                             l_pol='positive',
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

})



test_that("checking spectral matching functions (spectralMatching) library vs library", {
  print ("\n")
  print("########################################################")
  print("## Testing spectral matching (spectralMatching)  lvl  ##")
  print("########################################################")

  # q_purity=NA
  # q_ppmProd=10
  # q_ppmPrec=5
  # q_raThres=NA
  # q_pol=NA
  # q_instrumentTypes=NA
  # q_instrument=NA
  # q_sources=NA
  # q_spectraTypes=NA
  # q_pids=c(1,2,3)
  # q_rtrange=c(NA, NA)
  # q_spectraFilter=TRUE
  # q_xcmsGroups=NA
  # q_accessions=NA
  #
  # l_purity=NA
  # l_ppmProd=10
  # l_ppmPrec=5
  # l_raThres=NA
  # l_pol=NA
  # l_instrumentTypes=NA
  # l_instrument=NA
  # l_sources=NA
  # l_spectraTypes=NA
  # l_pids=c(1,2,3)
  # l_rtrange=c(NA, NA)
  # l_spectraFilter=FALSE
  # l_xcmsGroups=NA
  # l_accessions=NA
  # usePrecursors=TRUE
  # raW=0.5
  # mzW=2
  # rttol=NA
  # updateDb = TRUE
  # copyDb = TRUE
  # outPth = sm_out_pth

  td <- tempdir()

  rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
  sm_out_pth <- file.path(td, rid)


  q_dbPth <- system.file("extdata", "library_spectra", "library_spectra.db", package="msPurityData")
  l_dbPth <- system.file("extdata", "library_spectra", "library_spectra.db", package="msPurityData")

  result <- spectralMatching(q_dbPth=q_dbPth,
                             l_dbPth=l_dbPth,
                             q_pids = c(1,2,3),
                             q_spectraTypes = NA,
                             cores=1,
                             q_pol = NA,
                             l_pids = c(1,2,3),
                             l_spectraTypes = NA,
                             updateDb = TRUE,
                             copyDb = TRUE,
                             outPth = sm_out_pth)

  matched <- result$matchedResults

  expect_equal(matched$lpid, c(1,2,3))
  expect_equal(matched$dpc, c(1,1,1))
  expect_equal(matched$rdpc, c(1,1,1))
  expect_equal(matched$cdpc, c(1,1,1))
  expect_equal(matched$mcount, c(43,9,321))
  expect_equal(matched$allcount, c(43,9,321))
  expect_equal(matched$accession, c("CCMSLIB00000001577", "CCMSLIB00000001666", "CCMSLIB00000004593"))
  expect_equal(matched$name, c("Carmaphycin B", "CarmaphycinA_enone", "MaprotilineHCl"))
  expect_equal(matched$qpid, c(1,2,3))


})

test_that("checking spectral matching functions (spectralMatching) query vs query", {
  print ("\n")
  print("########################################################")
  print("## Testing spectral matching (spectralMatching)   qvq ##")
  print("########################################################")
  td <- tempdir()

  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="msPurity")
  l_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="msPurity")

  rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
  sm_out_pth <- file.path(td, rid)


  result <- spectralMatching(q_dbPth=q_dbPth,
                             l_dbPth=l_dbPth,
                             q_xcmsGroups = c(12, 27),
                             q_spectraTypes = 'av_all',
                             q_spectraFilter = TRUE,
                             q_pol = NA,
                             l_xcmsGroups = c(12, 27),
                             l_spectraTypes = 'av_all',
                             l_pol = NA,
                             l_spectraFilter = TRUE,
                             cores=1,
                             updateDb = TRUE,
                             copyDb = TRUE,
                             outPth = sm_out_pth)

  matched <- result$matchedResults
  expect_equal(matched$lpid, c(1666, 1670))
  expect_equal(matched$qpid, c(1666, 1670))
  expect_equal(matched$dpc, c(1,1))
  expect_equal(matched$rdpc, c(1,1))
  expect_equal(matched$cdpc, c(1,1))
  expect_equal(matched$mcount, c(3, 1))
  expect_equal(matched$allcount, c(3, 1))

  xcmsMatchedResults <- result$xcmsMatchedResults

  expect_equal(xcmsMatchedResults$pid, c(1666, 1670))
  expect_equal(xcmsMatchedResults$lpid, c(1666, 1670))
  expect_equal(xcmsMatchedResults$grpid, c(12, 27))

})
