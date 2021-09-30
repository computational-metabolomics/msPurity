context ("checking spectral matching functions (xcms v2 functions) ")


test_that("checking spectral matching functions (spectralMatching) query vs library (xcms v2 functions)", {
  print ("\n")
  print("############################################################################")
  print("## Testing spectral matching (spectralMatching)  qvl (xcms v2 functions)  ##")
  print("############################################################################")

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
  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example_OLD.sqlite", package="msPurity")

  rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
  sm_out_pth <- file.path(td, rid)

  result <- spectralMatching(q_dbPth, q_xcmsGroups = c(12, 27), cores=1, l_accessions=c('CCMSLIB00000577898','CE000616'),
                             q_spectraTypes = 'av_all',
                             updateDb = TRUE,
                             copyDb = TRUE,
                             outPth = sm_out_pth)

  expect_equal(result$xcmsMatchedResults$grpid, c(12,27))
  expect_equal(result$xcmsMatchedResults$library_accession, c("CCMSLIB00000577898", "CE000616"))
  expect_equal(result$xcmsMatchedResults$inchikey, c("ONIBWKKTOPOVIA-UHFFFAOYSA-N", "AGPKZVBTJJNPAG-UHFFFAOYSA-N"))

  expect_equal(round(as.numeric(result$xcmsMatchedResults$dpc),3), c( 0.879, 0.941))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$rdpc),3), c( 0.881,  0.941))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$cdpc),3), c( 0.857, 0.896))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$mcount),3), c( 2, 1))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$allcount),3), c( 23, 20))

})




test_that("checking spectral matching functions (spectralMatching) query vs query (xcms v2 functions)", {
  print ("\n")
  print("############################################################################")
  print("## Testing spectral matching (spectralMatching)   qvq (xcms v2 functions) ##")
  print("############################################################################")
  td <- tempdir()

  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example_OLD.sqlite", package="msPurity")
  l_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example_OLD.sqlite", package="msPurity")

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
                             usePrecursors=TRUE,
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

test_that("checking spectral matching functions (spectralMatching) query(scans) vs library (xcms v2 functions)", {
  print ("\n")
  print("######################################################################################")
  print("## Testing spectral matching (spectralMatching)   q(scan) v q  (xcms v2 functions)  ##")
  print("######################################################################################")
  td <- tempdir()

  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example_OLD.sqlite", package="msPurity")

  rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
  sm_out_pth <- file.path(td, rid)


  result <- spectralMatching(q_dbPth, q_xcmsGroups = c(12, 27), cores=1, l_accessions=c('CCMSLIB00000577898','CE000616'),
                             q_spectraTypes = 'scan',
                             updateDb = TRUE,
                             copyDb = TRUE,
                             outPth = sm_out_pth)

  matched <- result$matchedResults
  expect_equal(matched$lpid, c(5325, 5325, 5325, 5325, 5325, 5325, 53807, 53807, 53807, 53807, 53807, 53807))
  expect_equal(matched$qpid, c(226, 281, 336, 1110, 1166, 1055, 509, 394, 451, 1220, 1275, 1330))
  expect_equal(round(matched$dpc,2), c(0.14, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(round(matched$rdpc,2), c( 0.87 , 0.87 ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ,NaN ))
  expect_equal(round(matched$cdpc,2), c(0.13, 0.03, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(matched$mcount, c(1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  expect_equal(matched$allcount, c(34, 46, 56, 53, 37, 43, 34, 59, 30, 24, 40, 29))

  xcmsMatchedResults <- result$xcmsMatchedResults

  expect_equal(xcmsMatchedResults$grpid, c(12, 12, 12, 12, 12, 12, 27, 27, 27, 27, 27, 27))

})

