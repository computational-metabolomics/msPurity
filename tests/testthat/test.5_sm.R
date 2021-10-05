context ("checking spectral matching functions ")

test_that("checking spectral matching functions (spectralMatching) query vs library", {
  print ("\n")
  print("########################################################")
  print("## Testing spectral matching (spectralMatching)  qvl  ##")
  print("########################################################")

  # library(msPurity)
  # library(xcms)
  # library(magrittr)
  #
  # #### read in files and data
  # msPths <- dirname(list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = 'MSMS'))
  # msPths[1] <- file.path(msPths[1], 'LCMSMS_1.mzML')
  # msPths[2] <- file.path(msPths[2], 'LCMSMS_2.mzML')
  # ms_data = readMSData(msPths, mode = 'onDisk', msLevel. = 1)
  #
  # ##### subset the data to focus on retention times 30-90 seconds and m/z values between 100 and 200 m/z.
  # rtr = c(30, 90)
  # mzr = c(100, 200)
  # ms_data = ms_data %>%  filterRt(rt = rtr) %>%  filterMz(mz = mzr)
  #
  # ##### perform feature detection in individual files
  # cwp <- CentWaveParam(snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30))
  # xcmsObj <- xcms::findChromPeaks(ms_data, param = cwp)
  # xcmsObj@phenoData@data$class = c('sample', 'sample')
  # xcmsObj@phenoData@varMetadata = data.frame('labelDescription' = 'sampleNames', 'class')
  # pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = 0, bw = 5, binSize = 0.017)
  # xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)
  #
  # td <- tempdir()
  #
  # pa  <- purityA(msPths)
  # pa <- frag4feature(pa, xcmsObj)
  # pa <- filterFragSpectra(pa, allfrag = TRUE)
  # pa  <- averageAllFragSpectra(pa)
  # q_copyDbPth <- createDatabase(pa, xcmsObj = xcmsObj, outDir = td, metadata=list('polarity'='positive', 'instrument'='Q-Exactive'))

  td <- tempdir()
  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="msPurity")

  rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
  sm_out_pth <- file.path(td, rid)

  result <- spectralMatching(q_dbPth, q_xcmsGroups = c(53, 89, 410), cores=1, l_accessions = c('PR100407', 'ML005101', 'CCMSLIB00003740024'),
                             q_spectraTypes = 'av_all',
                             updateDb = TRUE,
                             copyDb = TRUE,
                             outPth = sm_out_pth)

  expect_equal(result$xcmsMatchedResults$grpid, c(53, 89, 410))
  expect_equal(result$xcmsMatchedResults$library_accession, c("ML005101", "CCMSLIB00003740024", "PR100407" ))
  expect_equal(result$xcmsMatchedResults$inchikey, c("DFPAKSUCGFBDDF-UHFFFAOYSA-N", "CKLJMWTZIZZHCS-UHFFFAOYSA-N", "YHHSONZFOIEMCP-UHFFFAOYSA-O"))

  expect_equal(round(as.numeric(result$xcmsMatchedResults$dpc),3), c(0.934, 0.507, 0.969))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$rdpc),3), c(  0.965, 0.507, 0.972))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$cdpc),3), c(0.781, 0.531, 0.963))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$mcount),3), c( 2, 2, 2))
  expect_equal(round(as.numeric(result$xcmsMatchedResults$allcount),3), c( 7, 34, 8))

})



test_that("checking spectral matching functions (spectralMatching) library vs library", {
  print ("\n")
  print("########################################################")
  print("## Testing spectral matching (spectralMatching)  lvl  ##")
  print("########################################################")
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
  expect_equal(matched$library_accession, c("CCMSLIB00000001577", "CCMSLIB00000001666", "CCMSLIB00000004593"))
  expect_equal(matched$library_entry_name, c("Carmaphycin B", "CarmaphycinA_enone", "MaprotilineHCl"))
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
                             q_xcmsGroups = c(89, 410),
                             q_spectraTypes = 'av_all',
                             q_spectraFilter = TRUE,
                             q_pol = NA,
                             l_xcmsGroups = c(89, 410),
                             l_spectraTypes = 'av_all',
                             l_pol = NA,
                             l_spectraFilter = TRUE,
                             cores=1,
                             updateDb = TRUE,
                             copyDb = TRUE,
                             usePrecursors=TRUE,
                             outPth = sm_out_pth)

  matched <- result$matchedResults
  expect_equal(matched$lpid, c(1661, 1677))
  expect_equal(matched$qpid, c(1661, 1677))
  expect_equal(matched$dpc, c(1,1))
  expect_equal(matched$rdpc, c(1,1))
  expect_equal(matched$cdpc, c(1,1))
  expect_equal(matched$mcount, c(2, 5))
  expect_equal(matched$allcount, c(2, 5))

  xcmsMatchedResults <- result$xcmsMatchedResults

  expect_equal(xcmsMatchedResults$grpid, c(89, 410))

})

test_that("checking spectral matching functions (spectralMatching) query(scans) vs library", {
  print ("\n")
  print("#################################################################")
  print("## Testing spectral matching (spectralMatching)   q(scan) v q ##")
  print("#################################################################")
  td <- tempdir()

  q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="msPurity")

  rid <- paste0(paste0(sample(LETTERS, 5, TRUE), collapse=""),  paste0(sample(9999, 1, TRUE), collapse=""), ".sqlite")
  sm_out_pth <- file.path(td, rid)


  # Note we get considerably better sm results when using averaged spectra for these xcms groups - this could just reflect the processing that was performed
  result <- spectralMatching(q_dbPth, q_xcmsGroups = c(53, 89, 410), cores=1, l_accessions = c('CCMSLIB00000479720', 'KNA00052', 'PR100407'),
                             q_spectraTypes = 'scan',
                             updateDb = TRUE,
                             copyDb = TRUE,
                             outPth = sm_out_pth)

  matched <- result$matchedResults
  expect_equal(matched$lpid, c(5203, 53964, 57016, 57016))
  expect_equal(matched$qpid, c( 391, 198, 205, 1034))
  expect_equal(round(matched$dpc,2), c(0.00, 0.00, 0.00, 0.76))

  expect_equal(matched$mcount, c(0, 0, 0, 2))
  expect_equal(matched$allcount, c(17, 10,  15, 10))

  xcmsMatchedResults <- result$xcmsMatchedResults

  expect_equal(xcmsMatchedResults$grpid, c(53,89, 410, 410))

})
