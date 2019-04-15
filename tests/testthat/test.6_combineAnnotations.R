context ("checking combine annotations based functions")

test_that("checking combine annotations based functions", {
  print ("\n")
  print("########################################################")
  print("## Checking combine annotations based functions       ##")
  print("########################################################")
  metfrag_resultPth <- system.file("extdata", "tests", "external_annotations", "metfrag.tsv", package="msPurity")
  sirius_csi_resultPth <- system.file("extdata","tests", "external_annotations", "sirus_csifingerid.tsv", package="msPurity")
  probmetab_resultPth <- system.file("extdata","tests", "external_annotations", "probmetab.tsv", package="msPurity")
  sm_resultPth <- system.file("extdata","tests", "sm", "spectralMatching_result.sqlite", package="msPurity")

  td <- tempdir()
  sm_resultPthCopy <- file.path(td, 'sm_result_tmp.sqlite')
  file.copy(sm_resultPth, sm_resultPthCopy)
  #msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  #xset <- xcms::xcmsSet(msmsPths)
  #xset <- xcms::group(xset)
  #xset <- xcms::retcor(xset)
  #xset <- xcms::group(xset)
  #pa  <- purityA(msmsPths)
  #pa <- frag4feature(pa, xset)
  #pa <- averageAllFragSpectra(pa)

  #db_path <- create_database(pa, xset = xset, out_dir = td)
  #result <- spectral_matching(db_path, spectra_type_q="av_all")


  combined <- combineAnnotations(sm_resultPthCopy, metfrag_resultPth, sirius_csi_resultPth, probmetab_resultPth)

  #write.table(combined ,'inst/extdata/tests/external_annotations/combined.tsv', sep='\t', row.names = FALSE, col.names = TRUE )


  combinedOLD <- read.table(system.file("extdata", "tests", "external_annotations", "combined.tsv", package="msPurity"),
                            header = TRUE, stringsAsFactors = FALSE)
  expect_equal(as.numeric(combined$grpid), as.numeric(combinedOLD$grpid))
  expect_equal(as.numeric(combined$sirius_id), as.numeric(combinedOLD$sirius_id))
  expect_equal(as.numeric(combined$sirius_wscore), as.numeric(combinedOLD$sirius_wscore))
  expect_equal(as.numeric(combined$metfrag_id), as.numeric(combinedOLD$metfrag_id))
  expect_equal(as.numeric(combined$metfrag_wscore), as.numeric(combinedOLD$metfrag_wscore))
  expect_equal(as.numeric(combined$metfrag_wscore), as.numeric(combinedOLD$metfrag_wscore))
  expect_equal(as.numeric(combined$wscore), as.numeric(combinedOLD$wscore))
  expect_equal(as.numeric(combined$probmetab_wscore), as.numeric(combinedOLD$probmetab_wscore))

})



