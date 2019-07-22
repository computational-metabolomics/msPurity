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
  localCompoundDbPth <- system.file("extdata", "tests", "db", "compounds_18July2019_0319.sqlite", package="msPurity")
  td <- tempdir()
  sm_resultPthCopy <- file.path(td, 'sm_result_tmp16.sqlite')
  file.copy(sm_resultPth, sm_resultPthCopy)

  sqlitePth <- sm_resultPthCopy

  combined <- combineAnnotations(sqlitePth, metfrag_resultPth, sirius_csi_resultPth, probmetab_resultPth,
                                 localCompoundDbPth=localCompoundDbPth)
  combined <- combineAnnotations(sqlitePth, metfrag_resultPth, sirius_csi_resultPth, probmetab_resultPth)


})
