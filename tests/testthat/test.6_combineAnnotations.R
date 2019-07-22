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
  weights=list('sm'=0.2,
               'metfrag'=0.15,
               'sirius_csifingerid'=0.15,
               'cfm'=0.15,
               'probmetab'=0.05,
               'biosim'=0.3
  )
  td <- tempdir()
  sm_resultPthCopy <- file.path(td, 'sm_result_tmp37.sqlite')
  file.copy(sm_resultPth, sm_resultPthCopy)

  sqlitePth <- sm_resultPthCopy

  combined <- combineAnnotations(sqlitePth,
                                 localCompoundDbPth=localCompoundDbPth,
                                 metfrag_resultPth = metfrag_resultPth,
                                 sirius_csi_resultPth = sirius_csi_resultPth,
                                 probmetab_resultPth = probmetab_resultPth,
                                 weights = weights
                                 )
  # likely to change to only providing limited testing
  expect_equal(combined$rank, c(1,2,3,4,4,4,4,4,1,2,1,2,3,4,1))
})
