context ("checking combine annotations based functions")

test_that("checking combine annotations based functions", {
  print ("\n")
  print("########################################################")
  print("## Checking combine annotations based functions       ##")
  print("########################################################")
  metfrag_resultPth <- system.file("extdata", "tests", "external_annotations", "metfrag.tsv", package="msPurity")
  sirius_csi_resultPth <- system.file("extdata","tests", "external_annotations", "sirus_csifingerid.tsv", package="msPurity")
  probmetab_resultPth <- system.file("extdata","tests", "external_annotations", "probmetab.tsv", package="msPurity")

  ms1_lookup_resultPth <- '/home/tomnl/Downloads/Galaxy810-[BEAMS_on_data_808,_data_160,_and_data_159__Summary_table_(Multiple_rows_and_separate_columns)].tsv'
  sm_resultPth <- system.file("extdata","tests", "sm", "spectralMatching_result.sqlite", package="msPurity")
  localCompoundDbPth <- '/media/rds-dma/external_data/pubchem_sqlite/metab_compound_test2.sqlite'

  #localCompoundDbPth <- system.file("extdata", "tests", "db", "compounds_18July2019_0319.sqlite", package="msPurity")
  weights=list('sm'=0.3,
               'metfrag'=0.2,
               'sirius_csifingerid'=0.2,
               'probmetab'=0,
               'ms1_lookup'=0.05,
               'biosim'=0.25
  )
  td <- tempdir()
  sm_resultPthCopy <- file.path(td, 'sm_result_tmp40.sqlite')
  file.copy(sm_resultPth, sm_resultPthCopy)

  sqlitePth <- sm_resultPthCopy

  combined <- combineAnnotations(sqlitePth,
                                 localCompoundDbPth=localCompoundDbPth,
                                 metfrag_resultPth = metfrag_resultPth,
                                 sirius_csi_resultPth = sirius_csi_resultPth,
                                 probmetab_resultPth = probmetab_resultPth,
                                 ms1_lookup_resultPth = ms1_lookup_resultPth,
                                 weights = weights,
                                 ms1_lookup_dbSource='hmdb',
                                 ms1_lookup_checkAdducts=FALSE,
                                 ms1_lookup_keepAdducts=c('[M+H]+', '[M-H]-'),
                                 )
  # likely to change to only providing limited testing
  expect_equal(combined$rank, c(1,2,3,4,4,4,4,4,1,2,1,2,3,4,1))
})
