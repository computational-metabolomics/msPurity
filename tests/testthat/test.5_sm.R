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


  expect_warning(

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
                               outPth = sm_out_pth),


    class = "deprecatedWarning")

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

test_that("bundled query and reference spectra reproduce the three library hits", {
    fixture <- function(dir, name)
        system.file("extdata", "tests", dir, name, package = "msPurity",
                    mustWork = TRUE)
    hits <- spectralMatching(fixture("db", "createDatabase_example.sqlite"),
        fixture("library", "mini_library.sqlite"), q_xcmsGroups = c(53, 89, 410),
        l_accessions = c("PR100407", "ML005101", "CCMSLIB00003740024"),
        q_spectraTypes = "av_all", cores = 1L)$xcmsMatchedResults
    expect_equal(hits$grpid, c(53, 89, 410))
    expect_equal(hits$library_accession,
                 c("ML005101", "CCMSLIB00003740024", "PR100407"))
    expect_equal(round(hits$dpc, 3), c(0.934, 0.507, 0.969))
})
