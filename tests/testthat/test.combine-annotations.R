test_that("checking lcms based functions", {
  print("########################################################")
  print("## Checking LCMS based class and functions            ##")
  print("########################################################")
  metfrag_resultPth <- system.file("extdata", "external_annotations", "metfrag.tsv", package="msPurity")
  sirius_csi_resultPth <- system.file("extdata", "external_annotations", "sirus_csifingerid.tsv", package="msPurity")
  probmetab_resultPth <- system.file("extdata", "external_annotations", "probmetab.tsv", package="msPurity")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  xset <- xcms::xcmsSet(msmsPths)
  xset <- xcms::group(xset)
  xset <- xcms::retcor(xset)
  xset <- xcms::group(xset)
  pa  <- purityA(msmsPths)
  pa <- frag4feature(pa, xset)
  pa <- averageAllFragSpectra(pa)
  td <- tempdir()
  db_path <- create_database(pa, xset = xset, out_dir = td)
  result <- spectral_matching(db_path, spectra_type_q="av_all")
  combined <- combineAnnotations(db_path, metfrag_resultPth, sirius_csi_resultPth, probmetab_resultPth)
  #combined$sirius_score <- round(as.numeric(combined$sirius_score), 2)
  #combined <- data.frame(lapply(combined , as.character), stringsAsFactors=FALSE)
  #write.table(combined ,'inst/extdata/external_annotations/combined.tsv', sep='\t', row.names = FALSE, col.names = TRUE )


  combinedOLD <- read.table(system.file("extdata", "external_annotations", "combined.tsv", package="msPurity"),
                            header = TRUE, stringsAsFactors = FALSE)
  combinedOLD <- data.frame(lapply(combinedOLD , as.character), stringsAsFactors=FALSE)
  expect_equal(combined, combinedOLD)

})



