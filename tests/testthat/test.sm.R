context ("checking spectral matching functions with averaging")

test_that("checking spectral matching functions with averaging", {
  print ("\n")
  print("########################################################")
  print("## Spectral matching functions                        ##")
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
  q_dbPth <- system.file("extdata", "create_database_example.sqlite", package="msPurity")

  q_copyDbPth <- file.path(td,paste('lcmsms_data', format(Sys.time(), "%Y-%m-%d-%I%M%S"), 'copy.sqlite', sep="-"))
  file.copy(from = q_dbPth, to=q_copyDbPth)

  result <- spectral_matching(q_copyDbPth, spectra_type_q = 'av_all', match_alg = 'dpc')

  expect_equal(result$xcms_summary_df$grpid, c(12,27,33,42,49,61,65,76,213, 351))
  expect_equal(result$xcms_summary_df$best_name, c("proline", "Isoleucine", "N-methylnicotinate", "Acetylcholine",
                                                   "Oxypurinol", "4-Coumaric acid", "alpha-Methyl-DL-histidine",
                                                   "Tyrosine", "Aspartame", "Doxazosin"))
  expect_equal(round(as.numeric(result$xcms_summary_df$best_median_score), 3),
               c(1.000,0.929,0.988,0.949,0.797,0.615,0.922,0.774,0.718,0.955))


})
