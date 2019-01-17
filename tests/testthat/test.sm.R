# test_that("checking spectral matching functions (massbank)", {
#  too time consuming to run on bioconductor
#   print("########################################################")
#   print("## Spectral matching functions                        ##")
#   print("########################################################")
#   library(msPurity)
#
#   msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#   xset <- xcms::xcmsSet(msmsPths)
#   xset <- xcms::group(xset)
#   xset <- xcms::retcor(xset)
#   xset <- xcms::group(xset)
#
#   pa  <- purityA(msmsPths)
#   td <- tempdir()
#   pa <- frag4feature(pa, xset, out_dir = td, create_db = TRUE)
#
#   scan_ids <- c(1120,  366, 1190, 601,  404,1281, 1323, 1289)
#
#   result <- spectral_matching(pa@db_path, scan_ids =scan_ids)
#
#   expect_equal(unname(unlist(result$xcms_summary_df)),
#                c("46", "Methionine", "0.605663747378067", "0.125", "3",
#                  "Methionine, Methionine, Methionine", "1190, 1190, 1190", "47280, 53635, 53643",
#                  "0.686, 0.606, 0.601", "1", "scans"))
#
#
#   con <- DBI::dbConnect(RSQLite::SQLite(), file.path(result$result_db_pth))
#
#   XLI <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups
#                          LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.grpid=c_peak_groups.grpid
#                          LEFT JOIN c_peaks on c_peaks.cid=cXg.cid
#                          LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.cid=c_peaks.cid
#                          LEFT JOIN s_peak_meta ON cXs.pid=s_peak_meta.pid
#                          LEFT JOIN matches ON matches.pid=s_peak_meta.pid
#                          LEFT JOIN library_meta ON matches.lid=library_meta.lid
#                          WHERE matches.score IS NOT NULL')
#
#   expect_equal(ncol(XLI), 80)
#   expect_equal(nrow(XLI), 58)
#   expect_equal(round(sum(XLI$mz),3),  9737.548)
#   expect_equal(round(sum(XLI$precursorMZ),3), 9737.997)
#
#   expect_equal(XLI$score, c(0.031426723, 0.033767212, 0.037243842, 0.124010669,
#                             0.192138067, 0.192896407, 0.196154318, 0.202619470,
#                             0.164473357, 0.168898447, 0.294436615, 0.347514858,
#                             0.366409140, 0.419947544, 0.476714032, 0.505022823,
#                             0.518372745, 0.562079805, 0.577939912, 0.597344198,
#                             0.601224155, 0.605663747, 0.686332609, 0.012929117,
#                             0.013388544, 0.016620808, 0.017613880, 0.033050314,
#                             0.038243196, 0.040964938, 0.047862922, 0.050382537,
#                             0.050597154, 0.060602179, 0.095390574, 0.096959927,
#                             0.407464481, 0.037185179, 0.038463747, 0.009790092,
#                             0.117525062, 0.144347531, 0.155452807, 0.282820410,
#                             0.318886279, 0.324863487, 0.327210614, 0.407980001,
#                             0.447448267, 0.450909940, 0.457047449, 0.505938929,
#                             0.536365177, 0.064598670, 0.004899517, 0.005344766,
#                             0.006063552, 0.006463533))
#
#   expect_equal(XLI$perc_mtch, c(0.07692308, 0.07692308, 0.07692308, 0.05555556, 0.07692308,
#                                  0.08333333, 0.05263158, 0.05882353, 0.09523810, 0.20000000,
#                                  0.10000000, 0.26666667, 0.25000000, 0.12500000, 0.28571429,
#                                  0.23529412, 0.25000000, 0.20000000, 0.23076923, 0.17391304,
#                                  0.12500000, 0.11538462, 0.26666667, 0.06250000, 0.05263158,
#                                  0.05263158, 0.12500000, 0.08333333, 0.10000000, 0.07692308,
#                                  0.09090909, 0.08333333, 0.10000000, 0.09090909, 0.06666667,
#                                  0.07142857, 0.13333333, 0.06666667, 0.05263158, 0.03125000,
#                                  0.11111111, 0.06250000, 0.18181818, 0.18181818, 0.20000000,
#                                  0.16666667, 0.21052632, 0.33333333, 0.36363636, 0.30769231,
#                                  0.12903226, 0.36363636, 0.36363636, 0.07692308, 0.03030303,
#                                  0.03225806, 0.03571429, 0.03448276))
#
# })


test_that("checking spectral matching functions with averaging", {
  print("########################################################")
  print("## Spectral matching functions                        ##")
  print("########################################################")
  library(msPurity)

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  xset <- xcms::xcmsSet(msmsPths)
  xset <- xcms::group(xset)
  xset <- xcms::retcor(xset)
  xset <- xcms::group(xset)

  pa  <- purityA(msmsPths)
  pa <- frag4feature(pa, xset)
  pa  <- averageAllFragSpectra(pa)
  td <- tempdir()
  db_path <- create_database(pa, xset = xset, out_dir = ".")

  result <- spectral_matching(db_path, spectra_type_q = 'av_all', match_alg = 'dpc')


  #db_path <- create_database(pa, xset = xset, out_dir = ".")
  #result2 <- spectral_matching(db_path, spectra_type_q = 'av_all', match_alg = 'mf')

  expect_equal(result$xcms_summary_df$grpid, c(12,27,33,42,49,61,65,76,213, 351))
  expect_equal(result$xcms_summary_df$best_name, c("proline", "Isoleucine", "N-methylnicotinate", "Acetylcholine",
                                                   "Oxypurinol", "4-Coumaric acid", "alpha-Methyl-DL-histidine",
                                                   "Tyrosine", "Aspartame", "Doxazosin"))
  expect_equal(round(as.numeric(result$xcms_summary_df$best_median_score), 3),
               c(1.000,0.929,0.988,0.949,0.797,0.615,0.922,0.774,0.718,0.975))


  db_path <- create_database(pa, xset = xset, out_dir = ".")

  result <- spectral_matching(db_path, spectra_type_q = 'av_all',ra_thres_q = )

})
