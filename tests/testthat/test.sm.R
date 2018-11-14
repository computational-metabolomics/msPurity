test_that("checking spectral matching functions (massbank)", {
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
  td <- tempdir()
  pa <- frag4feature(pa, xset, out_dir = td, create_db = TRUE)

  scan_ids <- c(1120,  366, 1190, 601,  404,1281, 1323, 1289)

  result <- spectral_matching(pa@db_path, scan_ids =scan_ids)

  expect_equal(unname(unlist(result$xcms_summary_df)),
               c("46", "Methionine", "0.603443951203275", "0.120192307692308", "3",
                 "Methionine, Methionine", "1190, 1190", "34329, 34337", "0.606, 0.601", "1"))


  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(result$result_db_pth))

  XLI <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups
                         LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.grpid=c_peak_groups.grpid
                         LEFT JOIN c_peaks on c_peaks.cid=cXg.cid
                         LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.cid=c_peaks.cid
                         LEFT JOIN s_peak_meta ON cXs.pid=s_peak_meta.pid
                         LEFT JOIN matches ON matches.pid=s_peak_meta.pid
                         LEFT JOIN library_meta ON matches.lid=library_meta.lid
                         WHERE matches.score IS NOT NULL')

  expect_equal(ncol(XLI), 80)
  expect_equal(nrow(XLI), 45)
  expect_equal(round(sum(XLI$mz),3), 7626.806)
  expect_equal(round(sum(XLI$precursorMZ),3), 7627.127)

  expect_equal(XLI$score, c(0.031426723, 0.033767212, 0.037243842, 0.124010669, 0.192138067,
                            0.192896407, 0.196154318, 0.202619470, 0.168898447, 0.347514858,
                            0.366409140, 0.419947544, 0.476714032, 0.505022823, 0.518372745,
                            0.562079805, 0.597344198, 0.601224155, 0.605663747, 0.012929117,
                            0.013388544, 0.016620808, 0.033050314, 0.038243196, 0.040964938,
                            0.047862922, 0.050382537, 0.050597154, 0.407464481, 0.037185179,
                            0.038463747, 0.155452807, 0.318886279, 0.324863487, 0.327210614,
                            0.407980001, 0.447448267, 0.450909940, 0.505938929, 0.536365177,
                            0.064598670, 0.004899517, 0.005344766, 0.006063552, 0.006463533))

  expect_equal(XLI$perc_mtch, c(0.07692308, 0.07692308, 0.07692308, 0.05555556,
                                0.07692308, 0.08333333, 0.05263158, 0.05882353,
                                0.20000000, 0.26666667, 0.25000000, 0.12500000,
                                0.28571429, 0.23529412, 0.25000000, 0.20000000,
                                0.17391304, 0.12500000, 0.11538462, 0.06250000,
                                0.05263158, 0.05263158, 0.08333333, 0.10000000,
                                0.07692308, 0.09090909, 0.08333333, 0.10000000,
                                0.13333333, 0.06666667, 0.05263158, 0.18181818,
                                0.20000000, 0.16666667, 0.21052632, 0.33333333,
                                0.36363636, 0.30769231, 0.36363636, 0.36363636,
                                0.07692308, 0.03030303, 0.03225806, 0.03571429,
                                0.03448276))

})

