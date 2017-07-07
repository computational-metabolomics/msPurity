test_that("checking spectral matching functions", {
  print("########################################################")
  print("## Spectral matching functions                        ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
  xset <- xcms::group(xset)
  xset <- xcms::retcor(xset)
  xset <- xcms::group(xset)

  pa  <- purityA(msmsPths, interpol = "linear")
  pa <- frag4feature(pa, xset)
  result <- spectral_matching(pa, xset, out_dir = tempdir())



  expect_equal(unname(unlist(result$xcms_summary_df)),
               c("46", "Methionine", "0.603443951203275", "0.120192307692308", "3", "Methionine, Methionine","1190, 1190",
                "39362, 38959", "0.606, 0.601",  "1"))

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(result$out_dir, result$db_name))

  XLI <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups
   LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.grpid=c_peak_groups.grpid
   LEFT JOIN c_peaks on c_peaks.cid=cXg.cid
   LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.cid=c_peaks.cid
   LEFT JOIN s_peak_meta ON cXs.pid=s_peak_meta.pid
   LEFT JOIN matches ON matches.pid=s_peak_meta.pid
   LEFT JOIN library_meta ON matches.lid=library_meta.lid
   WHERE matches.score IS NOT NULL')

  expect_equal(ncol(XLI), 67)
  expect_equal(nrow(XLI), 43)



})

