# context ("checking database functions")
#
# test_that("checking database functions", {
#   print ("\n")
#   print("########################################################")
#   print("## Checking database handles lcms and msms            ##")
#   print("########################################################")
#
#   msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#   library(xcms)
#
#   pa <- purityA(msmsPths)
#   pa <- frag4feature(pa, xset)
#   pa <- averageIntraFragSpectra(pa, snr = 100, rmp = T)
#   pa <- averageInterFragSpectra(pa, snr = 100, rmp = T)
#   pa <- averageAllFragSpectra(pa, snr = 100, rmp = F)
#
#   td <- tempdir()
#
#   db_pth = create_database(pa, xset, out_dir = td)
#
#   con <- DBI::dbConnect(RSQLite::SQLite(), file.path(db_pth))
#
#   cpg <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups')
#   expect_equal(nrow(cpg), 375)
#
#   cpgX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_c_peak_group')
#   expect_equal(nrow(cpgX), 780)
#
#   csX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_s_peak_meta')
#   expect_equal(nrow(csX), 270)
#
#   c_peaks <- DBI::dbGetQuery(con, 'SELECT * FROM c_peaks')
#   expect_equal(nrow(c_peaks), 780)
#
#
# })
