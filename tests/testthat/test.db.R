test_that("checking database functions", {
  print("########################################################")
  print("## Checking LCMS based class and functions            ##")
  print("########################################################")

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  library(xcms)

  xset <- xcmsSet(msmsPths)
  xset <- group(xset)

  pa <- purityA(msmsPths)
  pa <- frag4feature(pa, xset)

  td <- tempdir()
  db_pth = create_database(pa, xset, out_dir = td)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(db_pth))

  cpg <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups')
  expect_equal(nrow(cpg), 375)

  cpgX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_c_peak_group')
  expect_equal(nrow(cpgX), 780)

  csX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_s_peak_meta')
  expect_equal(nrow(csX), 270)

  c_peaks <- DBI::dbGetQuery(con, 'SELECT * FROM c_peaks')
  expect_equal(nrow(c_peaks), 780)

  ####################################
  # Check EIC database from purityX
  ####################################
  msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MS_")
  xset <- xcmsSet(msPths, method = "centWave")
  xset <- group(xset)

  px  <- purityX(xset, saveEIC = TRUE, sqlitePth = db_pth, plotP = TRUE, xgroups=c(1,2,3))



})
