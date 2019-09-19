context ("checking database functions")

test_that("checking create_database (old schema)", {
  print ("\n")
  print("########################################################")
  print("## Checking database (old schema)                     ##")
  print("########################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "9_averageAllFragSpectra_with_filter_pa.rds", package="msPurity"))
  xset <- readRDS(system.file("extdata","tests", "xcms", "msms_only_xset.rds", package="msPurity"))

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  pa@fileList[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
  pa@fileList[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
  xset@filepaths[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
  xset@filepaths[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]


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

  cpg <- DBI::dbGetQuery(con, 'SELECT * FROM av_peaks')
  expect_equal(nrow(cpg), 1886 )

  ####################################
  # Check EIC database from purityX
  ####################################
  px  <- purityX(xset, saveEIC = TRUE, sqlitePth = db_pth, plotP = TRUE, xgroups=c(1,2,3))
  eics <- DBI::dbGetQuery(con, 'SELECT * FROM eics')
  expect_equal(nrow(eics), 211)

})


test_that("checking createDatabase functions (new schema)", {
  print ("\n")
  print("#######################################################")
  print("## Checking database (new schema)                    ##")
  print("#######################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "9_averageAllFragSpectra_with_filter_pa.rds", package="msPurity"))
  xset <- readRDS(system.file("extdata","tests", "xcms", "msms_only_xset.rds", package="msPurity"))

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  pa@fileList[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
  pa@fileList[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
  xset@filepaths[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
  xset@filepaths[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]


  td <- tempdir()
  db_pth = createDatabase(pa, xset, outDir = td)

  con <- DBI::dbConnect(RSQLite::SQLite(), file.path(db_pth))

  cpg <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups')
  expect_equal(nrow(cpg), 375)

  cpgX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_c_peak_group')
  expect_equal(nrow(cpgX), 780)

  csX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_s_peak_meta')
  expect_equal(nrow(csX), 270)

  c_peaks <- DBI::dbGetQuery(con, 'SELECT * FROM c_peaks')
  expect_equal(nrow(c_peaks), 780)

  s_peak_meta <- DBI::dbGetQuery(con, 'SELECT * FROM s_peak_meta')
  expect_equal(nrow(s_peak_meta), 1919)


})

