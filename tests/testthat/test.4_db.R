context ("checking database functions")

test_that("checking createDatabase functions (new schema)", {
  print ("\n")
  print("#######################################################")
  print("## Checking database (new schema)                    ##")
  print("#######################################################")

  pa <- readRDS(system.file("extdata", "tests", "purityA", "9_averageAllFragSpectra_with_filter_pa.rds", package="msPurity"))

  msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
  pa@fileList[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
  pa@fileList[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]

  fns = c("msms_only_xcmsnexp.rds", "msms_only_xset.rds")
  for (fn in fns){
    xcmsObj <- readRDS(system.file("extdata","tests", "xcms", fn, package="msPurity"))
    if('XCMSnExp' == class(xcmsObj)){
      xcmsObj@phenoData@data[1,] = msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
      xcmsObj@phenoData@data[2,] = msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
      xcmsObj@processingData@files[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
      xcmsObj@processingData@files[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
    }else if('xcmsSet' == class(xcmsObj)){
      xcmsObj@filepaths[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
      xcmsObj@filepaths[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
    }

    td <- tempdir()
    db_pth = createDatabase(pa = pa, xcmsObj = xcmsObj, outDir = td)

    con <- DBI::dbConnect(RSQLite::SQLite(), file.path(db_pth))

    cpg <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups')
    expect_equal(nrow(cpg), 521)

    cpgX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_c_peak_group')
    expect_equal(nrow(cpgX), 799)

    csX <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_X_s_peak_meta')
    expect_equal(nrow(csX), 75)

    c_peaks <- DBI::dbGetQuery(con, 'SELECT * FROM c_peaks')
    expect_equal(nrow(c_peaks), 800)

    s_peak_meta <- DBI::dbGetQuery(con, 'SELECT * FROM s_peak_meta')
    expect_equal(nrow(s_peak_meta), 1747)

  }

})

