#' @title Validate precursor purity predictions using LC-MS and LC-MS/MS dataset
#'
#' @description
#' The method is used to validate the precursor purity predictions made
#' from an LC-MS dataset
#'
#' @aliases validate
#'
#' @param pa = purityA object
#' @param ppLCMS = purityX object
#' @return purityA object
#' @export
setMethod(f="validate", signature="purityA",  definition= function(pa, ppLCMS){
  check <- merge(pa@grped_df, ppLCMS@predictions, by = "grpid", sort = FALSE, all.x = TRUE)
  # count how many ms/ms spectra there are per feature
  check$grpid <- as.numeric(check$grpid)

  #msmsperfeature <- median(daply(check,.(grpid), function(x){nrow(x)}))

  check <- check[order(check$grpid),]

  pdiff <- abs(check$median-check$inPurity)

  predictP <- check$median

  head(pa@grped_df)
  pa@grped_df <- cbind(pa@grped_df, predictP, pdiff)

  return(pa)
})

# \dontshow{
# datapth <- system.file("extdata", "lcms", package="msPurityData")
# files <- c("LCMS_1.mzML","LCMS_2.mzML", "LCMSMS_1.mzML", "LCMSMS_2.mzML")
# mzdatafiles <- vector()
# for (i in 1:length(files)){
#    mzdatafiles[i] <- file.path(datapth, "mzML", files[i])
# }
# xset <- readRDS(file.path(datapth, "example_results", "xset.rds"))
# xset@filepaths <- mzdatafiles
# }
# pa <- purityA(mzdatafiles, cores=1, interpol = "linear")
# pa <- frag4feature(pa, xset)
# ppLCMS <- purityPL(xset, fileignore = c(3,4), cores = 1, xgroups = 6)
# pa <- validate(pa, ppLCMS)


