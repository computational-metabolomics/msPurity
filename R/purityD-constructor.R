#' @include purityD-class.R
NULL


######################################################################
# Constructor
######################################################################
#' @title Constructor for S4 class to represent a DI-MS purityD
#'
#' @description
#' The class used to predict purity from an DI-MS dataset.
#' @param .Object object = purityD object
#' @param fileList data.frame = created using GetFiles, data.frame with filepaths and sample class information
#' @param cores numeric = Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 1
#' @param mzML boolean = TRUE if mzML to be used FALSE if .csv file to be used
#' @return purityD object
#' @examples
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
#' @export purityD
setMethod("initialize", "purityD", function(.Object, fileList, cores=1, mzML=TRUE){

  .Object@fileList <- fileList

  for (i in 1:nrow(fileList)){
    file <- as.character(fileList$files[i])

  }

  .Object@sampleIdx <- as.numeric(rownames(fileList[fileList$sampleType=="sample",]))
  .Object@cores <- cores
  .Object@mzML <- mzML

  return(.Object)
})
