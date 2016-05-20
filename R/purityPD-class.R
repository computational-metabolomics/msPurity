######################################################################
# Constructor
######################################################################
#' @title An S4 class to represent a DI-MS purityPD
#'
#' @description
#' The class used to predict purity from an DI-MS dataset.
#'
#' @param .Object object = purityPD object
#' @param fileList data.frame = created using GetFiles, data.frame with filepaths and sample class information
#' @param cores numeric = Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 1
#' @param mzML boolean = TRUE if mzML to be used FALSE if .csv file to be used
#' @examples
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityPD(fileList=inDF, cores=1, mzML=TRUE)
#' @return purityPD object
#' @export purityPD
purityPD <- setClass(
  # Set the name for the class
  "purityPD",

  # Define the slots
  slots = c(
    fileList = "data.frame",
    sampleTypes = "vector",
    avPeaks  = "list",
    avParam = "list",
    sampleIdx = "vector",
    cores = "numeric",
    purityParam = "list",
    outFiles = "vector",
    mzML = "logical"
  )
)

######################################################################
# A few getters and Setters
######################################################################
#' @title Get peaklist for a purityPD object
#'
#' @description
#' output peak list for a purityPD object
#'
#' @aliases getP
#'
#' @param x object = purityPD object
#' @return peaks
#' @examples
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityPD(fileList=inDF, cores=1, mzML=TRUE)
#' peaks <- getP(ppDIMS)
#' @export
setMethod("getP", "purityPD", function(x) x@avPeaks)

setMethod("updatePeaks", "purityPD", function(x, newlist) {
  x@avPeaks <- newlist;
  validObject(x);
  return(x)
})


######################################################################
# Show method
######################################################################
#' @title Show method for purityPD
#'
#' @description
#' Show method for purityPD object
#'
#' @param object = purityPD object
#' @return a print statement of regarding object
#' @export
setMethod("show", "purityPD", function(object) {
  print("purityPD object for processing MS data")
})
