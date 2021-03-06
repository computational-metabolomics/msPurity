# msPurity R package for processing MS/MS data - Copyright (C)
#
# This file is part of msPurity.
#
# msPurity is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msPurity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msPurity.  If not, see <https://www.gnu.org/licenses/>.



######################################################################
# Constructor
######################################################################
#' @title An S4 class to represent a DI-MS purityD
#'
#' @description
#' The class used to assess anticipated purity from a DI-MS run
#'
#' @param .Object object; purityD object
#' @param fileList data.frame; Created using GetFiles, data.frame with filepaths and sample class information
#' @param cores numeric; Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 1
#' @param mzML boolean; TRUE if mzML to be used FALSE if .csv file to be used
#' @examples
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
#' @return purityD object
#' @export purityD
purityD <- setClass(
  # Set the name for the class
  "purityD",

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
    mzML = "logical",
    groupedPeaks = "data.frame",
    mzRback = 'character'
  )
)

######################################################################
# A few getters and Setters
######################################################################
#' @title Get peaklist for a purityD object
#'
#' @description
#' output peak list for a purityD object
#'
#' @aliases getP
#'
#' @param x object; purityD object
#' @return peaks
#' @examples
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
#' peaks <- getP(ppDIMS)
#' @export
setMethod("getP", "purityD", function(x) x@avPeaks)

setMethod("updatePeaks", "purityD", function(x, newlist) {
  x@avPeaks <- newlist;
  validObject(x);
  return(x)
})


######################################################################
# Show method
######################################################################
#' @title Show method for purityD
#'
#' @description
#' Show method for purityD object
#'
#' @param object = purityD object
#' @return a print statement of regarding object
#' @export
setMethod("show", "purityD", function(object) {
  print("purityD object for processing MS data")
})
