#' @include purityD-class.R
NULL

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
#' @title Constructor for S4 class to represent a DI-MS purityD
#'
#' @description
#' The class used to predict purity from an DI-MS dataset.
#' @param .Object object; purityD object
#' @param fileList data.frame; created using GetFiles, data.frame with filepaths and sample class information
#' @param cores numeric; Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 1
#' @param mzML boolean; TRUE if mzML to be used FALSE if .csv file to be used
#' @param mzRback character; backend to use for mzR parsing
#' @return purityD object
#' @examples
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityD(fileList=inDF, cores=1, mzML=TRUE)
#' @export purityD
setMethod("initialize", "purityD", function(.Object, fileList, cores=1, mzML=TRUE, mzRback='pwiz'){

  .Object@fileList <- fileList

  for (i in 1:nrow(fileList)){
    file <- as.character(fileList$files[i])

  }

  .Object@sampleIdx <- as.numeric(rownames(fileList[fileList$sampleType=="sample",]))
  .Object@cores <- cores
  .Object@mzML <- mzML
  .Object@mzRback <- 'pwiz'

  return(.Object)
})
