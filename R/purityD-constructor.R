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
#' @rdname msPurity-deprecated
#' @section initialize (purityD method):
#' The class used to predict purity from an DI-MS dataset.
#'
#' \strong{Value:} purityD object. The files and their sample information are also kept as an \code{MsExperiment} in
#' the experiment slot, and the averaged peak lists of later steps as \code{Spectra} in the avSpectra slot
#' (see \code{\link{averagedPeaks}})
#'
#' @param .Object object; purityD object
#' @param cores numeric; Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 1
#' @param mzRback character; deprecated and ignored. Raw data is read through Spectra, which uses the pwiz reader
#' @export purityD
setMethod("initialize", "purityD", function(.Object, fileList, cores=1, mzML=TRUE, mzRback='pwiz'){
  .msp_deprecate_dims("purityD")
  if (missing(fileList)){
    return(.Object)
  }

  # An MsExperiment gives the files and their sample information
  if (is(fileList, "MsExperiment")){
    .Object@experiment <- fileList
    fileList <- .pd_filelist(fileList)
  }else{
    .Object@experiment <- .pd_experiment(fileList)
  }

  .Object@fileList <- fileList
  .msp_deprecate_mzRback(mzRback)

  for (i in 1:nrow(fileList)){
    file <- as.character(fileList$files[i])

  }

  .Object@sampleIdx <- as.numeric(rownames(fileList[fileList$sampleType=="sample",]))
  .Object@cores <- cores
  .Object@mzML <- mzML
  .Object@mzRback <- 'pwiz'

  return(.Object)
})
