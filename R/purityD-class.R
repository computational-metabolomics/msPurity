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
#' @section purityD:
#' The class used to assess anticipated purity from a DI-MS run
#'
#' \strong{Value:} purityD object
#'
#' @param .Object object; purityD object
#' @param cores numeric; Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 1
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
    mzRback = 'character',
    # Averaged peak lists, one spectrum per file and stage.
    avSpectra = "Spectra",
    # The files and their sample information.
    experiment = "MsExperiment"
  ),
  prototype = prototype(
    avSpectra = Spectra::Spectra(Spectra::MsBackendMemory()),
    experiment = MsExperiment::MsExperiment()
  )
)

######################################################################
# A few getters and Setters
######################################################################
#' @rdname msPurity-deprecated
#' @section getP:
#' output peak list for a purityD object
#'
#' \strong{Value:} peaks
#'
#' @aliases getP
#'
#' @param x object; purityD object
#' @export
setMethod("getP", "purityD", function(x) {
  .msp_deprecate_dims("getP")
  averagedPeaks(x)
})

setMethod("updatePeaks", "purityD", function(x, newlist) {
  .msp_deprecate_dims("updatePeaks")
  x <- .pd_update(x)
  x@avPeaks <- newlist;
  x <- .pd_end(x)
  validObject(x);
  return(x)
})


######################################################################
# Show method
######################################################################
#' @rdname msPurity-deprecated
#' @section show (purityD method):
#' Show method for purityD object
#'
#' \strong{Value:} a print statement of regarding object
#'
#' @param object = purityD object
#' @export
setMethod("show", "purityD", function(object) {
  .msp_deprecate_dims("show")
  print("purityD object for processing MS data")
})
