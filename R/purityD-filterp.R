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



#' @title Filter out peaks based on intensity and RSD criteria
#'
#' @description
#' Uses a purityD object remove peaks from either (or both) samples and
#' blanks that are either below an intensity threshold
#' or greater than a Relative Standard Deviation (RSD) threshold
#'
#' @aliases filterp
#'
#' @param Object object; purityD object
#' @param thr numeric; intensity threshold
#' @param rsd numeric; rsd threshold
#' @param sampleOnly boolean; if only the sample (not blanks) should be filtered
#' @examples
#'
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#'
#' ppDIMS <- purityD(inDF, cores=1)
#' ppDIMS <- averageSpectra(ppDIMS)
#' ppDIMS <- filterp(ppDIMS, thr = 5000)
#' @return  purityD object
#' @export
setMethod(f="filterp", signature="purityD", definition= function(Object, thr = 5000, rsd=20, sampleOnly = TRUE) {
  if(sampleOnly){
    Object@avPeaks$processed[Object@sampleIdx] <- lapply(Object@avPeaks$processed[Object@sampleIdx],
                                                               function(x){ x[(x$i>thr & x$rsd<rsd), ]})

  }else{
    Object@avPeaks$processed <- lapply(Object@avPeaks$processed, function(x){ x[(x$i>thr & x$rsd<rsd), ]})
  }
  return(Object)
})
