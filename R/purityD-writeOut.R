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



#' @title Using purityD object, save peaks as text files
#'
#' @description
#' Uses a purityD object with references to multiple MS files.
#' Predicts the purity of the processed sample files
#'
#' @aliases writeOut
#'
#' @param Object object; purityD object
#' @param outDir character; Directory to save text files
#' @param original boolean; If the original (unprocessed) files are to be saved to text files
#' @return  purityD object
#' @export
setMethod(f="writeOut", signature="purityD",
          definition= function(Object, outDir, original) {

  dir.create(outDir)

  # write out processed spectra
  for (i in 1:length(Object@avPeaks$processed)){
    # get file info
    fileinfo <- Object@fileList[i,]
    name <- as.character(fileinfo$name)

    # Filename
    fn1 <- 'processed.csv'
    fn2 <- 'original.csv'
    nopeaks <- 'nopeaks.txt'

    fname1 <- file.path(outDir, paste(name, fn1, sep="_"))
    Object@outFiles[i] <- fname1

    if(nrow(Object@avPeaks$processed[[i]])==0){
      fnopeak <- file.path(outDir, paste(name, sep="_NO_PEAKS.csv"))
      write("No peaks after processing, possibly the parameters are to harsh for RSD or min scan",fnopeak)
    }else{
      write.csv(Object@avPeaks$processed[[i]], fname1, row.names = FALSE)
    }

    if(original){
      fname2 <- file.path(outDir, paste(name, fn2, sep="_"))
      write.csv(Object@avPeaks$orig[[i]], fname2, row.names = FALSE)
    }


  }

  # write out meta information
  return(Object)
})
