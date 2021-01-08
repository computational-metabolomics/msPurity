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



#' @title Gaussian normalisation for isolation window efficiency
#'
#' @description
#' Creates a function based on a gaussian curve shape that will normalise any intensity values within
#' a defined isolation window.
#'
#' The function that is created will output a value between 0 to 1 based on the position between
#' the minOff and maxOff params. (The value 1.0 being equivalent to 100% efficient)
#'
#' @param sdlim numerical; Standard deviation limit for gaussian curve
#' @param minOff numerical; Offset to the 'left' for the precursor range. (Should be negative)
#' @param maxOff character; Offset to the 'left' for the precursor range. (Should be positive)
#' @examples
#'
#' iwNormFun <- iwNormGauss(minOff=-0.5, maxOff=0.5)
#' pm <- data.frame(mz=c(99.5, 99.9, 100, 100.1, 100.5),i=c(1000, 1000, 1000, 1000, 1000))
#' mzmax = 100.5
#' mzmin = 99.5
#' middle <- mzmax-(mzmax-mzmin)/2
#' adjustmz = pm$mz-middle
#'
#' # normalise the intensities
#' pm$normi = pm$i*iwNormFun(adjustmz)
#'
#'
#' @return normalisation function for selected range.
#' @export
iwNormGauss <- function(sdlim=3, minOff=-0.5, maxOff=+0.5){

  # get a gaussian curve
  x <- seq(0-sdlim, 0+sdlim, 0.05)
  y <- dnorm(x, mean = 0)

  # linear scaling to range 0 to 1
  y <- (y-min(y))/(max(y)-min(y))

  # scaling to the min and max of the window
  x <- seq(minOff, maxOff, length.out = length(y))

  # Create function that outputs a 0 to 1 value depending on
  # the position of the x (i.e. mz)
  f <- approxfun(x, y)
}



#' @title Q-Exactive +/- 0.5 range, normalisation for isolation window efficiency
#'
#' @description
#' Creates a function based on a previous experimental analysis of a Q-Exactive at +/- 0.5
#' isolation window efficiency. See http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b04358
#'
#' The function that is created will output a value between 0 to 1 based on the position between
#' the minOff and maxOff params
#'
#' NOTE: The resulting function will work for values greater that 0.5 and less than -0.5.
#'
#' This is because (on our instrument tested at least) when using a window of +/- 0.5,
#' the isolation is NOT confined to the +/-0.5 Da window. Resulting in ions from outside the window
#' being isolated. For this reason the function can normalise values outside of the
#' the +/- 1 Da range. Please see above paper figure 3 for more details.
#'
#' @examples
#' iwNormFun <- iwNormQE.5()
#' pm <- data.frame(mz=c(99.5, 99.9, 100, 100.1, 100.5),i=c(1000, 1000, 1000, 1000, 1000))
#' mzmax = 100.5
#' mzmin = 99.5
#' middle <- mzmax-(mzmax-mzmin)/2
#' adjustmz = pm$mz-middle
#'
#' # normalise the intensities
#' pm$normi = pm$i*iwNormFun(adjustmz)
#'
#' @return normalisation function for +/- 0.5 range for Q-Exactive
#' @export
iwNormQE.5 <- function(){
    y <- c(0.0000, 0.0000, 0.0000, 0.0550, 0.2336, 0.4437, 0.6509, 0.8210,
           0.9339, 0.9915, 0.9975, 0.9555, 0.8694, 0.7428, 0.5805, 0.3986,
           0.2208, 0.0710, 0.0000, 0.0000, 0.0000)
    x <- seq(-1, 1, 0.1)
    f <- approxfun(x, y)
}

#' @title Raised cosine normalisation for isolation window efficiency
#'
#' @description
#' Creates a function based on a rasied cosine curve shape that will normalise any intensity values within
#' a defined isolation window
#'
#' The function that is created will output a value between 0 to 1 based on the position between
#' the minOff and maxOff params
#'
#' @param minOff numerical; Offset to the 'left' for the precursor range. (Should be negative)
#' @param maxOff character; Offset to the 'left' for the precursor range. (Should be positive)
#' @examples
#' iwNormFun <- iwNormRcosine()
#' pm <- data.frame(mz=c(99.5, 99.9, 100, 100.1, 100.5),i=c(1000, 1000, 1000, 1000, 1000))
#' mzmax = 100.5
#' mzmin = 99.5
#' middle <- mzmax-(mzmax-mzmin)/2
#' adjustmz = pm$mz-middle
#'
#' # normalise the intensities
#' pm$normi = pm$i*iwNormFun(adjustmz)
#' @return normalisation function for selected range
#' @export
iwNormRcosine <- function(minOff = -0.5, maxOff = +0.5){
   # s <- sapa::taper(type="raised cosine")
   # y <- as.vector(s)
   # Taken from above function from package sapa. Generates a raised cosine curve
   #
   y <- c(3e-04, 0.001, 0.0023, 0.0041, 0.0064, 0.0091, 0.0123, 0.016, 0.02, 0.0244, 0.0291,
          0.0341, 0.0393, 0.0448, 0.0504, 0.0562, 0.062, 0.0679, 0.0738, 0.0796, 0.0853, 0.0909,
          0.0962, 0.1014, 0.1062, 0.1108, 0.115, 0.1188, 0.1222, 0.1252, 0.1277, 0.1297, 0.1313, 0.1323,
          0.1328, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329,
          0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329,
          0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1329, 0.1328, 0.1323, 0.1313, 0.1297, 0.1277,
          0.1252, 0.1222, 0.1188, 0.115, 0.1108, 0.1062, 0.1014, 0.0962, 0.0909, 0.0853, 0.0796, 0.0738,
          0.0679, 0.062, 0.0562, 0.0504, 0.0448, 0.0393, 0.0341, 0.0291, 0.0244, 0.02, 0.016, 0.0123, 0.0091,
          0.0064, 0.0041, 0.0023, 0.001, 3e-04)

   x <- seq(minOff, maxOff, length.out = length(y))

   # linear scaling to range 0 to 1
   y <- (y-min(y))/(max(y)-min(y))
   f <- approxfun(x, y)

}




