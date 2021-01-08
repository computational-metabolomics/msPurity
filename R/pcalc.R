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



#' @title Perform purity calculation on a peak matrix
#'
#' @description
#' This is the main purity calculation that is performed in purityX,
#' purityD and purityA.
#' \itemize{
#'  \item{Takes in a matrix of peaks}
#'  \item{gets isolation window based on mzmin mzmax}
#'  \item{locates the mz target in the peak matrix}
#'  \item{removes isotopic peaks}
#'  \item{removes any peaks below limit (percentage of target peak intensity)}
#'  \item{normalises}
#'  \item{Calculates purity: Divides the target peak intensity by the total peak intensity for
#'     the isolation window}
#' }
#'
#' @param peaks matrix; Matrix of peaks consisting of 2 columns: mz and i
#' @param mzmin numeric; Isolation window (min)
#' @param mzmax numeric; Isolation window (max)
#' @param mztarget numeric; The mz window to target in the isolation window
#' @param ppm numeric; PPM tolerance for the target mz value. If NA will presume targetMinMZ and targetMaxMZ will be used
#' @param targetMinMZ numeric; Range to look for the mztarget (min)
#' @param targetMaxMZ numeric; Range to look for the mztarget (max)
#' @param iwNorm boolean; If TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5% (0.05)
#' @param isotopes boolean; TRUE if isotopes are to be removed
#' @param im matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#'
#' @return a vector of the purity score and the number of peaks in the window e.g c(purity, pknm)
#'
#' @examples
#' pm <- rbind(c(100, 1000),c(101.003, 10))
#' pcalc(pm, mzmin = 98, mzmax = 102, mztarget=100, ppm=5)
#' pcalc(pm, mzmin = 98, mzmax = 102, mztarget=100, ppm=5, isotopes = TRUE)
#'
#' @export
pcalc <- function(peaks, mzmin, mzmax, mztarget, ppm=NA, iwNorm=FALSE,
                  iwNormFun=NULL, ilim=0, targetMinMZ=NA, targetMaxMZ=NA,
                  isotopes=FALSE, im=NULL){
  if(is.null(peaks)){
    return(c(NA, NA))
  }else if(nrow(peaks)==0){
    return(c(NA, NA))
  }
  if(is.na(mztarget)){
    return(c(NA, NA))
  }

  # Make sure the input is a matrix
  if(is.vector(peaks)){
    peaks <- matrix(peaks, ncol = length(peaks), nrow = 1)
  }else if(is.data.frame(peaks)){
    peaks <- peaks[,c("mz","i")]
    peaks <- data.matrix(peaks)
  }

  # Get the isolation window range
  subp <- peaks[(peaks[,1]>mzmin) & (peaks[,1]<mzmax),]

  if(is.vector(subp)){
    subp <- matrix(subp, ncol = length(subp), nrow = 1)
  }

  # Check if the mztarget is to be identified based on a ppm tolerance or
  # pre deterimend range
  if(!is.na(ppm)){
    # check if the target mz in the list
    mzLo <- round(mztarget - ((mztarget*0.000001)*ppm),10)
    mzUp <- round(mztarget + ((mztarget*0.000001)*ppm),10)
    mtch <- subp[(subp[,1]>mzLo) & (subp[,1]<mzUp),]
    cmatch <- cleanUp(mtch, mztarget)

  }else if(!is.na(targetMinMZ) || !is.na(targetMaxMZ)){
    # pre determined range
    mtch <- subp[(subp[,1]>=targetMinMZ) & (subp[,1]<=targetMaxMZ), ]
    cmatch <- cleanUp(mtch, mztarget)

  }else{
    # As no ppm or target range given we presume the mz to have an
    # exact match within the subp peak list
    cmatch <- unlist(subp[subp[,1]==unlist(mztarget),])
  }



  mtchmz <- unname(cmatch[1])
  mtchi <- unname(cmatch[2])



  # if nothing was matched then the purity will be zero
  if(identical(mtchi, numeric(0))){
    return(c(NA, NA))
  }else if(is.na(mtchi)){
    return(c(NA, NA))
  }else if (mtchi == 0){
    return(c(0, NA))
  }



  if(isotopes){
    if(is.null(im)){
      # matrix of isotopes to look at. by default this is just C12/C13 (only single, double and triple charge)
      # row composed of:
      # 'isotope_id','mass diff', 'abundance of isotope', 'ppm tol for mz', 'abundance buffer',
      # 'charge', 'relative atomic mass (int)', 'xflag'
      # Notes: The estimated max and min number of atoms in is calculated from the relative atomic mass. This
      #        a very broad estimate as such the number should really just be a int
      #        The xflag indicates if the larger (mass) isotope is the most abundant or less abundant. e.g.
      #        for c12 to c13, the c13 is less abundant so we flag as 1
      #        for Li6 to Li7, the Li7 is more abundant so we would flag as 0
      im = rbind(c(1, 1.003355, 1.07, 4, 0.1, 1, 12, 1), # C13
                 c(2, 1.003355/2, 1.07, 4, 0.1, 2, 12, 1), # C13 double charge
                 c(3, 1.003355/3, 1.07, 4, 0.1, 3, 12, 1))  # C13 triple charge

      # example below is if the C14 isotope is to be looked for as well. Not relevant in most cases though
#       im = rbind(c(1, 1.003355, 1.07, 4, 0.1, 1, 12, 1), # C13
#                  c(2, 1.003355/2, 1.07, 4, 0.1, 2, 12, 1), # C13 double charge
#                  c(3, 1.003355/3, 1.07, 4, 0.1, 3, 12, 1),  # C13 triple charge
#                  c(4, 2.003242, 0.1, 4, 0.1, 1, 12, 1), # C14
#                  c(5, 2.003242/2, 0.1, 4, 0.1, 2, 12, 1), # C14 double charge
#                  c(6, 2.003242/3, 0.1, 4, 0.1, 3, 12, 1)) # C14 triple charge

    }

    subp <- removeIsotopes(subp, im, mtchmz, mtchi)
  }

  # Normalise based on the peak shape of ion contribution
  if(iwNorm){
    # If the isolation window is to be normalised based on the contribution
    # of ions in relation to position in isolation window
    middle <- mzmax-(mzmax-mzmin)/2

    # Get the mz position as distances away from the centre position (0 Da)
    adjustmz = subp[,1]-middle
    subp[,2] = subp[,2]*iwNormFun(adjustmz)

    # get the adjusted intensity for the matched intensity
    # (bit of duplication here, but easiest way to update code for isotopes without
    # having to change too much)
    mtchi <- mtchi*iwNormFun(mtchmz-middle)
  }

  # only count intensity from peaks that are above x percentage of the target
  ilimit <- mtchi*ilim
  subpl <- subp[subp[,2]>ilimit,]


  # Get the updated total intensity and peak number
  if((!is.na(mtchi)) && ((is.null(subpl)) || (length(subpl)==0))){
    return(c(1, 1))
  }else if(is.vector(subpl)){
    li <- subpl[2]
    pknm = 1
  }else{
    li <- sum(subpl[,2])
    # peak count
    pknm = nrow(subpl)
  }

  # Perform the purity calculation
  purity <- unlist(mtchi)/unlist(li)

  return(c(purity, pknm))
}

cleanUp <- function(mtch, mztarget){
  if(is.null(mtch)){
    mtch <- c(NA, NA)
  }else if(is.vector(mtch)){
    mtch <- unlist(mtch)
  }else if(nrow(mtch)>1){
    # multiple mz's in the same region of ppm tolerance.
    #Get the closest match for mz and use that intensity
    mtch <- mtch[which.min(abs(mtch[,1] - mztarget)),]
  }else if(nrow(mtch)==1) {
    # just one match so get the single intensity value
    mtch <- unlist(mtch)
  }else{
    mtch <- c(NA,NA)
  }

  return(mtch)

}


removeIsotopes <- function(peaks, im, target_mz, target_i, writeout=FALSE){
  # loop through all the isotopes in the im list  tomn

  if (!is.data.frame(peaks)){
    peaks <- data.frame(peaks)
    colnames(peaks) <- c('mz', 'i')
  }

  peaks$id <- seq(1, nrow(peaks))

  l = apply(im, 1, im_tag, peaks=peaks, target_mz=target_mz, target_i=target_i)
  idf <- data.frame(matrix(unlist(l), nrow=length(l)*2, byrow=TRUE))

  to_remove <- idf[,2][!is.na(idf[,2])]

  if (identical(numeric(0), to_remove)){
  }else{
    if (writeout){
      rst <- MHmakeRandomString()
      colnames(idf) <- c('iid', 'peakid', 'charge')
      pref <- paste('/tmp/', rst, round(target_mz,4), round(target_i, 4), sep="_")
      write.csv(idf, paste(pref, 'idf.csv', sep="_"), row.names = FALSE)
      write.csv(peaks, paste(pref, 'peaks.csv', sep="_"), row.names = FALSE)
    }
  }

  updated_peaks <- as.matrix(peaks[!peaks$id %in% to_remove,])[,1:2]

  if(is.vector(updated_peaks)){
    updated_peaks <- matrix(updated_peaks, nrow=1, ncol=2)
  }

  return(updated_peaks)


}

MHmakeRandomString <- function(n=1, lenght=12){
  # taken from https://ryouready.wordpress.com/2008/12/18/generate-random-string-name/
  # used for testing only
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    lenght, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}



im_tag <- function(x, peaks, target_mz, target_i){

  iid <- x[1]
  mdiff <- x[2]
  adiff <- x[3]
  mtol <- x[4]
  atol <- x[5]
  charge <- x[6]
  ram <- x[7]
  xflag <- x[8]

  # get the difference between the target and all the other peaks
  mr = target_mz+mdiff
  ml = target_mz-mdiff

  peaks$ppmR <- sapply(peaks$mz, ppm_error, MZiso=mr)
  peaks$ppmL <- sapply(peaks$mz, ppm_error, MZiso=ml)

  # Check if contaminating peaks are the M+1.. isotopes
  # filter the peaks so we are only looking at those with less than the
  # determined ppm tolerance
  peakr <- peaks[which(peaks$ppmR<mtol),]

  # Get intensity ratio tolerances
  if(xflag==1){
    i_range_r <- get_iso_intensity_range(rp=FALSE, target_mz, target_i, adiff, ram, atol)
  }else{
    i_range_r <- get_iso_intensity_range(rp=TRUE, target_mz, target_i, adiff, ram, atol)
  }

  inten_min_r <- i_range_r[1]
  inten_max_r <- i_range_r[2]

  peakr <- peakr[peakr$i<=inten_max_r & peakr$i>=inten_min_r,]

  peakr <- peakr[order(peakr$ppmR, decreasing = FALSE),]

  # Check best match
  bmatchr <- peakr[1,]

  # Check if contaminating peaks are the more abundant isotope (e.g. C12)
  # filter the peaks so we are only looking at those with less than the
  # determined ppm tolerance
  peakl <- peaks[which(peaks$ppmL<mtol),]

  if(xflag==1){
    i_range_l <- get_iso_intensity_range(rp=TRUE, target_mz, target_i, adiff, ram, atol)
  }else{
    i_range_l <- get_iso_intensity_range(rp=FALSE, target_mz, target_i, adiff, ram, atol)
  }

  inten_min_l <- i_range_l[1]
  inten_max_l <- i_range_l[2]

  peakl <- peakl[peakl$i<=inten_max_l & peakl$i>=inten_min_l,]
  peakl <- peakl[order(peakl$ppmL, decreasing = FALSE),]

  bmatchl <- peakl[1,]

  ilog <- list(c(iid, bmatchr$id, charge),
               c(iid, bmatchl$id, charge-1))
  return(ilog)

}

get_iso_intensity_range <- function(rp=FALSE, target_mz, target_i, adiff, ram, atol){
  numE  <- abs(round(target_mz / ram)) # max. number of element  (e.g. C) in molecule
  if(rp){
    inten_max <- (target_i / (1 * adiff - atol)) * 100 # lowest possible intensity (for most abundant isotope e.g C12)
    inten_min <- (target_i / (numE * adiff + atol)) * 100 # highest possible intensity (for most abundant isotope e.g C12)
  }else{
    inten_max <- ((numE * adiff + atol) / 100) * target_i # highest possible intensity (for M+1(or more)  isotope)
    inten_min <- ((1 * adiff - atol) / 100) * target_i # lowest possible intensity (for M+1(or more) isotope)
  }
  return(c(inten_min, inten_max))
}

ppm_error <- function(MZcont, MZiso){
  # MZcont = Observered MZcontaminating peak
  # MZiso = Theoritical MZisotopic peak
  abs(1e6*(MZcont-abs(MZiso))/MZiso)
}
