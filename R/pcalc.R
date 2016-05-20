################################
# pcalc: The purity calculation
################################
# This is the main purity calculation that is performed in purityPL,
# purityPD and purityA.
#  - Takes in a matrix of peaks
#  - gets isolation window based on mzmin mzmax
#  - normalises (if required)
#  - locates the mz target in the peak matrix
#  - removes any peaks below limit (percentage of target peak intensity)
#  - Divides the target peak intensity by the total peak intensity for
#      the isolation window
pcalc <- function(peaks, mzmin, mzmax, mztarget, ppm=NA, iwNorm=FALSE,
                  iwNormFun=NULL, ilim=0, targetMinMZ=NA, targetMaxMZ=NA){

  # Make sure the input is a matrix
  if(is.vector(peaks)){
    peaks <- matrix(peaks, ncol = length(peaks), nrow = 1)
  }else if(is.data.frame(peaks)){
    peaks <- data.matrix(peaks)
  }


  # Get the isolation window range
  subp <- peaks[(peaks[,1]>mzmin) & (peaks[,1]<mzmax),]

  if(is.vector(subp)){
    subp <- matrix(subp, ncol = length(subp), nrow = 1)
  }

  # Normalise based on the peak shape of ion contribution
  if(iwNorm){
    # If the isolation window is to be normalised based on the contribution
    # of ions in relation to position in isolation window
    adjustmz = (subp[,1]-mzmin)/(mzmax-mzmin)
    subp[,2] = subp[,2]*iwNormFun(adjustmz)
  }

  # Check if the mztarget is to be identified based on a ppm tolerance or
  # pre deterimed range
  if(!is.na(ppm)){
    # check if the target mz in the list
    mzLo <- round(mztarget - ((mztarget*0.000001)*ppm),10)
    mzUp <- round(mztarget + ((mztarget*0.000001)*ppm),10)
    mtch <- subp[(subp[,1]>mzLo) & (subp[,1]<mzUp),]
    mtchi <- cleanUp(mtch, mztarget)

  }else if(!is.na(targetMinMZ) || !is.na(targetMaxMZ)){
    # pre determined range
    mtch <- subp[(subp[,1]>=targetMinMZ) & (subp[,1]<=targetMaxMZ), ]
    mtchi <- cleanUp(mtch, mztarget)

  }else{
    # As no ppm or target range given we presume the mz to have an
    # exact match within the subp peak list
    mtchi <- unlist(subp[subp[,1]==unlist(mztarget),2])
  }

  # if nothing was matched then the purity will be zero
  if(is.na(mtchi)){
    return(c(NA, NA))
  }else if (mtchi == 0){
    return(c(0, NA))
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

  if(is.vector(mtch)){
    mtchi <- unlist(mtch[2])
  }else if(nrow(mtch)>1){
    # multiple mz's in the same region of ppm tolerance.
     #Get the closest match for mz and use that intensity
    mtcht <- mtch[which.min(abs(mtch[,1] - mztarget)),]
    mtchi <- unlist(mtcht[2])
  }else if(nrow(mtch)==1) {
    # just one match so get the single intensity value
    mtchi <- unlist(mtch[2])
  }else{
    mtchi <- NA
  }
  return(mtchi)

}
