#' @import plyr
#' @import mzR
#' @import foreach
#' @import parallel
#' @import doSNOW
#' @import stringr
#' @import reshape2
#' @import fastcluster

#' @title Assessing anticipated purity of XCMS features from an LC-MS run
#'
#' @description
#' Constructor for the purityX class.
#'
#' Given an XCMS object get the anticipated precursor purity of the grouped peaks
#'
#' @param xset object = xcms object
#' @param cores numeric = number of cores to use
#' @param purityType character = Area and average used for the purity predictions. Options are
#'                               "purityFWHMmedian", "purityFWmedian", "purityFWHMmean", "purityFWmean"
#' @param offsets vector = vector of the isolation window upper and lower offsets
#' @param fileignore vector = vector of files to ignore for the prediction calculation
#' @param xgroups vector = vector of xcms groups to perform prediction on
#' @param iwNorm boolean = if TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function = A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric = All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5\% (0.05)
#' @param plotP boolean = TRUE if plot of the EIC of feature and associated contamination is the be save to the working directory
#' @param mzRback character = backend to use for mzR parsing
#' @param isotopes boolean = TRUE if isotopes are to be removed
#' @param im matrix = Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#' @param singleFile numeric = If just a single file for purity is to be calculated (rather than the grouped XCMS peaks). Uses the index of the files in xcmsSet object. If zero this is ignored.
#'
#' @return a purityX object containing a dataframe of predicted purity scores
#' @examples
#' msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
#' xset <- xcms::xcmsSet(msPths)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#' ppLCMS <- purityX(xset, cores = 1, xgroups = c(1, 2))
#'
#' @export
purityX <- function(xset, purityType="purityFWHMmedian", offsets=c(0.5, 0.5),
                    fileignore=NULL, cores=1, xgroups=NULL,
                    iwNorm=FALSE, iwNormFun=NULL, ilim=0.05, plotP=FALSE, mzRback='pwiz', isotopes=FALSE, im=NULL,
                    singleFile=0){

  if (singleFile>0){

    ppLCMS <- xcmsSinglePurity(xset, fileidx=singleFile, offsets=offsets, iwNorm=iwNorm,
                               iwNormFun=iwNormFun, ilim=ilim, plotP=plotP, mzRback=mzRback,
                               isotopes=isotopes, im=im)
  }else{
    ppLCMS <- xcmsGroupPurity(xset, offsets=offsets, fileignore=fileignore, purityType=purityType,
                              cores=cores, xgroups=xgroups, iwNorm=iwNorm, iwNormFun=iwNormFun,
                              ilim=ilim, plotP=plotP, mzRback=mzRback, isotopes=isotopes, im=im)
  }


  return(ppLCMS)

}

xcmsSinglePurity <- function(xset, fileidx, offsets, iwNorm, iwNormFun, ilim, plotP, mzRback, isotopes, im){

  # get the filepaths from the xcms object
  filepth <- xset@filepaths[fileidx]

  # get xcms peaklist
  peaklist <- xset@peaks

  peaklist <- peaklist[peaklist[,'sample'] == fileidx,]

  scanpeaks <- list()
  scanpeaks[[fileidx]] <- getscans(filepth, mzRback)


  # Create a purityX object
  ppLCMS <- new("purityX")


  maxscan <- minscan <- rep(NA, nrow(peaklist))
  peaklist <- cbind(peaklist, minscan, maxscan)

  for(i in 1:nrow(peaklist)){
    peak <- peaklist[i,]
    sid <- as.numeric(peak['sample'])
    raw <- xset@rt$raw[[sid]]
    corrected <- xset@rt$corrected[[sid]]

    rtmed <- as.numeric(peak['rt'])
    rtmin <- as.numeric(peak['rtmin'])
    rtmax <- as.numeric(peak['rtmax'])
    rtmedidx <- which(corrected==rtmed)
    rtminidx <- which(corrected==rtmin)
    rtmaxidx <- which(corrected==rtmax)

    peaklist[i,]['minscan'] <- rtminidx
    peaklist[i,]['maxscan'] <- rtmaxidx

  }
  dfp <- data.frame(peaklist)
  dfp$id <- seq(1, nrow(peaklist))

  if(plotP){
    dir.create(file.path(getwd(), "purityXplots"), showWarnings = FALSE)
  }

  sgrp <- plyr::ddply(dfp, ~ id, pp4file, scanpeaks,
                      rtmed=NA, offsets=offsets, iwNorm=iwNorm, iwNormFun=iwNormFun,
                      ilim=ilim, plotP=plotP, isotopes=isotopes, im=im, singleCheck=FALSE)

  ppLCMS@predictions <- sgrp


  return(ppLCMS)

}



xcmsGroupPurity <- function(xset, purityType, offsets,
                            fileignore, cores, xgroups,
                            iwNorm, iwNormFun, ilim, plotP, mzRback, isotopes, im){
  # get the filepaths from the xcms object
  filepths <- xset@filepaths

  # get xcms peaklist
  peaklist <- xset@peaks

  # Create a purityX object
  ppLCMS <- new("purityX")

  # assign a blank grpid column
  grpid <- rep(0, nrow(peaklist))
  peaklist <- cbind(peaklist, grpid)

  # Add the associated grpid to each peak
  for(i in 1:length(xset@groupidx)){
    if(is.vector(peaklist[xset@groupidx[[i]],])){
      peaklist[xset@groupidx[[i]],]['grpid'] <- i
    }else{
      peaklist[xset@groupidx[[i]],][,'grpid'] <- i
    }
  }

  # remove peaks that have not been grouped together (i.e. only found in 1 file)
  grouplist <- peaklist[peaklist[,'grpid']>0,]
  
  # remove peaks that do not have SN value (means that they will be created from 
  #  the 'fillpeaks' function
  grouplist <- grouplist[!is.na(grouplist[,'grpid']),]

  # Remove files that we do not want to look at
  if(!is.null(fileignore)){
    grouplist <- grouplist[!grouplist[,'sample'] %in% fileignore, ]
  }


  # Select which groups to perform the predictions on
  if(!is.null(xgroups)){
    grouplist <- grouplist[grouplist[,'grpid'] %in% xgroups, ]
  }

  if(is.vector(grouplist)){
    grouplist <- matrix(grouplist, nrow = 1)
    colnames(grouplist) <- colnames(peaklist)
  }

  #print(grouplist)

  # create blank columns for additional useful group info
  rtmedscan <- maxscan <- minscan <- rtmaxraw <- rtminraw <- rtraw <- rep(NA, nrow(grouplist))
  grouplist <- cbind(grouplist, rtraw, rtminraw, rtmaxraw, minscan, maxscan, rtmedscan)

  # Need to get the raw retention time and scans
  # (i.e. the times prior to retention time correction)
  for(i in 1:nrow(grouplist)){
    peak <- grouplist[i,]
    sid <- as.numeric(peak['sample'])
    raw <- xset@rt$raw[[sid]]
    corrected <- xset@rt$corrected[[sid]]

    rtmed <- as.numeric(peak['rt'])
    rtmin <- as.numeric(peak['rtmin'])
    rtmax <- as.numeric(peak['rtmax'])
    rtmedidx <- which(corrected==rtmed)
    rtminidx <- which(corrected==rtmin)
    rtmaxidx <- which(corrected==rtmax)

    grouplist[i,]['rtraw'] <- raw[rtmedidx]
    grouplist[i,]['rtminraw'] <- raw[rtminidx]
    grouplist[i,]['rtmaxraw'] <- raw[rtmaxidx]
    grouplist[i,]['rtmedscan'] <- rtmedidx
    grouplist[i,]['minscan'] <- rtminidx
    grouplist[i,]['maxscan'] <- rtmaxidx

  }

  # Trn into dataframe for ease of use with plyr
  grouplist <- data.frame(grouplist)

  # get in order (for visual checking)
  grouplist <- grouplist[order(grouplist$grpid),]

  # get all the peaks from scans
  scanpeaks <- getscans(filepths, mzRback)

  # Check if it is going to be multi-core
  if(cores<=1){
    para = FALSE
  }else{
    cl<-parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    para = TRUE
  }

  # if iwNorm is TRUE and iwNormFun is NULL then a gaussian model of the
  # isolation window will be used to normalise intensity
  if(is.null(iwNormFun)){
    # Using a gaussian curve 3 SD either side
    iwNormFun <- iwNormGauss(minOff=-offsets[1], maxOff=offsets[2])
  }

  if(plotP){
    dir.create(file.path(getwd(), "purityXplots"), showWarnings = FALSE)
  }


  # perform predictions
  purityPredictions <- plyr::dlply(grouplist,
                                   ~ grpid,
                                   .parallel = para,
                                   predictPurityLCMSloop,
                                   average="median",
                                   scanpeaks=scanpeaks,
                                   offsets=offsets,
                                   iwNorm=iwNorm,
                                   iwNormFun=iwNormFun,
                                   ilim=ilim,
                                   plotP=plotP,
                                   isotopes=isotopes,
                                   im=im)

  if(cores>1){
    parallel::stopCluster(cl)
  }

  # Extract the choosen metric
  dataout <- plyr::ldply(purityPredictions, function(x){ x[,purityType] })

  if(purityType=="purityFWHMmedian"){

    pknmout <- plyr::ldply(purityPredictions, function(x){ x[,"pknmFWHMmedian"] })
  }else{
    pknmout <- plyr::ldply(purityPredictions, function(x){ x[,"pknmFWmedian"] })
  }

  dataout$pknm <- pknmout$median

  # get the median intensity for each group, useful for comparisons
  xcmsgrpi <- xcms::groupval(xset, value="into")
  xcmsgrpmz <- xset@groups[,'mzmed']

  if(!is.null(xgroups)){
    xcmsgrpi <- xcmsgrpi[xgroups,]
    xcmsgrpmz <- xcmsgrpmz[xgroups]
  }

  if((is.null(xgroups)) || (length(xgroups)>1)){
    dataout$i <- apply(xcmsgrpi, 1, median, na.rm = TRUE)
  }else{
    dataout$i <- median(xcmsgrpi, na.rm = TRUE)
  }

  dataout$mz <- xcmsgrpmz

  # save to object
  ppLCMS@predictions <- dataout
  ppLCMS@purityType <- purityType
  ppLCMS@cores <- cores
  if(is.null(fileignore)){
    ppLCMS@fileignore <- NA
  }else{
    ppLCMS@fileignore <- fileignore
  }

  return(ppLCMS)


}




predictPurityLCMSloop <- function(grp, average="median", scanpeaks,
                                  offsets, iwNorm, iwNormFun, ilim, plotP, isotopes, im){

  # Need to loop through for each file in each group
  grp <- data.frame(grp)
  rtmed <- median(grp$rt)

  sgrp <- plyr::ddply(grp, ~ sample, pp4file, scanpeaks,
                      rtmed, offsets, iwNorm, iwNormFun, ilim, plotP, isotopes, im)

  puritySummary <- apply(sgrp[,2:ncol(sgrp)], 2, function(x){
    x <- na.omit(x)
    c("mean"=mean(x), "median"=median(x), "sd"=sd(x), "stde"=stde(x), "RSD"=CV(x))
  })


}

pp4file <- function(grpi, scanpeaks, rtmed, offsets, iwNorm, iwNormFun, ilim,
                    plotP, isotopes, im, singleCheck=TRUE){

  # Sometimes XCMS groups more than 1 peak from the same file in the XCMS
  # grouping stage.
  #we get the peak closet to the median retention time if this happens
  if(nrow(grpi)>1 | singleCheck==TRUE){
    mtchidx <- which(abs(grpi$rt-rtmed)==min(abs(grpi$rt-rtmed)))
    # if two rt the same then pick the one with the highest intensity
    if(length(mtchidx)>1){
      mtchidx <- which(abs(grpi$rt-rtmed)==min(abs(grpi$rt-rtmed)))

    }
    # if still more than one, (i.e. intensitys are the same,
    #just pick the first in list)
    if(length(mtchidx)>1){
      mtchidx <- mtchidx[1]
    }
    target <- grpi[mtchidx, ]
  }else{
    target <- grpi

  }

  # Get the peaks from each scan of the region of interest (ROI)
  roi_scns <- scanpeaks[[target$sample]][target$minscan:target$maxscan]

  mzmax <- target$mz + offsets[1]
  mzmin <- target$mz - offsets[2]

  #get purity for that region
  dfp <- plyr::ldply(roi_scns, pcalc,
                     mzmin=mzmin,
                     mzmax=mzmax,
                     mztarget=target$mz,
                     iwNorm=iwNorm,
                     iwNormFun=iwNormFun,
                     ilim=ilim,
                     targetMinMZ=target$mzmin,
                     targetMaxMZ=target$mzmax,
                     isotopes=isotopes,
                     im=im)

  colnames(dfp) <- c("purity", "pknm")

  scan <- seq(target$minscan, target$maxscan)
  dfp <- cbind(dfp, scan)

  intensity <- plyr::ldply(roi_scns, getEic, target)

  dfp <- cbind(dfp, "intensity"=intensity[,1])

  # Calculate FWHM
  fwhm <- calculateFWHM(dfp)

  if(plotP){
    if (singleCheck){
      tid = target$grpid
    }else{
      tid = target$id
    }
    plotnm <- paste(paste(tid, "sample", target$sample, "mz",
                          round(target$mz, 3), "rt", round(target$rt, 3), sep="_"),
                    "png", sep=".")

    tic <- plyr::ldply(roi_scns, getTic, target = target, minOff = offsets[1],
                       maxOff = offsets[2])
    contamination <- tic$V1-dfp$intensity
    maxi <- max(c(max(contamination), max(dfp$intensity)))

    fpth <- file.path(getwd(),"purityXplots", plotnm)
    png(fpth)


    plot(dfp$scan, dfp$intensity, type = "l", col="red", lwd=2.5,
         xlab="scan number", ylab="intensity", ylim=c(0, maxi))

    lines(dfp$scan, contamination, col="black", lwd=2.5)
    legend("topright", # places a legend at the appropriate place
           lwd=c(2.5,2.5),
           c("EIC of feature", "contamination"),
           col=c("red", "black"))
    abline(v=fwhm[1])
    abline(v=fwhm[2])
    dev.off()
  }



  purityFWHMmedian <- median(dfp$purity[dfp$scan>=fwhm[1] & dfp$scan<=fwhm[2]],
                             na.rm = TRUE)

  purityFWmedian <- median(dfp$purity, na.rm = TRUE)


  pknmFWHMmedian <- median(dfp$pknm[dfp$scan>=fwhm[1] & dfp$scan<=fwhm[2]],
                           na.rm = TRUE)

  pknmFWmedian <- median(dfp$pknm, na.rm = TRUE)

  purityMedianRT <- dfp[dfp$time==target$rtmedscan,]$purity





  return(c("purityFWHMmedian"=purityFWHMmedian, "purityFWmedian"=purityFWmedian,
           "purityMedianRT"=purityMedianRT, "pknmFWHMmedian"=pknmFWHMmedian,
           "pknmFWmedian"=pknmFWmedian))

}


calculateFWHM <- function(df){
  # Calculate FWHM
  xmax <- df$scan[df$i==max(df$intensity)]
  x1 <- df$scan[df$scan < xmax][which.min(abs(df$intensity[df$scan < xmax]-max(df$intensity)/2))]
  x2 <- df$scan[df$scan > xmax][which.min(abs(df$intensity[df$scan > xmax]-max(df$intensity)/2))]
  # http://stackoverflow.com/questions/25015410/r-find-full-width-at-half-maximum-for-a-gausian-density-distribution
  # stack overflow explanation as to performing FWHM for any kind of peak shape (not specifically for gaussian)
  return(c(x1, x2))

}

getEic <- function(roi_scn, target){
  roi_scn <- data.frame(roi_scn)

  sub <- roi_scn[(roi_scn[,1]>=target$mzmin) & (roi_scn[,1]<=target$mzmax),]

  if(nrow(sub)>1){
    # Use the mz value nearest to the
    closeMtch <- sub[which.min(abs(sub[,1] - target$mz)),]
    return(unlist(closeMtch[2]))
  }else if(nrow(sub)==0){
    return(0)
  }else{
    return(sub[,2])
  }
}


getTic <- function(roi_scn, target, minOff, maxOff ){
  roi_scn <- data.frame(roi_scn)

  sub <- roi_scn[(roi_scn[,1]>=target$mz-minOff) & (roi_scn[,1]<=target$mzmax+maxOff),]

  sum(sub[,2])
}



stde <- function(x) sd(x)/sqrt(length(x))
CV <- function(x) (sd(x)/mean(x))*100
