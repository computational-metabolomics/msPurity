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
#' @param xset object; xcms object
#' @param cores numeric; number of cores to use
#' @param purityType character; Area and average used for the purity predictions. Options are
#'                               "purityFWHMmedian", "purityFWmedian", "purityFWHMmean", "purityFWmean"
#' @param offsets vector; vector of the isolation window upper and lower offsets
#' @param fileignore vector; vector of files to ignore for the prediction calculation
#' @param xgroups vector; vector of xcms groups to perform prediction on
#' @param iwNorm boolean; if TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5% (0.05)
#' @param plotP boolean; TRUE if plot of the EIC of feature and associated contamination is the be save to the working directory
#' @param mzRback character; backend to use for mzR parsing
#' @param isotopes boolean; TRUE if isotopes are to be removed
#' @param im matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#' @param rtrawColumns boolean; TRUE if the rt_raw values are included as additional columns in the @peaks slot (only required if using the obiwarp)
#' @param singleFile numeric; If just a single file for purity is to be calculated (rather than the grouped XCMS peaks). Uses the index of the files in xcmsSet object. If zero this is ignored.
#' @param saveEIC boolean; If True extracted ion chromatograms will be saved to SQLite database
#' @param sqlitePth character; If saveEIC True, then a path to sqlite database can be used. If NULL then a database will be created in the working directory called eics
#'
#' @return a purityX object containing a dataframe of predicted purity scores
#' @examples
#' msPths <- list.files(system.file("extdata", "lcms", "mzML",
#'                     package="msPurityData"), full.names = TRUE,
#'                     pattern = "LCMS_")
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
                    singleFile=0, rtrawColumns=FALSE, saveEIC=FALSE, sqlitePth=NULL){

  if (singleFile>0){

    ppLCMS <- xcmsSinglePurity(xset, fileidx=singleFile, offsets=offsets, iwNorm=iwNorm,
                               iwNormFun=iwNormFun, ilim=ilim, plotP=plotP, mzRback=mzRback,
                               isotopes=isotopes, im=im, rtrawColumns=rtrawColumns)
  }else{
    ppLCMS <- xcmsGroupPurity(xset, offsets=offsets, fileignore=fileignore, purityType=purityType,
                              cores=cores, xgroups=xgroups, iwNorm=iwNorm, iwNormFun=iwNormFun,
                              ilim=ilim, plotP=plotP, mzRback=mzRback, isotopes=isotopes, im=im,
                              rtrawColumns=rtrawColumns, saveEIC=saveEIC, sqlitePth=sqlitePth)
  }


  return(ppLCMS)

}

xcmsSinglePurity <- function(xset, fileidx, offsets, iwNorm, iwNormFun, ilim, plotP, mzRback, isotopes, im,
                             rtrawColumns){

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
  id <- seq(1, nrow(peaklist))
  peaklist <- cbind(peaklist, minscan, maxscan, id)

  # remove peaks that do not have SN value (means that they will be created from
  #  the 'fillpeaks' function)
  peaklist <- peaklist[!is.na( peaklist[,'sn']),]

  msLevelTracking <- get_mslevel_tracking(filepth)

  for(i in 1:nrow(peaklist)){
    peak <- peaklist[i,]
    lidx <- get_rt_idx(peak, xset, rtrawColumns, msLevelTracking)

    peaklist[i,]['minscan'] <- lidx$rtminidx
    peaklist[i,]['maxscan'] <- lidx$rtmaxidx

  }
  dfp <- data.frame(peaklist)


  if(plotP){
    dir.create(file.path(getwd(), "purityXplots"), showWarnings = FALSE)
  }


  sgrp <- plyr::ddply(dfp, ~ id, pp4file, scanpeaks,
                      rtmed=NA, offsets=offsets, iwNorm=iwNorm, iwNormFun=iwNormFun,
                      ilim=ilim, plotP=plotP, isotopes=isotopes, im=im, singleCheck=FALSE,
                      msLevelTracking=msLevelTracking)

  ppLCMS@predictions <- sgrp


  return(ppLCMS)

}



xcmsGroupPurity <- function(xset, purityType, offsets,
                            fileignore, cores, xgroups,
                            iwNorm, iwNormFun, ilim, plotP, mzRback, isotopes, im, rtrawColumns,
                            saveEIC, sqlitePth){
  # get the filepaths from the xcms object
  filepths <- xset@filepaths

  # get xcms peaklist
  peaklist <- xset@peaks

  peaklist <- cbind(peaklist, 'c_peak_id'=1:nrow(peaklist))

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
  #  the 'fillpeaks' function)
  grouplist <- grouplist[!is.na(grouplist[,'sn']),]

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
  msLevelTracking <- get_mslevel_tracking(filepths)


  for(i in 1:nrow(grouplist)){
    peak <- grouplist[i,]

    lidx <- get_rt_idx(peak, xset, rtrawColumns, msLevelTracking)

    grouplist[i,]['rtraw'] <- lidx$rtraw
    grouplist[i,]['rtminraw'] <- lidx$rtminraw
    grouplist[i,]['rtmaxraw'] <- lidx$rtmaxraw
    grouplist[i,]['rtmedscan'] <- lidx$rtmedidx
    grouplist[i,]['minscan'] <- lidx$rtminidx
    grouplist[i,]['maxscan'] <- lidx$rtmaxidx


  }

  #print(nrow(grouplist))

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
                                   im=im,
                                   xset=xset,
                                   saveEIC=saveEIC,
                                   sqlitePth=sqlitePth,
                                   msLevelTracking=msLevelTracking)

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


get_rt_idx <- function(peak, xset, rtrawColumns, msLevelTracking){

  sid <- as.numeric(peak['sample'])
  raw <- xset@rt$raw[[sid]]
  corrected <- xset@rt$corrected[[sid]]

  if (rtrawColumns){
    rtmedraw <- as.numeric(peak['rt_raw'])
    rtminraw <- as.numeric(peak['rtmin_raw'])
    rtmaxraw <- as.numeric(peak['rtmax_raw'])
    rtmedidx <- msLevelTracking$scan[msLevelTracking$retentionTime==rtmedraw]
    rtminidx <- msLevelTracking$scan[msLevelTracking$retentionTime==rtminraw]
    rtmaxidx <- msLevelTracking$scan[msLevelTracking$retentionTime==rtmaxraw]


  }else{
    rtmed <- as.numeric(peak['rt'])
    rtmin <- as.numeric(peak['rtmin'])
    rtmax <- as.numeric(peak['rtmax'])
    rtmedraw <- raw[which(corrected==rtmed)]
    rtminraw <- raw[which(corrected==rtmin)]
    rtmaxraw <- raw[which(corrected==rtmax)]
    rtmedidx <- msLevelTracking$scan[msLevelTracking$retentionTime==rtmedraw]
    rtminidx <- msLevelTracking$scan[msLevelTracking$retentionTime==rtminraw]
    rtmaxidx <- msLevelTracking$scan[msLevelTracking$retentionTime==rtmaxraw]
  }

  return(list('rtmedidx'=rtmedidx, 'rtminidx'=rtminidx, 'rtmaxidx'=rtmaxidx,
              'rtraw'=rtmedraw, 'rtminraw'=rtminraw,
              'rtmaxraw'=rtmaxraw))

}

predictPurityLCMSloop <- function(grp, average="median", scanpeaks,
                                  offsets, iwNorm, iwNormFun, ilim,
                                  plotP, isotopes, im, xset, saveEIC,
                                  sqlitePth, msLevelTracking){

  # Need to loop through for each file in each group
  grp <- data.frame(grp)
  rtmed <- median(grp$rt)


  sgrp <- plyr::ddply(grp, ~ sample, pp4file, scanpeaks,
                      rtmed, offsets, iwNorm, iwNormFun, ilim, plotP, isotopes, im, xset,
                      singleCheck=TRUE, saveEIC=saveEIC, sqlitePth=sqlitePth, msLevelTracking)

  puritySummary <- apply(sgrp[,2:ncol(sgrp)], 2, function(x){
    x <- na.omit(x)
    c("mean"=mean(x), "median"=median(x), "sd"=sd(x), "stde"=stde(x), "RSD"=CV(x))
  })


}

pp4file <- function(grpi, scanpeaks, rtmed, offsets, iwNorm, iwNormFun, ilim,
                    plotP, isotopes, im, xset, singleCheck=TRUE, saveEIC=FALSE,
                    sqlitePth=NULL, msLevelTracking){

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

  # we wan't to ignore any ms2 or above scans for any of the proceeding calculations
  msLevelTrackingShrt <- msLevelTracking[msLevelTracking$sample==target$sample,]

  dfp <- merge(x = dfp, y = msLevelTrackingShrt, by = "scan")
  dfp <- dfp[dfp$msLevel==1, ]

  mslevel1_indx <- as.numeric(rownames(dfp))

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
    tic <- tic[mslevel1_indx,]

    contamination <- tic-dfp$intensity
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

  if (saveEIC){
    if (!is.null(sqlitePth)){
      con <- DBI::dbConnect(RSQLite::SQLite(), sqlitePth)
    }else{
      con <- DBI::dbConnect(RSQLite::SQLite(),'eics.sqlite')
    }

    colnames(dfp)[colnames(dfp)=='retentionTime'] = 'rt_raw'
    #dfp$rt_raw <- xset@rt$raw[[target$sample]][dfp$scan]

    dfp$rt_corrected <- xset@rt$corrected[[target$sample]][which(xset@rt$raw[[target$sample]] %in% dfp$rt_raw)]

    dfp$grpid <- target$grpid
    dfp$c_peak_id <- target$c_peak_id
    dfp$eicidi <- 1:nrow(dfp)
    DBI::dbWriteTable(con, name='eics', value=dfp, row.names=FALSE, append=TRUE)
    #custom_dbWriteTable(name_pk = 'eicid', fks=NA, table_name = 'eic', df=dfp, con=con)
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
  if(is.null(roi_scn)){
    return(0)
  }else if(nrow(roi_scn)==0){
    return(0)
  }else{
    roi_scn <- data.frame(roi_scn)
  }

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

get_mslevel_tracking <- function(filepths){
  msLevelTracking <- getmrdf_standard_all(filepths)
  msLevelTracking <- msLevelTracking[,c('seqNum','sample', 'msLevel', 'acquisitionNum', 'retentionTime')]
  colnames(msLevelTracking)[1] <- 'scan'
  return(msLevelTracking)
}

getmrdf_standard_all <- function(filepths, backend=NULL){
  mrdf <- plyr::adply(filepths, 1, getmrdf_standard, backend=backend)
  colnames(mrdf)[1] <- 'sample'
  return(mrdf)
}

# get the standard mrdf
getmrdf_standard <- function(filepth, backend=NULL){
  if(is.null(backend)){
    mr <- mzR::openMSfile(filepth)
  }else{
    mr <- mzR::openMSfile(filepth, backend=backend)
  }
  return(mzR::header(mr))

}
