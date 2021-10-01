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


#' @title Filter fragmentation spectra associated with an XCMS feature
#' @aliases filterFragSpectra
#' @description
#' **General**
#'
#' Flag and filter features based on signal-to-noise ratio, relative abundance, intensity threshold and purity of the precursor ion.
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> (optionally) xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'    +  --- *Older versions of XCMS* --- (mzML files) -> xcms.xcmsSet -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> frag4feature -> **filterFragSpectra** -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' @param pa object; purityA object
#' @param ilim numeric; min intensity of a peak
#' @param plim numeric; min precursor ion purity of the associated precursor for fragmentation spectra scan
#' @param ra numeric; minimum relative abundance of a peak
#' @param snr numeric; minimum signal-to-noise of a peak within each file
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged.
#' @param snmeth character; Method to calculate signal to noise ration (either median or mean)
#' @param allfrag boolean; Whether to filter on all fragmentation spectra or just the fragmentation spectra grouped to XCMS feature
#'
#' @examples
#' #====== XCMS =================================
#' ## Read in MS data
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#' #           package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #ms_data = readMSData(msmsPths, mode = 'onDisk', msLevel. = 1)
#'
#' ## Find peaks in each file
#' #cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3, 30))
#' #xcmsObj  <- xcms::findChromPeaks(ms_data, param = cwp)
#'
#' ## Optionally adjust retention time
#' #xcmsObj  <- adjustRtime(xcmsObj , param = ObiwarpParam(binSize = 0.6))
#'
#' ## Group features across samples
#' #pdp <- PeakDensityParam(sampleGroups = c(1, 1), minFraction = 0, bw = 30)
#' #xcmsObj <- groupChromPeaks(xcmsObj , param = pdp)
#'
#' #====== msPurity ============================
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xcmsObj)
#' #pa <- filterFragSpectra(pa)
#'
#' ## Run from previously generated data
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "2_frag4feature_pa.rds", package="msPurity"))
#' pa <- filterFragSpectra(pa)
#'
#' @return Returns a purityA object with the pa@@grped_msms spectra matrices are updated with the following columns
#'
#' * snr: Signal to noise ratio (calculated at scan level)
#' * ra: Relative abundance (calculated at scan level)
#' * purity_pass_flag: Precursor ion purity flag (1 pass, 0 fail)
#' * intensity_pass_flag: Intensity flag (1 pass, 0 fail)
#' * snr_pass_flag: Signal-to-noise pass flag (1 pass, 0 fail)
#' * ra_pass_flag: Relative abundance pass flag (1 pass, 0 fail)
#' * pass_flag: Overall pass flag, all flags must pass for this to pass (1 pass, 0 fail)
#' @md
#' @export
setMethod(f="filterFragSpectra", signature="purityA",
          definition = function(pa, ilim=0, plim=0.8, ra=0, snr=3, rmp=FALSE, snmeth='median', allfrag=FALSE){

            if ((!length(pa@filter_frag_params)==0) && (pa@filter_frag_params$rmp)){
              message('Fragmentation peaks have been previously filtered and removed - function can\'t be performed')
              return(pa)
            }

            filter_frag_params = list()
            filter_frag_params$ilim = ilim
            filter_frag_params$plim = plim
            filter_frag_params$ra = ra
            filter_frag_params$snr = snr
            filter_frag_params$snmeth = snmeth
            filter_frag_params$rmp = rmp
            filter_frag_params$allfrag = allfrag

            pa@filter_frag_params <- filter_frag_params

            # reset (incase filterFrag has already been run)
            pa@grped_ms2 <- lapply(pa@grped_ms2, lapply, resetGrpedMS2)

            # Calculate and add flags to matrix
            # Add the purity flag
            pa@grped_df$purity_pass_flag <-  pa@grped_df$inPurity > plim
            pa@grped_ms2 <- plyr::dlply(pa@grped_df, ~grpid, setPuritySpectraGrp, pa)

            # calculate snr, ra and combine all flags
            pa@grped_ms2 <- lapply(pa@grped_ms2, lapply, setFlagMatrix, filter_frag_params=filter_frag_params)

            # if allfrag is TRUE calculate the ms2 for all fragmentation spectra (note this is duplicating the
            # processing at the momement - for sake of clarity of the code.
            # So Ideally should be updated but waiting on the backend refactor to update though)
            if (allfrag){
              # add to purityA object
              scanpeaksFrag <- getScanPeaks(pa)

              # get purity flag
              pa@puritydf$purity_pass_flag <- pa@puritydf$inPurity > plim
              scanpeaksFrag <- merge(pa@puritydf[,c('pid', 'inPurity', 'purity_pass_flag')], scanpeaksFrag, by='pid')
              scanpeaksFrag <- scanpeaksFrag[,c("pid","sid","fileid", "scan", "mz", "i", "type",  "purity_pass_flag"), drop=FALSE]


              pa@all_frag_scans <- plyr::ddply(scanpeaksFrag, ~pid, setFlagMatrix, filter_frag_params=filter_frag_params)

            }

            return(pa)

          }
)




getScanPeaks <- function(pa){
  scaninfo <- pa@puritydf
  fileList <- pa@fileList
  filedf <- data.frame(filepth=fileList,
                       fileid=seq(1, length(fileList)))
  scanpeaksFrag <- plyr::ddply(filedf, ~ fileid, scanPeaksFromfiledf)
  comb <- paste(scanpeaksFrag[,1], scanpeaksFrag[,2], sep=' ')
  scanpeaksFrag <- cbind(1:nrow(scanpeaksFrag), cumsum(!duplicated(comb)), scanpeaksFrag)
  colnames(scanpeaksFrag) <- c('sid','pid', 'fileid', 'scan', 'mz', 'i')
  scanpeaksFrag$type <- 'scan'
  return(scanpeaksFrag)
}

scanPeaksFromfiledf <- function(x){
  mr <- mzR::openMSfile(as.character(x$filepth))
  scanpeaks <- mzR::peaks(mr)
  scans <- mzR::header(mr)
  names(scanpeaks) <- seq(1, length(scanpeaks))
  scanpeaks_df <- plyr::ldply(scanpeaks[scans$seqNum[scans$msLevel>1]], .id=TRUE)
}




resetGrpedMS2 <- function(m){
  m[,1:2, drop=FALSE]
}

setPuritySpectraGrp <- function(x, pa){
  grpid <- as.character(unique(x$grpid))

  msms_l <- pa@grped_ms2[as.character(grpid)][[1]]

  purity_pass_flag <- x$purity_pass_flag

  x$subsetid <- 1:nrow(x)
  result <- plyr::dlply(x, ~subsetid, setPurityFlagMatrix, msms_l=msms_l)
  result <- unname(result)
  attr(result, "split_type") <- NULL
  attr(result, "split_labels") <- NULL
  return(result)
}



setPurityFlagMatrix <- function(grpinfo, msms_l){
  m <- cbind(msms_l[[grpinfo$subsetid]], grpinfo$purity_pass_flag)
  return(m)
}

setFlagMatrix <- function(x, filter_frag_params){

  snmeth <- filter_frag_params$snmeth
  i_thre <- filter_frag_params$ilim
  ra_thre <- filter_frag_params$ra
  snr_thre <- filter_frag_params$snr
  rmp <- filter_frag_params$rmp

  # check if columns defined for intensity
  if('i' %in% colnames(x)){
    intensityIdx <- which(colnames(x)=='i')
  }else{
    intensityIdx <- 2
  }

  if (snmeth=="median"){
    snr <- x[,intensityIdx]/median(x[,intensityIdx])
  }else if(snmeth=="mean"){
    snr <- x[,intensityIdx]/mean(x[,intensityIdx])
  }

  ra <- x[,intensityIdx]/max(x[,intensityIdx])*100

  intensity_pass_flag <- x[,intensityIdx]>i_thre
  ra_pass_flag <- ra>ra_thre
  snr_pass_flag <- snr>snr_thre

  x <- cbind(x, snr, ra, intensity_pass_flag, ra_pass_flag, snr_pass_flag)


  # check if columns available for sid and pid
  if('pid' %in% colnames(x)){
    col_order <-   c('sid', 'pid', 'fileid', 'scan', 'mz', 'i', 'snr', 'ra','type', 'purity_pass_flag', 'intensity_pass_flag','ra_pass_flag', 'snr_pass_flag', 'pass_flag')

  }else{
    colnames(x) <- c('mz', 'i','purity_pass_flag', 'snr', 'ra', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')
    col_order <-   c('mz', 'i', 'snr', 'ra', 'purity_pass_flag', 'intensity_pass_flag','ra_pass_flag', 'snr_pass_flag', 'pass_flag')
  }




  if (nrow(x)==1){
    pass_flag <- sum(unlist(x[,c('purity_pass_flag', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')]))==4
  }else{
    pass_flag <- rowSums(x[,c('purity_pass_flag', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')])==4
  }

  x <- cbind(x, pass_flag)
  # reoder so it is easier to read


  if (rmp){
    x <- x[x[,'pass_flag']==1,,drop=FALSE]
    if (nrow(x)==0){
      return(NULL)
    }
  }

  return(x[,col_order, drop=FALSE])
}

