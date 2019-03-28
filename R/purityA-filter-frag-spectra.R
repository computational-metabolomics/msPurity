#' @title Filter fragmentations spectra associated with an XCMS feature
#' @aliases filterFragSpectra
#' @description
#'
#' Flag and filter features based on signal-to-noise ratio, relative abundance, intensity threshold and
#' precursor ion purity of precursor.
#'
#' The grp_spectra slot to add the columns snr, ra, snr_flag, ra_flag, purity_flag, intensity_flag pass_flag
#'
#' This filtering occurs at the scan level (i.e. should be run prior to any averaging)
#'
#' @param pa object; purityA object
#' @param ilim numeric; min intensity of a peak
#' @param plim numeric; min precursor ion purity of the associated precursor for fragmentation spectra scan
#' @param ra numeric; minimum relative abundance of a peak
#' @param snr numeric; minimum signal-to-noise of a peak  peak within each file
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged.
#' @param snmeth character; Method to calculate signal to noise ration (either median or mean)
#' @param allfrag boolean; Whether to filter on all fragmentation spectra or or just the fragmentation spectra grouped to XCMS feature
#'
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- frag4feature(pa, xset)
#' pa <- filterFragSpectra(pa)
#'
#' @return purityA object with additional columns for the fragmentation spectra in the grped_msms slot
#'
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

