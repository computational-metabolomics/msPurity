#' @title Filter fragmentations spectra associated with an XCMS feature
#'
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
#'
#' @export
setMethod(f="filterFragSpectra", signature="purityA",
          definition = function(pa, ilim=0, plim=0.8, ra=0, snr=3, rmp=FALSE, snmeth='median'){

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

            pa@filter_frag_params <- filter_frag_params

            # reset (incase filterFrag has already been run)
            pa@grped_ms2 <- lapply(pa@grped_ms2, lapply, reset_grped_ms2)

            # Calculate and add flags to matrix
            # Add the purity flag
            pa@grped_df$purity_pass_flag <-  pa@grped_df$inPurity > plim
            pa@grped_ms2 <- plyr::dlply(pa@grped_df, ~grpid, set_purity_spectra, pa)

            # calculate snr, ra and combine all flags
            pa@grped_ms2 <- lapply(pa@grped_ms2, lapply, set_flag_matrix, filter_frag_params=filter_frag_params)

            # add to purityA object

            return(pa)

          }
)

reset_grped_ms2 <- function(m){
  m[,1:2, drop=FALSE]
}

set_purity_spectra <- function(x, pa){
  grpid <- as.character(unique(x$grpid))

  msms_l <- pa@grped_ms2[as.character(grpid)][[1]]

  purity_pass_flag <- x$purity_pass_flag

  x$subsetid <- 1:nrow(x)
  result <- plyr::dlply(x, ~subsetid, set_purity_flag_matrix, msms_l=msms_l)
  result <- unname(result)
  attr(result, "split_type") <- NULL
  attr(result, "split_labels") <- NULL
  return(result)
}



set_purity_flag_matrix <- function(grpinfo, msms_l){
  m <- cbind(msms_l[[grpinfo$subsetid]], grpinfo$purity_pass_flag)
  return(m)
}

set_flag_matrix <- function(x, filter_frag_params){

  snmeth <- filter_frag_params$snmeth
  i_thre <- filter_frag_params$ilim
  ra_thre <- filter_frag_params$ra
  snr_thre <- filter_frag_params$snr
  rmp <- filter_frag_params$rmp

  # Add colnames


  if (snmeth=="median"){
    snr <- x[,2]/median(x[,2])
  }else if(snmeth=="mean"){
    snr <- x[,2]/mean(x[,2])
  }
  ra <- x[,2]/max(x[,2])*100

  intensity_pass_flag <- x[,2]>i_thre
  ra_pass_flag <- ra>ra_thre
  snr_pass_flag <- snr>snr_thre

  x <- cbind(x, snr, ra, intensity_pass_flag, ra_pass_flag, snr_pass_flag)


  colnames(x) <- c('mz', 'i','purity_pass_flag', 'snr', 'ra', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')

  if (nrow(x)==1){
    pass_flag <- sum(x[,c('purity_pass_flag', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')])==4
  }else{
    pass_flag <- rowSums(x[,c('purity_pass_flag', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')])==4
  }

  x <- cbind(x, pass_flag)
  # reoder so it is easier to read
  col_order <- c('mz', 'i', 'snr', 'ra', 'purity_pass_flag', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag', 'pass_flag')

  if (rmp){
    x <- x[x[,'pass_flag']==1,,drop=FALSE]
    if (nrow(x)==0){
      return(NULL)
    }
  }

  return(x[,col_order, drop=FALSE])
}

