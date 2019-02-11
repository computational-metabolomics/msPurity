#' @title Average and filter fragmentation spectra for each XCMS feature within a MS data file
#'
#' @description
#'
#' Average and filter fragmentation spectra for each XCMS feature within a MS data file.
#'
#' The default approach is to use hierarchical clustering where peaks within a set ppm tolerance will be clustered.
#'
#' The clustered peaks are then averaged (or summed) and filtered.
#'
#'
#' @aliases averageIntraFragSpectra
#'
#'
#'
#' @param pa object; purityA object
#' @param cores numeric; Number of cores for multiprocessing
#' @param plim numeric; min purity of precursor for fragmentation spectra scan to be included
#' @param ppm numeric; ppm threshold to average within each file
#' @param minnum numeric; minimum number of times peak is present across fragmentation spectra within each file
#' @param minfrac numeric; minimum ratio of the peak fraction (peak count / total peaks) within each file
#' @param ra numeric; minimum relative abundance of the peak within each file
#' @param snr numeric; minimum signal-to-noise of the peak within each file
#' @param snr_pre numeric;  minimum signal-to-noise prior to averaging
#' @param ra_pre numeric;  minimum relative abundance prior to averaging
#'
#' @param av character; type of averaging to use (median or mean)
#' @param sum_i boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param remove_peaks boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#'
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths, interpol = "linear")
#' pa <- frag4feature(pa, xset)
#' pa <- averageIntraFragSpectra(pa)
#'
#' @export
setMethod(f="averageIntraFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                                snr_pre=0, ra_pre=0, av='median', sum_i=TRUE,  plim=0.5, remove_peaks=FALSE, cores=1
                                ){

            pa@av_intra_params$minfrac = minfrac
            pa@av_intra_params$minnum = minnum
            pa@av_intra_params$ppm = ppm
            pa@av_intra_params$snr = snr
            pa@av_intra_params$ra = ra

            pa@av_intra_params$av = av
            pa@av_intra_params$sum_i = sum_i
            pa@av_intra_params$plim = plim

            pa@av_intra_params$ra_pre = ra_pre
            pa@av_intra_params$snr_pre = snr_pre

            pa@av_intra_params$cores = cores
            pa@av_intra_params$remove_peaks = remove_peaks

            return(average_xcms_grouped_msms(pa, "intra"))

          }
)


#' @title Average and filter fragmentation spectra for each XCMS feature across MS data files
#'
#' @description
#'
#' Average and filter fragmentation spectra for each XCMS feature accross MS data files.
#'
#' The default approach is to use hierarchical clustering where peaks within a set ppm tolerance will be clustered.
#'
#' The clustered peaks are then averaged (or summed) and filtered.
#'
#'
#' @aliases averageInterFragSpectra
#'
#'
#'
#' @param pa object; purityA object
#' @param cores numeric; Number of cores for multiprocessing
#' @param plim numeric; min purity of precursor for fragmentation spectra scan to be included
#' @param ppm numeric; ppm threshold to average across files
#' @param minnum numeric; minimum number of times peak is present across fragmentation spectra across files
#' @param minfrac numeric; minimum ratio of the peak fraction (peak count / total peaks) across files
#' @param ra numeric; minimum relative abundance of the peak across files
#' @param snr numeric; minimum signal-to-noise of the peak across files
#'
#' @param av character; type of averaging to use (median or mean)
#' @param sum_i boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param remove_peaks boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#'
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths, interpol = "linear")
#' pa <- frag4feature(pa, xset)
#' pa <- averageIntraFragSpectra(pa)
#' pa <- averageInterFragSpectra(pa)
#'
#' @export
setMethod(f="averageInterFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                                av='median', sum_i=TRUE,  plim=0.5, remove_peaks=FALSE, cores=1
          ){

            pa@av_inter_params$minfrac = minfrac
            pa@av_inter_params$minnum = minnum
            pa@av_inter_params$ppm = ppm
            pa@av_inter_params$snr = snr
            pa@av_inter_params$ra = ra

            pa@av_inter_params$av = av
            pa@av_inter_params$sum_i = sum_i
            pa@av_inter_params$plim = plim

            pa@av_inter_params$cores = cores
            pa@av_inter_params$remove_peaks = remove_peaks

            if (is.null(pa@av_spectra[[names(pa@grped_ms2)[1]]][["av_intra"]])){
              stop("Apply averageIntraFragSpectra first")
            }

            return(average_xcms_grouped_msms(pa, "inter"))

          }
)


#' @title Average and filter fragmentation spectra for each XCMS feature within and accross MS data files (ignoring intra and inter relationships)
#'
#' @description
#'
#' Average and filter fragmentation spectra for each XCMS feature within and accross MS data files (ignoring intra and inter relationships).
#'
#' The default approach is to use hierarchical clustering where peaks within a set ppm tolerance will be clustered.
#'
#' The clustered peaks are then averaged (or summed) and filtered.
#'
#'
#' @aliases averageAllFragSpectra
#'
#'
#'
#' @param pa object; purityA object
#' @param cores numeric; Number of cores for multiprocessing
#' @param plim numeric; min purity of precursor for fragmentation spectra scan to be included
#' @param ppm numeric; ppm threshold to average across all scans (ignoring intra and inter relationships)
#' @param minnum numeric; minimum number of times peak is present across all fragmentation spectra (ignoring intra and inter relationships)
#' @param minfrac numeric;minimum ratio of the peak fraction (peak count / total peaks) across all (ignoring intra and inter relationships)
#' @param ra numeric; minimum relative abundance of the peak fraction across all (ignoring intra and inter relationships)
#' @param snr numeric;  minimum signal-to-noise of the peak across all (ignoring intra and inter relationships)
#' @param snr_pre numeric;  minimum signal-to-noise prior to averaging
#' @param ra_pre numeric;  minimum relative abundance prior to averaging
#'
#' @param av character; type of averaging to use (median or mean)
#' @param sum_i boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param remove_peaks boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#'
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths, interpol = "linear")
#' pa <- frag4feature(pa, xset)
#' pa <- averageAllFragSpectra(pa)
#'
#' @export
setMethod(f="averageAllFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                                snr_pre=0, ra_pre=0, av='median', sum_i=TRUE,  plim=0.5, remove_peaks=FALSE, cores=1
          ){

            pa@av_all_params$minfrac = minfrac
            pa@av_all_params$minnum = minnum
            pa@av_all_params$ppm = ppm
            pa@av_all_params$snr = snr
            pa@av_all_params$ra = ra

            pa@av_all_params$av = av
            pa@av_all_params$sum_i = sum_i
            pa@av_all_params$plim = plim

            pa@av_all_params$ra_pre = ra_pre
            pa@av_all_params$snr_pre = snr_pre

            pa@av_all_params$cores = cores
            pa@av_all_params$remove_peaks = remove_peaks

            return(average_xcms_grouped_msms(pa, "all"))

          }
)



average_xcms_grouped_msms <- function(pa, av_level){

  if(pa@cores>1){
    cl <- parallel::makeCluster(pa@cores)
    doSNOW::registerDoSNOW(cl)
    para = TRUE
  }else{
    para = FALSE
  }

  av_spectra <- plyr::alply(names(pa@grped_ms2), 1, average_xcms_grouped_msms_indiv, pa=pa, av_level=av_level, .parallel = para)
  names(av_spectra) <- names(pa@grped_ms2)

  pa@av_spectra <- av_spectra

  return(pa)

}



average_xcms_grouped_msms_indiv <- function(grp_idx, pa, av_level){

  ##############################################################################
  # Get the appropiate details for the xcms grouped feature from purityA object
  ##############################################################################
  grped_info <- pa@grped_df[pa@grped_df==as.numeric(grp_idx),]
  grped_spectra <- pa@grped_ms2[as.character(grp_idx)][[1]]

  grped_info$index <- 1:nrow(grped_info)
  names(grped_spectra) <- 1:length(grped_spectra)

  ##############################################################################
  # Create a new dataframe with only the valid info and frag spectra
  ##############################################################################
  grped_spectra <- plyr::llply(grped_spectra, data.frame)
  df <- data.frame(do.call("rbind", grped_spectra))

  colnames(df) <- c('mz', 'i')
  df$index <-   rep(seq_along(grped_spectra), sapply(grped_spectra, nrow))

  spectra_to_average <- merge(df, grped_info[, c('grpid', 'sample', 'cid', 'index', 'inPurity')], by = "index")

  # Set return variable to empty list or already existing results
  if (!is.null(pa@av_spectra[[as.character(grp_idx)]][["av_intra"]])){
    av_intra = pa@av_spectra[[as.character(grp_idx)]][["av_intra"]]
  } else {
    av_intra = NULL
  }

  # Set return variable to empty list or already existing results
  if (!is.null(pa@av_spectra[[as.character(grp_idx)]][["av_inter"]])){
    av_inter = pa@av_spectra[[as.character(grp_idx)]][["av_inter"]]
  } else {
    av_inter = NULL
  }

  # Set return variable to empty list or already existing results
  if (!is.null(pa@av_spectra[[as.character(grp_idx)]][["av_all"]])){
    av_all = pa@av_spectra[[as.character(grp_idx)]][["av_all"]]
  } else {
    av_all = NULL
  }

  # filter out peaks below precursor ion purity thres
  if (av_level=="intra"){
    plim = pa@av_intra_params$plim
  } else if (av_level=="inter"){
    plim = pa@av_inter_params$plim
  } else if (av_level=="all"){
    plim = pa@av_all_params$plim
  } else {
    stop("Incorrect av_level for averaging fragmentation spectra; use intra, inter or all")
  }


  spectra_to_average <- spectra_to_average[spectra_to_average$inPurity>plim, ]

  ##############################################################################
  # Performing averaging
  ##############################################################################
  # clustering requires data to be in order of mz
  spectra_to_average <- spectra_to_average[order(spectra_to_average$mz),]

  if (av_level=="intra"){
    # Average by sample (file)

    av_intra <- plyr::dlply(spectra_to_average, ~sample, average_spectra,
                                            cores=1,
                                            ppm=pa@av_intra_params$ppm,
                                            minnum=pa@av_intra_params$minnum,
                                            sum_i=pa@av_intra_params$sum_i,
                                            minfrac=pa@av_intra_params$minfrac,
                                            snthr=pa@av_intra_params$snr,
                                            rathr=pa@av_intra_params$ra,
                                            rathr_pre= pa@av_intra_params$ra_pre,
                                            snrthr_pre= pa@av_intra_params$snr_pre)
    if (pa@av_intra_params$remove_peaks){
      av_intra  <- plyr::llply(av_intra , function(x){x[x$pass_flag,]})
    }

  } else if (av_level=="inter") {

    av_intra_df <- plyr::ldply(av_intra, .id = 'sample', function(x){x[x$pass_flag,]})

    # Average the averaged spectra across files
    av_inter <- average_spectra(av_intra_df,
                                 indx='sample',
                                 cores=1,
                                 ppm=pa@av_inter_params$ppm,
                                 minnum=pa@av_inter_params$minnum,
                                 sum_i=pa@av_inter_params$sum_i,
                                 minfrac=pa@av_inter_params$minfrac,
                                 snthr=pa@av_inter_params$snr,
                                 rathr=pa@av_inter_params$ra
                               )

    if (pa@av_inter_params$remove_peaks){
      av_inter <- av_inter[av_inter$pass_flag,]
    }

  } else if (av_level=="all") {

    # add additional column used later if filtering applied
    # Average the averaged spectra across everything (ignore intra and inter )
    av_all <- average_spectra(spectra_to_average,
                                cores=1,
                                ppm=pa@av_all_params$ppm,
                                minnum=pa@av_all_params$minnum,
                                sum_i=pa@av_all_params$sum_i,
                                minfrac=pa@av_all_params$minfrac,
                                snthr=pa@av_all_params$snr,
                                rathr=pa@av_all_params$ra,
                                rathr_pre= pa@av_all_params$ra_pre,
                                snrthr_pre= pa@av_all_params$snr_pre)
    if (pa@av_all_params$remove_peaks){
      av_all <- av_all[av_all$pass_flag,]
    }
  } else {

    stop("Incorrect av_level for averaging fragmentation spectra; use intra, inter or all")

  }

  return(list('av_intra'=av_intra , 'av_inter'=av_inter, 'av_all'=av_all))
}



average_spectra <- function(spectra, indx='index', ppm, cores, minnum, sum_i,
                            minfrac, snthr, snmeth='median', rathr, rathr_pre=NULL, snrthr_pre=NULL){
  if (nrow(spectra)==0){
    return(NULL)
  }

  if (snmeth=="median"){
    spectra$snr <- spectra$i/median(spectra$i)
  }else if(snmeth=="mean"){
    spectra$snr <- spectra$i/mean(spectra$i)
  }

  spectra$ra <- spectra$i/max(spectra$i)*100


  if (!is.null(rathr_pre)){
    spectra <- spectra[spectra$ra>rathr_pre, ]
  }

  if (!is.null(snrthr_pre)){
    spectra <- spectra[spectra$snr>snrthr_pre, ]
  }

  if (nrow(spectra)==0){
    return(NULL)
  }


  mz <- spectra$mz


  # Cluster the peaks togther
  spectra$cl <- clustering(mz, clustType = 'hc', cores = cores, ppm = ppm)

  averages <- plyr::ddply(spectra, ~ cl,
                          averageCluster, av="median", minnum=minnum,
                          missingV="ignore", totalScans=length(unique(spectra[,indx])), normTIC=FALSE,
                          sumI=sum_i)


  averages$frac <- averages$count/averages$total

  averages$snr_pass_flag <- averages$snr > snthr

  averages$minfrac_pass_flag <- averages$frac > minfrac

  averages$ra_pass_flag <- averages$ra > rathr

  averages$pass_flag <- (averages$minfrac_pass_flag & averages$snr_pass_flag & averages$ra_pass_flag)

  return(averages)
}


