#' @title Average fragemntation spectra across XCMS features
#'
#' @description
#'
#' Average fragemntation spectra across XCMS features. Averaging can be done within file (intra), across files (inter)
#' and independantly of the files (all).
#'
#' The default approach is to use hierarchical clustering where peaks within a set ppm tolerance will be clustered.
#'
#' The clustered peaks are then averaged (or summed) and filtered.
#'
#'
#' @aliases averageFragmentation
#'
#'
#'
#' @param pa object; purityA object
#' @param cores numeric; Number of cores for multiprocessing
#' @param plim numeric; min purity of precursor for fragmentation spectra scan to be included
#' @param ppm_intra numeric; ppm threshold to average within each file
#' @param ppm_inter numeric; ppm threshold to average across files
#' @param ppm_all numeric; ppm threshold to average across all scans (ignoring intra and inter relationships)
#' @param minnum_intra numeric; minimum number of times peak is present across fragmentation spectra within each file
#' @param minnum_inter numeric; minimum number of times peak is present across fragmentation spectra across files
#' @param minnum_all numeric; minimum number of times peak is present across all fragmentation spectra (ignoring intra and inter relationships)
#' @param minfrac_intra numeric; minimum ratio of the peak fraction (peak count / total peaks) within each file
#' @param minfrac_inter numeric; minimum ratio of the peak fraction (peak count / total peaks) across files
#' @param minfrac_all numeric;minimum ratio of the peak fraction (peak count / total peaks) across all (ignoring intra and inter relationships)
#' @param av character; type of averaging to use (median or mean)
#' @param sum_i boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param remove_peaks boolean; TRUE if peaks are to be removed that do not meet the minfrac criteria. Otherwise they will just be flagged
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
#' pa <- averageFragmentation(pa)
#'
#' @export
setMethod(f="averageFragmentation", signature="purityA",
          definition = function(pa, minfrac_intra=0.5, minfrac_inter=0.5,  minnum_intra=1,
                                minnum_inter=1, ppm_intra=5, ppm_inter=5, minfrac_all=0.5, minnum_all=1,
                                ppm_all=5, av='median', sum_i=TRUE,  plim=0.5, remove_peaks=FALSE, cores=1){

            pa@av_params$minfrac_intra = minfrac_intra
            pa@av_params$minfrac_inter = minfrac_inter
            pa@av_params$minnum_intra = minnum_intra
            pa@av_params$minnum_inter = minnum_inter
            pa@av_params$ppm_intra = ppm_intra
            pa@av_params$ppm_inter = ppm_inter

            pa@av_params$minfrac_all = minfrac_all
            pa@av_params$minnum_all = minnum_all
            pa@av_params$ppm_all = ppm_all

            pa@av_params$av = av
            pa@av_params$sum_i = sum_i
            pa@av_params$plim = plim

            pa@av_params$cores = cores
            pa@av_params$remove_peaks = remove_peaks

            return(average_xcms_grouped_msms_all(pa))

          }
)

average_xcms_grouped_msms_all <- function(pa){

  if(pa@cores>1){
    cl <- parallel::makeCluster(pa@cores)
    doSNOW::registerDoSNOW(cl)
    para = TRUE
  }else{
    para = FALSE
  }

  av_spectra <- plyr::alply(names(pa@grped_ms2), 1, average_xcms_grouped_msms_indiv, pa=pa, .parallel = para)

  names(av_spectra) <- names(pa@grped_ms2)

  pa@av_spectra <- av_spectra

  return(pa)

}



average_xcms_grouped_msms_indiv <- function(grp_idx, pa){
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

  # filter out peaks below precursor ion purity thres
  spectra_to_average <- spectra_to_average[spectra_to_average$inPurity>pa@av_params$plim, ]

  ##############################################################################
  # Performing averaging
  ##############################################################################
  # clustering requires data to be in order of mz
  spectra_to_average <- spectra_to_average[order(spectra_to_average$mz),]

  # Average by sample (file)
  av_intra <- plyr::dlply(spectra_to_average, ~sample, average_spectra,
                                          cores=1,
                                          ppm=pa@av_params$ppm_intra,
                                          minnum=pa@av_params$minnum_intra,
                                          sum_i=pa@av_params$sum_i,
                                          minfrac=pa@av_params$minfrac_intra)


  # Average the averaged spectra across files
  av_inter <- average_spectra( plyr::ldply(av_intra, function(x){x[x$minfrac_pass_flag,]}),
                            indx='sample',
                            cores=1,
                            ppm=pa@av_params$ppm_inter,
                            minnum=pa@av_params$minnum_inter,
                            sum_i=pa@av_params$sum_i,
                            minfrac=pa@av_params$minfrac_inter)

  # add additional column used later if filtering applied

  # Average the averaged spectra across everything (ignore intra and inter )
  av_all <- average_spectra(spectra_to_average,
                              cores=1,
                              ppm=pa@av_params$ppm_all,
                              minnum=pa@av_params$minnum_all,
                              sum_i=pa@av_params$sum_i,
                              minfrac=pa@av_params$minfrac_all)

  if (pa@av_params$remove_peaks){
    av_intra  <- plyr::llply(av_intra , function(x){x[x$minfrac_pass_flag,]})
    av_all <- av_all[av_all$minfrac_pass_flag,]
    av_inter <- av_inter[av_inter$minfrac_pass_flag,]
  }


  return(list('av_intra'=av_intra , 'av_inter'=av_inter, 'av_all'=av_all))
}



average_spectra <- function(spectra, indx='index', ppm, cores, minnum, sum_i, minfrac){
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
  averages$minfrac_pass_flag <- averages$frac > minfrac

  return(averages)
}


