#' @title Using a purityA object, average and filter fragmentation spectra for each XCMS feature within a MS data file
#' @aliases averageIntraFragSpectra
#' @description
#'
#' **General**
#'
#' Average and filter fragmentation spectra for each XCMS feature within a MS data file.
#'
#' The averaging is performed using hierarchical clustering of the m/z values of each peaks, where m/z values within a set ppm tolerance will be clustered. The clustered peaks are then averaged (or summed).
#'
#' The fragmentation can be filtered on the averaged spectra (with the arguments snr, rsd, minfrac and ra)
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
#'  * Fragmentation processing
#'    + (xset, pa) -> frag4feature -> filterFragSpectra -> **averageIntraFragSpectra** -> averageIntraFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#'
#' @param pa object; purityA object
#' @param ppm numeric; ppm threshold to average within each file
#' @param minnum numeric; minimum number of times peak is present across fragmentation spectra within each file
#' @param minfrac numeric; minimum ratio of the peak fraction (peak count / total peaks) within each file
#' @param ra numeric; minimum relative abundance of the peak within each file
#' @param snr numeric; minimum signal-to-noise of the peak within each file
#' @param av character; type of averaging to use (median or mean)
#' @param sumi boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#' @param cores numeric; Number of cores for multiprocessing
#'
#' @return Returns a purityA object (pa) with the following slots now with data
#'
#' * pa@@av_spectra: the average spectra is recorded here stored as a list. e.g. "pa@av_spectra$`1`$av_intra$`1`" would give the average spectra for grouped feature 1 and for file 1.
#' * pa@@av_intra_params: The parameters used are recorded here
#'
#' Each spectra in the av_spectra list contains the following columns:
#'
#' * cl: id of clustered (averaged) peak
#' * mz: average m/z
#' * i: average intensity
#' * snr: average signal to noise ratio
#' * rsd: relative standard deviation
#' * count: number of clustered peaks
#' * total: total number of potential scans to be used for averaging
#' * inPurity: average precursor ion purity
#' * ra: average relative abundance
#' * frac: the fraction of clustered peaks (e.g. the count/total)
#' * snr_pass_flag: TRUE if snr threshold criteria met
#' * minfrac_pass_flag: TRUE if minfrac threshold criteria
#' * ra_pass_flag: TRUE if ra threshold criteria met
#' * pass_flag: TRUE if all threshold criteria met
#'
#' @examples
#'
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' #xset <- xcms::group(xset)
#' #xset <- xcms::retcor(xset)
#' #xset <- xcms::group(xset)
#'
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xset)
#' pa <- readRDS(system.file("extdata", "tests", "purityA", "2_frag4feature_pa.rds", package="msPurity"))
#' pa <- averageIntraFragSpectra(pa)
#' @md
#' @export
setMethod(f="averageIntraFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                               av='median', sumi=TRUE, rmp=FALSE, cores=1
                                ){

            pa@av_intra_params$minfrac = minfrac
            pa@av_intra_params$minnum = minnum
            pa@av_intra_params$ppm = ppm
            pa@av_intra_params$snr = snr
            pa@av_intra_params$ra = ra

            pa@av_intra_params$av_type = av
            pa@av_intra_params$sumi = sumi


            pa@av_intra_params$cores = cores
            pa@av_intra_params$rmp = rmp

            return(average_xcms_grouped_msms(pa, "intra"))

          }
)


#' @title Using a purityA object, average and filter fragmentation spectra for each XCMS feature across multiple MS data files
#' @aliases averageInterFragSpectra
#' @description
#'
#' **General**
#'
#' Average and filter fragmentation spectra for each XCMS feature across MS data files. This can only be run after averageIntraFragSpectra has been used.
#'
#' The averaging is performed using hierarchical clustering of the m/z values of each peaks, where m/z values within a set ppm tolerance will be clustered. The clustered peaks are then averaged (or summed).
#'
#' The fragmentation can be filtered on the averaged spectra (with the arguments snr, rsd, minfrac and ra)
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
#'  * Fragmentation processing
#'    + (xset, pa) -> frag4feature -> filterFragSpectra -> averageIntraFragSpectra -> **averageInterFragSpectra** -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#'
#'
#' @param pa object; purityA object
#' @param cores numeric; Number of cores for multiprocessing
#' @param ppm numeric; ppm threshold to average across files
#' @param minnum numeric; minimum number of times peak is present across fragmentation spectra across files
#' @param minfrac numeric; minimum ratio of the peak fraction (peak count / total peaks) across files
#' @param ra numeric; minimum relative abundance of the peak across files
#' @param snr numeric; minimum signal-to-noise of the peak across files
#'
#' @param av character; type of averaging to use (median or mean)
#' @param sumi boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#'
#' @return Returns a purityA object (pa) with the following slots now with data
#'
#' * pa@@av_spectra: the average spectra is recorded here stored as a list. e.g. "pa@@av_spectra$`1`$av_inter" would give the average spectra for grouped feature 1
#' * pa@@av_intra_params: The parameters used are recorded here
#'
#' Each spectra in the av_spectra list contains the following columns:
#' *
#' * cl: id of clustered (averaged) peak
#' * mz: average m/z
#' * i: average intensity
#' * snr: average signal to noise ratio
#' * rsd: relative standard deviation
#' * count: number of clustered peaks
#' * total: total number of potential scans to be used for averaging
#' * inPurity: average precursor ion purity
#' * ra: average relative abundance
#' * frac: the fraction of clustered peaks (e.g. the count/total)
#' * snr_pass_flag: TRUE if snr threshold criteria met
#' * minfrac_pass_flag: TRUE if minfrac threshold criteria
#' * ra_pass_flag: TRUE if ra threshold criteria met
#' * pass_flag: TRUE if all threshold criteria met
#'
#' @examples
#'
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' #xset <- xcms::group(xset)
#' #xset <- xcms::retcor(xset)
#' #xset <- xcms::group(xset)
#'
#' #pa  <- purityA(msmsPths, interpol = "linear")
#' #pa <- frag4feature(pa, xset)
#' #pa <- averageIntraFragSpectra(pa)
#' pa <- readRDS(system.file("extdata", "tests", "purityA", "4_averageIntraFragSpectra_no_filter_pa.rds", package="msPurity"))
#' pa <- averageInterFragSpectra(pa)
#' @md
#' @export
setMethod(f="averageInterFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                                av='median', sumi=TRUE,  rmp=FALSE, cores=1
          ){

            pa@av_inter_params$minfrac = minfrac
            pa@av_inter_params$minnum = minnum
            pa@av_inter_params$ppm = ppm
            pa@av_inter_params$snr = snr
            pa@av_inter_params$ra = ra

            pa@av_inter_params$av_type = av
            pa@av_inter_params$sumi = sumi


            pa@av_inter_params$cores = cores
            pa@av_inter_params$rmp = rmp

            if (is.null(pa@av_spectra[[names(pa@grped_ms2)[1]]][["av_intra"]])){
              stop("Apply averageIntraFragSpectra first")
            }

            return(average_xcms_grouped_msms(pa, "inter"))

          }
)


#' @title Using a purityA object, average and filter MS/MS spectra for each XCMS feature within
#' and across MS data files (ignoring intra and inter relationships)
#' @aliases averageAllFragSpectra
#' @description
#'
#' **General**
#'
#' Average and filter fragmentation spectra for each XCMS feature within and across MS data files (ignoring intra and inter relationships).
#'
#' The averaging is performed using hierarchical clustering of the m/z values of each peaks, where m/z values within a set ppm tolerance will be clustered. The clustered peaks are then averaged (or summed).
#'
#' The fragmentation can be filtered on the averaged spectra (with the arguments snr, rsd, minfrac, ra)
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
#'  * Fragmentation processing
#'    + (xset, pa) -> frag4feature -> filterFragSpectra -> **averageAllFragSpectra** -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#'
#' @param pa object; purityA object
#' @param cores numeric; Number of cores for multiprocessing
#' @param ppm numeric; ppm threshold to average across all scans (ignoring intra and inter relationships)
#' @param minnum numeric; minimum number of times peak is present across all fragmentation spectra (ignoring intra and inter relationships)
#' @param minfrac numeric;minimum ratio of the peak fraction (peak count / total peaks) across all (ignoring intra and inter relationships)
#' @param ra numeric; minimum relative abundance of the peak fraction across all (ignoring intra and inter relationships)
#' @param snr numeric;  minimum signal-to-noise of the peak across all (ignoring intra and inter relationships)
#' @param av character; type of averaging to use (median or mean)
#' @param sumi boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#'
#' @return Returns a purityA object (pa) with the following slots now with data
#'
#' * pa@@av_spectra: the average spectra is recorded here stored as a list. E.g. pa@@av_spectra$`1`$av_all would give the average spectra for grouped feature 1.
#' * pa@@av_all_params: The parameters used are recorded here
#'
#' Each spectra in the av_spectra list contains the following columns:
#'
#' * cl: id of clustered (averaged) peak
#' * mz: average m/z
#' * i: average intensity
#' * snr: average signal to noise ratio
#' * rsd: relative standard deviation
#' * count: number of clustered peaks
#' * total: total number of potential scans to be used for averaging
#' * inPurity: average precursor ion purity
#' * ra: average relative abundance
#' * frac: the fraction of clustered peaks (e.g. the count/total)
#' * snr_pass_flag: TRUE if snr threshold criteria met
#' * minfrac_pass_flag: TRUE if minfrac threshold criteria
#' * ra_pass_flag: TRUE if ra threshold criteria met
#' * pass_flag: TRUE if all threshold criteria met
#'
#'
#' @examples
#'
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' #xset <- xcms::group(xset)
#' #xset <- xcms::retcor(xset)
#' #xset <- xcms::group(xset)
#'
#' #pa  <- purityA(msmsPths, interpol = "linear")
#' #pa <- frag4feature(pa, xset)
#' #pa <- filterFragSpectra(pa)
#' pa <- readRDS(system.file("extdata", "tests", "purityA", "3_filterFragSpectra_pa.rds", package="msPurity"))
#' pa <- averageAllFragSpectra(pa)
#' @md
#' @export
setMethod(f="averageAllFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                                 av='median', sumi=TRUE, rmp=FALSE, cores=1
          ){

            pa@av_all_params$minfrac = minfrac
            pa@av_all_params$minnum = minnum
            pa@av_all_params$ppm = ppm
            pa@av_all_params$snr = snr
            pa@av_all_params$ra = ra

            pa@av_all_params$av_type = av
            pa@av_all_params$sumi = sumi
            pa@av_all_params$cores = cores
            pa@av_all_params$rmp = rmp

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

  colnames(df)[1:2] <- c('mz', 'i')
  df$index <-   rep(seq_along(grped_spectra), sapply(grped_spectra, nrow))

  if (!length(pa@filter_frag_params)==0){
    # if prior filtering performed only use those that have passed
    df<-df[df$pass_flag==1,]
  }

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
  if (!(av_level %in% c("intra", "inter", "all"))){
    stop("Incorrect av_level for averaging fragmentation spectra; use intra, inter or all")
  }


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
                                            sumi=pa@av_intra_params$sumi,
                                            minfrac=pa@av_intra_params$minfrac,
                                            snthr=pa@av_intra_params$snr,
                                            rathr=pa@av_intra_params$ra,
                                            rathr_pre= pa@av_intra_params$ra_pre,
                                            snrthr_pre= pa@av_intra_params$snr_pre,
                                            av_type=pa@av_intra_params$av_type)

    if (pa@av_intra_params$rmp){
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
                                 sumi=pa@av_inter_params$sumi,
                                 minfrac=pa@av_inter_params$minfrac,
                                 snthr=pa@av_inter_params$snr,
                                 rathr=pa@av_inter_params$ra,
                                 av_type=pa@av_inter_params$av_type
                               )

    if (pa@av_inter_params$rmp){
      av_inter <- av_inter[av_inter$pass_flag,]
    }

  } else if (av_level=="all") {

    # add additional column used later if filtering applied
    # Average the averaged spectra across everything (ignore intra and inter )
    av_all <- average_spectra(spectra_to_average,
                                cores=1,
                                ppm=pa@av_all_params$ppm,
                                minnum=pa@av_all_params$minnum,
                                sumi=pa@av_all_params$sumi,
                                minfrac=pa@av_all_params$minfrac,
                                snthr=pa@av_all_params$snr,
                                rathr=pa@av_all_params$ra,
                                rathr_pre= pa@av_all_params$ra_pre,
                                snrthr_pre= pa@av_all_params$snr_pre,
                                av_type=pa@av_all_params$av_type)
    if (pa@av_all_params$rmp){
      av_all <- av_all[av_all$pass_flag,]
    }
  } else {

    stop("Incorrect av_level for averaging fragmentation spectra; use intra, inter or all")

  }

  return(list('av_intra'=av_intra , 'av_inter'=av_inter, 'av_all'=av_all))
}



average_spectra <- function(spectra, indx='index', ppm, cores, minnum, sumi,
                            minfrac, snthr, snmeth='median', rathr, rathr_pre=NULL, snrthr_pre=NULL, av_type='median'){
  if (nrow(spectra)==0){
    return(NULL)
  }


  if (indx=='index'){
    # calculate metrics per scan (if using inter, the index will be sample and the snr and ra will
    # have already have been calculated
    # these will have already been calculated if filterFragSpectra has already been applied
    if ((!'snr' %in% colnames(spectra)) & (!'ra' %in% colnames(spectra))){
      spectra <- ddply(spectra, indx, set_snr_ra)
    }

  }


  if (!is.null(rathr_pre)){
    spectra <- spectra[spectra$ra>rathr_pre, ]
  }

  if (!is.null(snrthr_pre)){
    spectra <- spectra[spectra$snr>snrthr_pre, ]
  }

  if (nrow(spectra)==0){
    return(NULL)
  }

  # ensure ordered by mz
  spectra <- spectra[order(spectra$mz),]

  mz <- spectra$mz


  # Cluster the peaks togther
  spectra$cl <- clustering(mz, clustType = 'hc', cores = cores, ppm = ppm)

  averages <- plyr::ddply(spectra, ~ cl,
                          averageCluster, av=av_type, minnum=1,
                          missingV="ignore", totalScans=length(unique(spectra[,indx])), normTIC=FALSE,
                          sumI=sumi)




  averages$frac <- averages$count/averages$total

  averages$snr_pass_flag <- averages$snr > snthr

  averages$minnum_pass_flag <- averages$count >= minnum

  averages$minfrac_pass_flag <- averages$frac >= minfrac

  averages$ra_pass_flag <- averages$ra > rathr

  averages$pass_flag <- (averages$minfrac_pass_flag & averages$snr_pass_flag & averages$ra_pass_flag & averages$minnum_pass_flag)

  return(averages)
}


set_snr_ra <- function(x, snmeth='median'){

  if (snmeth=="median"){
    x$snr <- x$i/median(x$i)
  }else if(snmeth=="mean"){
    x$snr <- x$i/mean(x$i)
  }
  x$ra <- x$i/max(x$i)*100

   return(x)
}

