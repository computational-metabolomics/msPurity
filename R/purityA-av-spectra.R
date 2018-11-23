#' @title Average fragemntation spectra across XCMS features
#'
#' @description
#'
#' Average fragemntation spectra across XCMS features.
#'
#'
#' @aliases frag4feature
#'
#' @param pa object; purityA object
#' @param xset object; XCMS object derived from the same files as the puritydf
#' @param ppm numeric; ppm tolerance between precursor mz and feature mz
#' @param plim numeric; min purity of precursor to be included
#' @param intense boolean; If the most intense precursor or the centered precursor is used
#' @param use_group boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If retention time correction has been used in XCMS set this to TRUE
#' @param create_db boolean; SQLite database will be created of the results
#' @param db_name character; If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @param out_dir character; Path where database will be created
#' @param grp_peaklist dataframe [optional]; Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @return purityA object with slots for fragmentation-XCMS links
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
                                minnum_inter=1, ppm_intra=5, ppm_inter=5, av='median', sum_i=TRUE,
                                purity_thres=0.5){

            pa@av_params$minfrac_intra = minfrac_intra
            pa@av_params$minfrac_inter = minfrac_inter
            pa@av_params$minnum_intra = minnum_intra
            pa@av_params$minnum_inter = minnum_inter
            pa@av_params$ppm_intra = ppm_intra
            pa@av_params$ppm_inter = ppm_inter
            pa@av_params$av = av
            pa@av_params$sum_i = sum_i
            pa@av_params$purity_thres = purity_thres

            return(average_xcms_grouped_msms_all(pa))

          }
)

average_xcms_grouped_msms_all <- function(pa){

  av_spectra <- plyr::alply(names(pa@grped_ms2), 1, average_xcms_grouped_msms_indiv, pa=pa)

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

  spectra_to_average <- merge(df, grped_info[, c('grpid', 'sample', 'cid', 'index')], by = "index")


  ##############################################################################
  # Performing averaging
  ##############################################################################
  # clustering requires data to be in order of mz
  spectra_to_average <- spectra_to_average[order(spectra_to_average$mz),]

  # Average by sample (file)
  av_sample <- plyr::ddply(spectra_to_average, ~sample, average_spectra)

  # Average the averaged spectra across files
  av_all <- average_spectra(av_sample, indx='sample')

  # add additional column used later if filtering applied

  av_all$keep = TRUE

  return(av_all)
}

average_spectra <- function(spectra, indx='index', ppm, cores, minnum){
  mz <- spectra$mz

  # Cluster the peaks togther
  spectra$cl <- clustering(mz, clustType = 'hc', cores = 1, ppm = 5)

  averages <- plyr::ddply(spectra, ~ cl,
                          averageCluster, av="median", minnum=1,
                          missingV="ignore", totalScans=length(unique(spectra[,indx])), normTIC=FALSE)
}


