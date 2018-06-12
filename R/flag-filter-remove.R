#' @title Flag and remove unwanted peaks
#'
#' @description
#' On an xcmsSet object, filter flag and remove unwanted peaks. When the peaks are removed, the
#' the xcmsSet object can be regrouped using xcms::group. The function then checks if any blank
#' peaks are still present and the process is repeated.
#'
#' The output is a list of the updated xcmsSet object, grouped peaklist and the blank removed peaks
#'
#' @param xset object; xcmsSet object
#' @param pol str; polarity (just used for naming purpose for files being saved) [positive, negative, NA]
#' @param rsd_i_blank numeric; RSD threshold for the blank
#' @param minfrac_blank numeric; minimum fraction of files for features needed for the blank
#' @param rsd_rt_blank numeric; RSD threshold for the RT of the blank
#' @param ithres_blank numeric; Intensity threshold for the blank
#' @param s2b numeric; fold change (sample/blank) needed for sample peak to be allowed. e.g.
#'                     if s2b set to 10 and the recorded sample 'intensity' value was 100 and blank = 10.
#'                     1000/10 = 100 so sample has fold change higher than the threshold and the peak is not considered a blank
#' @param ref.class str; A string representing the class that will be used for the blank.
#' @param egauss_thr numeric; Threshold for filtering out non gaussian shaped peaks. Note this only works
#'                            if the verbose option was set for XCMS;
#' @param rsd_i_sample numeric; RSD threshold for the sample
#' @param minfrac_sample numeric; minimum fraction of files for features needed for the sample
#' @param rsd_rt_sample numeric; RSD threshold for the RT of the sample
#' @param ithres_sample numeric; Intensity threshold for the sample
#' @param grp_rm_ids vector; vector of grouped_xcms peaks to remove (coresponds to the row from xcms::group output)
#'
#' @param remove_spectra bool; TRUE if flagged spectra is to be removed
#'
#' @param minfrac_xcms numeric; minfrac for xcms  grouping
#' @param mzwid numeric; xcms grouping parameter
#' @param bw numeric; xcms grouping parameter
#' @param out_dir str; out directory
#' @param temp_save boolean; Assign True if files for each step saved (for testing purpsoses)
#'
#' @return list(xset, grp_peaklist, removed_peaks)
#' @examples
#'
#' msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)
#' xset <- xcms::xcmsSet(msPths, BPPARAM = BiocParallel::SnowParam())
#' xset@phenoData[,1] <- c('blank', 'blank', 'sample', 'sample')
#' xset <- xcms::group(xset)
#' fr = flag_remove(xset)
#' @export
flag_remove <- function(xset, pol=NA, rsd_i_blank=NA, minfrac_blank=0.5,
                        rsd_rt_blank=NA, ithres_blank=NA, s2b=10, ref.class='blank',
                        egauss_thr=NA, rsd_i_sample=NA, minfrac_sample=0.7,
                        rsd_rt_sample=NA, ithres_sample=NA, minfrac_xcms=0.7,
                        mzwid=0.025, bw=5, out_dir='.', temp_save=FALSE, remove_spectra=TRUE, grp_rm_ids=NA){

  ################################
  # Get enriched peaklist
  ################################
  # Get a peaklist with additional information e.g. RSD, coverage, 'valid' columns
  grp_peaklist <- create_enriched_peaklist(xset, rsd_i_blank, minfrac_blank,
                                           rsd_rt_blank, ithres_blank, s2b, ref.class,
                                           rsd_i_sample, minfrac_sample, rsd_rt_sample, ithres_sample)

  # Save the original peaklist for reference
  grp_peaklist_orig <- grp_peaklist
  grp_peaklist <- get_full_peak_width(grp_peaklist_orig, xset)

  ##################################
  # Remove blank and invalid peaks
  ##################################
  xset.count <- 1

  if(remove_spectra){
    print('########## REMOVING FLAGGED PEAKS ###########')

    # Remove any grouped peaks that have been specifically selected
    if (!anyNA(grp_rm_ids)){
      grp_peaklist <- cbind(grp_peaklist, 'grp_remove'=rep(0, nrow(grp_peaklist)))
      grp_peaklist[grp_rm_ids,'grp_remove'] <- 1
      grp_peaklist[grp_rm_ids,'all_sample_valid'] <- 0
    }

    removed_peaks = data.frame()
    while(sum(grp_peaklist[,paste(ref.class, '_valid', sep='')])>0){
      # Remove blank peaks and invalid sample peaks from xcms object (then regroup)
      print(paste('xset', xset.count,sep=''))

      rms <- remove_spectra(xset, grp_peaklist, rclass=ref.class, rm_peak_out=TRUE)
      xset <- rms[[1]]
      remPeaksTemp <- data.frame(rms[[2]])
      remPeaksTemp$xcmsIDx <- xset.count


      if(xset.count==1){
        removed_peaks <- remPeaksTemp
      }else{
        removed_peaks <- rbind(remPeaksTemp, removed_peaks)
      }
      # need to regroup the xcms object
      xset <- xcms::group(xset, bw = bw, mzwid = mzwid, minfrac = minfrac_xcms)

      # Update the grouped_peaklist
      grp_peaklist <- create_enriched_peaklist(xset, rsd_i_blank, minfrac_blank,
                                               rsd_rt_blank, ithres_blank, s2b, ref.class,
                                               rsd_i_sample, minfrac_sample, rsd_rt_sample, ithres_sample)

      grp_peaklist <- get_full_peak_width(grp_peaklist, xsa = xset)

      xset.count <- xset.count + 1

      if(temp_save){
        saveRDS(xset, file.path(out_dir, paste(pol, 'xset', xset.count, '.rds', sep='_')))
      }

    }

    # Remove all individual where the egauss (RSM of the Gaussian fit to curve) is too high
    if (!is.na(egauss_thr)){
      xset@peaks <- xset@peaks[(xset@peaks[,'egauss']<=egauss_thr) & (!is.na(xset@peaks[,'egauss'])) ,]
      # need to regroup the xcms object
      xset <- xcms::group(xset, bw = bw, mzwid = mzwid, minfrac = minfrac_xcms)
      grp_peaklist <- create_enriched_peaklist(xset, rsd_i_blank, minfrac_blank,
                                               rsd_rt_blank, ithres_blank, s2b, ref.class,
                                               rsd_i_sample, minfrac_sample, rsd_rt_sample, ithres_sample)
      grp_peaklist <- get_full_peak_width(grp_peaklist, xset)
    }
    if (nrow(removed_peaks)>0){
      write.csv(removed_peaks, file.path(out_dir, paste('removed_peaks_', pol, '.csv', sep='')))
    }

  }else{
    removed_peaks = NA
  }

  if(temp_save){
    write.csv(grp_peaklist, file.path(out_dir, paste(pol, 'grp_peaklist_blanksRemoved.csv', sep='_')))
  }

  return(list(xset, grp_peaklist, removed_peaks))


}



