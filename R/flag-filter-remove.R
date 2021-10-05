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



#' @title Flag and remove unwanted peaks
#'
#' @description
#' Filter, flag and remove unwanted peaks from xcms object (xcmsObj) of class XCMSnExp, xcmsSet or xsAnnotate.
#' When the peaks are removed, the xcmsObj object can be regrouped (originally using xcms::group, now using xcms::groupChromPeaks).
#' The function then checks if any blank peaks are still present and the process is repeated.
#'
#' The output is a list object containing: 1) the updated xcms object, 2) the grouped peaklist and 3) the blank removed peaks
#'
#' @param xcmsObj object; XCMSnExp, xcmsSet or xsAnnotate object
#' @param pol str; polarity (just used for naming purpose for files being saved) \[positive, negative, NA\]
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
#' @param grp_rm_ids vector; vector of grouped_xcms peaks to remove (corresponds to the row from xcms::group output)
#'
#' @param remove_spectra_bool bool; TRUE if flagged spectra is to be removed
#'
#' @param minfrac_xcms numeric; minfrac for xcms  grouping
#' @param mzwid numeric; xcms grouping parameter (corresponds to variable 'binSize' in XCMS3)
#' @param bw numeric; xcms grouping parameter
#' @param out_dir str; out directory
#' @param temp_save boolean; Assign True if files for each step saved (for testing purpsoses)
#' @param xset object, DEPRECATED; xcmsSet object
#'
#' @return list(xset, grp_peaklist, removed_peaks)
#' @examples
#' library(xcms)
#' library(MSnbase)
#' library(magrittr)
#' #read in files and data
#' msPths <-list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)
#' ms_data = readMSData(msPths, mode = 'onDisk', msLevel. = 1)
#'
#' #subset the data to focus on retention times 30-90 seconds and m/z values between 100 and 200 m/z.
#' rtr = c(30, 90)
#' mzr = c(100, 200)
#' ms_data = ms_data %>%  filterRt(rt = rtr) %>%  filterMz(mz = mzr)
#'
#' ##### perform feature detection in individual files
#' cwp <- CentWaveParam(snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30))
#' xcmsObj <- findChromPeaks(ms_data, param = cwp)
#' xcmsObj@phenoData@data$class = c('blank', 'blank', 'sample', 'sample')
#' xcmsObj@phenoData@varMetadata = data.frame('labelDescription' = 'sampleNames', 'class')
#' pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = 0, bw = 5, binSize = 0.017)
#' xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)
#'
#' #### flag, filter and remove peaks, returning an updated xcmsObj (XCMSnExp or xcmsSet class), grouped_peaklist (data.frame) and removed_peaks (data.frame)
#' fr <- flag_remove(xcmsObj)
#'
#' ##### load from existing data
#' xcmsObj = readRDS(system.file("extdata", "tests", "purityA", "10_input_filterflagremove.rds", package="msPurity"))
#'
#'
#'
#'
#' @export
flag_remove <- function(xcmsObj, pol=NA, rsd_i_blank=NA, minfrac_blank=0.5,
                        rsd_rt_blank=NA, ithres_blank=NA, s2b=10, ref.class='blank',
                        egauss_thr=NA, rsd_i_sample=NA, minfrac_sample=0.7,
                        rsd_rt_sample=NA, ithres_sample=NA, minfrac_xcms=0.7,
                        mzwid=0.017, bw=5, out_dir='.', temp_save=FALSE, remove_spectra_bool=TRUE,
                        grp_rm_ids=NA, xset=NA){

  if(!is.na(xset)){
    xcmsObj = xset
  }

  if(is(xcmsObj, 'XCMSnExp')){
    XCMSnExp_bool = TRUE
  }else if(is(xcmsObj, 'xcmsSet')){
    XCMSnExp_bool = FALSE
  }else if(is(xcmsObj, 'xsAnnotation')){
    XCMSnExp_bool = FALSE
    xcmsObj = xcmsObj@xcmsSet
  }else{
    stop('unrecognised class for xcmsObj object')
  }

  ################################
  # Get enriched peaklist
  ################################
  # Get a peaklist with additional information e.g. RSD, coverage, 'valid' columns
  grp_peaklist <- create_enriched_peaklist(xcmsObj, rsd_i_blank, minfrac_blank,
                                           rsd_rt_blank, ithres_blank, s2b, ref.class,
                                           rsd_i_sample, minfrac_sample, rsd_rt_sample, ithres_sample,
                                           XCMSnExp_bool)

  # Save the original peaklist for reference
  grp_peaklist_orig <- grp_peaklist
  grp_peaklist <- get_full_peak_width(grp_peaklist_orig, xcmsObj)

  ##################################
  # Remove blank and invalid peaks
  ##################################
  obj.count <- 1

  if(remove_spectra_bool){
    message('########## REMOVING FLAGGED PEAKS ###########')

    # Remove any grouped peaks that have been specifically selected
    if (!anyNA(grp_rm_ids)){
      grp_peaklist <- cbind(grp_peaklist, 'grp_remove'=rep(0, nrow(grp_peaklist)))
      grp_peaklist[grp_rm_ids,'grp_remove'] <- 1
      grp_peaklist[grp_rm_ids,'all_sample_valid'] <- 0
    }

    removed_peaks = data.frame()
    while(sum(grp_peaklist[,paste(ref.class, '_valid', sep='')])>0){
      # Remove blank peaks and invalid sample peaks from xcms object (then regroup)
      message(paste('xcmsObj', obj.count,sep=''))

      rms <- remove_spectra(xcmsObj, grp_peaklist, rclass=ref.class, rm_peak_out=TRUE, XCMSnExp_bool = XCMSnExp_bool)
      xcmsObj <- rms[[1]]
      remPeaksTemp <- data.frame(rms[[2]])
      remPeaksTemp$xcmsIDx <- obj.count


      if(obj.count==1){
        removed_peaks <- remPeaksTemp
      }else{
        removed_peaks <- rbind(remPeaksTemp, removed_peaks)
      }


      # need to regroup the xcms object
      if(XCMSnExp_bool){
        pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = minfrac_xcms, bw = bw, binSize = mzwid)
        xcmsObj <- xcms::groupChromPeaks(xcmsObj, param = pdp)
      }else{
        xcmsObj <- xcms::group(xcmsObj, bw = bw, mzwid = mzwid, minfrac = minfrac_xcms)
      }

      # Update the grouped_peaklist
      grp_peaklist <- create_enriched_peaklist(xcmsObj, rsd_i_blank, minfrac_blank,
                                               rsd_rt_blank, ithres_blank, s2b, ref.class,
                                               rsd_i_sample, minfrac_sample, rsd_rt_sample, ithres_sample,
                                               XCMSnExp_bool = XCMSnExp_bool)

      grp_peaklist <- get_full_peak_width(grp_peaklist, xcmsObj = xcmsObj)

      obj.count <- obj.count + 1

      if(temp_save){
        saveRDS(xcmsObj, file.path(out_dir, paste(pol, 'xcmsnexp', obj.count, '.rds', sep='_')))
      }

    }

    # Remove all individual where the egauss (RSM of the Gaussian fit to curve) is too high
    if (!is.na(egauss_thr)){
      if(XCMSnExp_bool){
        if('egauss' %in% colnames(xcms::chromPeaks(xcmsObj))){
          chromPeaks(xcmsObj) <- chromPeaks(xcmsObj)[(chromPeaks(xcmsObj)[,'egauss']<=egauss_thr) & (!is.na(chromPeaks(xcmsObj)[,'egauss'])) ,]
          # need to regroup the xcms object
          pdp <- PeakDensityParam(sampleGroups = xcmsObj@phenoData@data$class, minFraction = minfrac_xcms, bw = bw, binSize = mzwid)
          xcmsObj <- xcms::groupChromPeaks(xcmsObj, param = pdp)
        }else{
          message('"egauss_thr" specified but no column "egauss" in xcms::chromPeaks(xcmsObj) - continuing without filtering based on egauss' )
        }

      }else{
        if('egauss' %in% colnames(xcmsObj@peaks)){
          xcmsObj@peaks <- xcmsObj@peaks[(xcmsObj@peaks[,'egauss']<=egauss_thr) & (!is.na(xcmsObj@peaks[,'egauss'])) ,]
          xcmsObj <- xcms::group(xcmsObj, bw = bw, mzwid = mzwid, minfrac = minfrac_xcms)
        }else{
          message('"egauss_thr" specified but no column "egauss" in xcmsObj@peaks (xcmsSet) - continuing without filtering based on egauss' )
        }
      }

      grp_peaklist <- create_enriched_peaklist(xcmsObj, rsd_i_blank, minfrac_blank,
                                               rsd_rt_blank, ithres_blank, s2b, ref.class,
                                               rsd_i_sample, minfrac_sample, rsd_rt_sample, ithres_sample,
                                               XCMSnExp_bool = XCMSnExp_bool)
      grp_peaklist <- get_full_peak_width(grp_peaklist, xcmsObj = xcmsObj)
    }

    if (nrow(removed_peaks)>0){
      if (temp_save){
        write.csv(removed_peaks, file.path(out_dir, paste('removed_peaks_', pol, '.csv', sep='')))
      }

    }

  }else{
    removed_peaks = NA
  }

  if(temp_save){
    write.csv(grp_peaklist, file.path(out_dir, paste(pol, 'grp_peaklist_blanksRemoved.csv', sep='_')))
  }


  grp_peaklist <- data.frame(cbind(grp_peaklist, 'grp_names'=xcms::groupnames(xcmsObj)))

  return(list('xcmsObj'= xcmsObj, 'grp_peaklist'= grp_peaklist, 'removed_peaks'= removed_peaks))


}


create_enriched_peaklist <- function(xcmsObj,rsd_i_blank, minfrac_blank,rsd_rt_blank,
                                     ithres_blank, s2b, ref.class, rsd_i_sample,
                                     minfrac_sample, rsd_rt_sample, ithres_sample, XCMSnExp_bool){
  ###########################################
  # Create enriched peaklist
  ###########################################
  # Description:
  #   Create a peaklist of the grouped xcms features but with additional information. The peaklist
  #   is flagged with valid and invalid peaks based on the params
  #
  # Args:
  #   xset: xcmsSet grouped object
  #   ...
  #
  # Returns:
  #   enriched peaklist
  #

  # Add all the summary stats for each class
  grp_peaklist <- sum_calc_peaklist(xcmsObj, XCMSnExp_bool)




  ##################################################
  # THIS IS FOR THE BLANK
  # Flag valid peaks to be used for the blank
  #################################################

  if(XCMSnExp_bool){
    classes <- as.character(unique(xcms::phenoData(xcmsObj)@data$class))
  }else{
    classes <- as.character(unique(xcms::phenoData(xcmsObj)$class))
  }

  grp_peaklist <- flag_peaks(peaklist = grp_peaklist,
                             RSD_I_filter = rsd_i_blank,
                             minfrac_filter = minfrac_blank,
                             RSD_RT_filter = rsd_rt_blank,
                             i_thre_filter = ithres_blank,
                             s2b = s2b,
                             fclass = ref.class,
                             classes = classes)


  ##################################################
  # THIS IS FOR ALL OTHER NON BLANK CLASSES
  #################################################
  # A Loop through the sample classes updating the associated
  # 'valid' column
  valids = ''
  for(c in classes){
    if(c==ref.class){
      next
    }
    grp_peaklist <- flag_peaks(peaklist = grp_peaklist,
                               RSD_I_filter = rsd_i_sample,
                               minfrac_filter = minfrac_sample,
                               RSD_RT_filter = rsd_rt_sample,
                               i_thre_filter = ithres_sample,
                               s2b = 0,
                               fclass = c,
                               classes = classes)

    valids <- c(valids, paste(c, 'valid', sep='_'))
  }
  valids <- valids[-1]

  # Check if any of the peaks are valid for each sample. If they are we create a
  # summary flag to keep this peak (unless we remove it later at the subtraction stage)
  temp <- apply(as.array(grp_peaklist[,valids]), 1, sum)

  grp_peaklist <- cbind(grp_peaklist, 'all_sample_valid'=rep(0, nrow(grp_peaklist)))
  grp_peaklist[temp>0,][,'all_sample_valid']=1

  return(grp_peaklist)


}

sum_calc_peaklist <- function(xcmsObj, XCMSnExp_bool){
  ###########################################
  # Calculates summary information for XCMS peaks
  ###########################################
  # Description:
  #   Create a peaklist of the grouped xcms features but with additional information
  #   including:
  #         * RSD for intensity [blank, samples]
  #         * RSD for RT [blank, samples, all]
  #         * median intensity [blank, samples]
  #         * coverage [blank, samples]
  # Args:
  #   xset = xcmsSet object (has to have been grouped prior)
  #
  # Returns:
  #   updated peaklist

  # Get the group peaklist
  if(XCMSnExp_bool){
    grp_peaklist <- data.frame(xcms::featureDefinitions(xcmsObj))
  }else{
    grp_peaklist <- data.frame(xcms::groups(xcmsObj))
    grp_peaklist[,'ms_level'] = xcms::mslevel(xcmsObj)
  }



  rownames(grp_peaklist) <- seq(1, nrow(grp_peaklist))

  # Get the intensity values
  if(XCMSnExp_bool == TRUE){
    indiv_peaks <- xcms::featureValues(xcmsObj, 'medret', intensity = 'into', value = 'into')
  }else{
    indiv_peaks <- xcms::groupval(xcmsObj, 'medret', intensity="into", value='into')
  }

  grp_peaklist <- cbind(grp_peaklist, indiv_peaks)

  # Calculate summary stats for each class
  classes = as(phenoData(xcmsObj), 'data.frame')
  classes$fileid  <- rownames(classes)
  classes$class <- as.character(classes$class)

  # Create a matrix for each of the summary stats that use intensity
  cnt = 1
  colnm <- length(unique(classes$class))
  med_i_m <- matrix(nrow = nrow(grp_peaklist), ncol=colnm)
  colnames(med_i_m) <- 1:colnm
  rsd_i_m <- med_i_m
  coverage_m <- med_i_m

  # Intensity: RSD and median (and coverage)

  for (c in unique(classes$class)){

    x <- classes[classes$class==c,]
    files <- x$fileid

    coverage_m[,cnt] <- grp_peaklist[,c]/length(files)
    colnames(coverage_m)[cnt] <- paste(c, "coverage", sep="_")

    subseti <- grp_peaklist[ , which(colnames(grp_peaklist) %in% files)]
    rsdi <- apply(subseti, 1, rsd)
    mediani <- apply(subseti, 1, median, na.rm=TRUE)

    med_i_m[,cnt] <- mediani
    colnames(med_i_m)[cnt] <- paste(c, "median_I", sep="_")

    rsd_i_m[,cnt] <- rsdi
    colnames(rsd_i_m)[cnt] <- paste(c, "RSD_I", sep="_")
    cnt = cnt+1
  }

  # Combine the summarys stats with the group peaklist
  grp_peaklist <- cbind(grp_peaklist, med_i_m, rsd_i_m, coverage_m)


  # Get the retention time values

  if(XCMSnExp_bool){
    grouprt <- xcms::featureValues(xcmsObj, method='medret', value = 'rt')
  }else{
    grouprt <- xcms::groupval(xcmsObj, method='medret', value = 'rt')
  }

  # retention time RSD
  for (c in unique(classes$class)){

    x <- classes[classes$class==c,]
    files <- x$fileid

    # rt RSD
    subsetT <- grouprt[ , which(colnames(grouprt) %in% files)]
    rsdT <- apply(subsetT, 1, rsd)
    grp_peaklist <- cbind(grp_peaklist, rsdT)
    colnames(grp_peaklist)[ncol(grp_peaklist)] = paste(c, "RSD_RT", sep="_")
  }

  # retention time rsd of all the peaks (see to what extent the times have been aligned)
  rsd_all_RT <- apply(grouprt, 1, rsd)
  grp_peaklist <- cbind(grp_peaklist, rsd_all_RT, 'grpid'=seq(1, nrow(grp_peaklist)))

  # Add a column as to whether the feature is considered valid (set to 0) now. We
  # can update this later
  valid_col <- rep(0, nrow(grp_peaklist))
  for (c in unique(classes$class)){
    grp_peaklist <- cbind(grp_peaklist, valid_col)
    colnames(grp_peaklist)[ncol(grp_peaklist)] <- paste(c, "valid", sep="_")
  }

  return(grp_peaklist)
}

flag_peaks <- function(peaklist, RSD_I_filter, minfrac_filter, RSD_RT_filter, i_thre_filter, fclass, s2b, classes){
  ###########################################
  # Flag peaks that pass the 'valid' criteria
  ###########################################
  # Description:
  #   Flag a peak list, check that a value passes RSD, minfrac, intensity threshold. Also checks the sample2blank
  #   intensity threshold
  #
  # Returns:
  #   updated peaklist

  rsd_i_nm <- paste(fclass, "RSD", "I", sep="_")
  rsd_rt_nm <- paste(fclass, "RSD", "RT", sep="_")
  median_i_nm <- paste(fclass, "median", "I", sep="_")
  coverage_nm <- paste(fclass, "coverage", sep="_")
  valid_col_nm <- paste(fclass, "valid", sep="_")


  if (anyNA(minfrac_filter)){
    min_samp_bool <- rep(TRUE, nrow(peaklist))
  }else{
    min_samp_bool <- peaklist[,coverage_nm]>=minfrac_filter
  }

  if (anyNA(RSD_I_filter)){
    rsd_i_bool_thr <- rep(TRUE, nrow(peaklist))
    rsd_i_bool_na <- rep(TRUE, nrow(peaklist))
  }else{
    rsd_i_bool_thr <- peaklist[,rsd_i_nm]<RSD_I_filter
    rsd_i_bool_na <- !is.na(peaklist[,rsd_i_nm])
  }


  if (anyNA(RSD_RT_filter)){
    rsd_rt_bool_thr  <- rep(TRUE, nrow(peaklist))
    rsd_rt_bool_na <- rep(TRUE, nrow(peaklist))
  }else{
    rsd_rt_bool_thr <- peaklist[,rsd_rt_nm ]<RSD_RT_filter
    rsd_rt_bool_na <- !is.na(peaklist[,rsd_rt_nm])
  }

  if (anyNA(i_thre_filter)){
    i_thre_bool <- rep(TRUE, nrow(peaklist))
    i_thre_bool_na <- rep(TRUE, nrow(peaklist))
  }else{
    i_thre_bool <- peaklist[,median_i_nm]>i_thre_filter
    i_thre_bool_na <- !is.na(peaklist[,median_i_nm ])
  }


  if(s2b>0){
    # FOR THE BLANK
    # if the intensity is n times less than any other sample class intensity
    # we do not use (for blank subtraction)
    blankI <- peaklist[,paste(fclass, 'median_I', sep="_")]
    #classes <- as.character(unique(xset@phenoData$class))
    classes <- classes[!classes==fclass]
    s2b_sum <- sapply(classes, function(c){
      targetI <- peaklist[,paste(c, 'median_I', sep="_")]
      return((targetI/blankI > s2b) & (!is.na(targetI/blankI)))
    })

    # if sums to 0 it means that the sample intensity is NOT n times greater than the s2b value
    s2b_check <- apply(s2b_sum, 1, sum)

    # valid blank peaks (i.e. sample is not 'n' times greater than the blank for any of the sample classes)
    s2b_bool = s2b_check==0

  }else{
    # NOT THE BLANK
    s2b_bool = rep(TRUE, nrow(peaklist))
  }


  sub <- peaklist[min_samp_bool &
                    rsd_i_bool_na & rsd_i_bool_thr &
                    rsd_rt_bool_na & rsd_rt_bool_thr &
                    i_thre_bool_na & i_thre_bool &
                    s2b_bool,,drop=FALSE
                  ]


  if(is.vector(sub)){
    peaklist[sub['grpid'],][valid_col_nm] = 1
  }else if(nrow(sub)==1){
    peaklist[rownames(sub),][valid_col_nm] = 1
  }else if(nrow(sub)<1){
    message(paste('Warning (ignore this statement if blank peaks have been remove). No peaks for ', fclass))
  }else{
    peaklist[rownames(sub),][,valid_col_nm] = 1
  }

  return(peaklist)

}

rsd <- function(x){ (sd(x,na.rm = TRUE)/mean(x, na.rm=TRUE))*100 }

remove_spectra <- function(xcmsObj, peaklist, rclass, rm_peak_out=FALSE, XCMSnExp_bool){
  # remove the valid peaks of one selected class from an xcms object

  peaklist = get_full_peak_width(peaklist, xcmsObj)

  valid_blank_bool <- peaklist[,paste(rclass, 'valid', sep='_')]==1
  invalid_sample_bool <- peaklist[,'all_sample_valid']==0


  removed_peaks <- peaklist[(valid_blank_bool | invalid_sample_bool), ,drop=FALSE]

  grp_ids_rm <- unlist(removed_peaks[,'grpid'])

  if(XCMSnExp_bool){
    peak_ids_rm <- unlist(xcms::featureDefinitions(xcmsObj)$peakidx[grp_ids_rm])
    xcms::chromPeaks(xcmsObj) <- xcms::chromPeaks(xcmsObj)[-peak_ids_rm,]
  }else{
    peak_ids_rm <- unlist(xcmsObj@groupidx[grp_ids_rm])
    xcmsObj@peaks <- xcmsObj@peaks[-peak_ids_rm,]
  }

  if( rm_peak_out){
    return(list(xcmsObj, removed_peaks))
  }else{
    return(xcmsObj)
  }

}

get_full_peak_width <- function(peaklist, xcmsObj){
  ###########################################
  # Get full peak width
  ###########################################
  # Args:
  #   peaklist: the peak list generated from either XCMS or CAMERA.
  #              Use the CAMERA peak list for this piplein
  #   xcmsObj: xcms object of class XCMSnExp, xcmsSet or xsAnnotate
  #
  # Returns:
  #   An updated peaklist with the full retention window ranges (and full mz ranges)
  #
  # See also:
  #   full_minmax, getpeaks, ldply (from the plyr library)

  message("Getting full peak widths")
  # Get 'peaks' (xcms features) from the xcmsSet object stored
  # in the camera annotation object

  message("Get 'individual' peaks from camera-xcms object")

  if(is(xcmsObj,'XCMSnExp')){
    XCMSnExp_bool = TRUE
  }else if(is(xcmsObj, 'xcmsSet')){
    XCMSnExp_bool = FALSE
  }else if(is(xcmsObj, 'xsAnnotate')){
    XCMSnExp_bool = FALSE
    xcmsObj = xcmsObj@xcmsSet
  }else{
    stop('unrecognised class for "xcmsObj", should be either "XCMSnExp", "xcmsSet", "xsAnnotate"')

  }

  if(XCMSnExp_bool && (is(xcmsObj,'XCMSnExp'))){
    rt.min = xcms::featureValues(xcmsObj, method = "medret", value = "rtmin", intensity = "into")
    rt.max = xcms::featureValues(xcmsObj, method = "medret", value = "rtmax", intensity = "into")
    mz.min = xcms::featureValues(xcmsObj, method = "medret", value = "mzmin", intensity = "into")
    mz.max = xcms::featureValues(xcmsObj, method = "medret", value = "mzmax", intensity = "into")
  }else if (XCMSnExp_bool==FALSE && (is(xcmsObj, 'xcmsSet'))){
    rt.min = xcms::groupval(xcmsObj, method = "medret", value = "rtmin", intensity = "into")
    rt.max = xcms::groupval(xcmsObj, method = "medret", value = "rtmax", intensity = "into")
    mz.min = xcms::groupval(xcmsObj, method = "medret", value = "mzmin", intensity = "into")
    mz.max = xcms::groupval(xcmsObj, method = "medret", value = "mzmax", intensity = "into")
  }

  rt.min = apply(rt.min, 1, min, na.rm = TRUE)
  rt.max = apply(rt.max, 1, max, na.rm = TRUE)
  mz.min = apply(mz.min, 1, min, na.rm = TRUE)
  mz.max = apply(mz.max, 1, max, na.rm = TRUE)


  peaklist_full = cbind(peaklist, "mzmin_full" = mz.min, "mzmax_full" = mz.max, "rtmin_full" = rt.min, "rtmax_full" = rt.max)
  return(peaklist_full)

}

full_minmax <- function(grouplist,peaks,groupdf){
  # Use the the plyr package to go through a list() of peaks and output a
  # dataframe of actual rtmin, rtmax, mzmin and mzmax
  full_info <- plyr::ldply(grouplist,getpeaks,peaks=peaks)

  # Update users original dataframe
  groupdf$mzmin_full <- full_info$mzmin
  groupdf$mzmax_full <- full_info$mzmax
  groupdf$rtmin_full <- full_info$rtmin
  groupdf$rtmax_full <- full_info$rtmax

  # return the updated dataframe
  return(groupdf)

}

getpeaks <- function(row,peaks) {
  # In case no peaks found. Should be worried if this happens though
  if (is.null(nrow(row)) || nrow(row)==0 ) {
    return(c(mzmin = NA, mzmax = NA, rtmin = NA, rtmax = NA))
  }

  # get the associated peaks for the group
  fullpeak <- peaks[row$peakID,]

  # output a vector of the actual min max values
  return(c(mzmin = min(fullpeak[,"mzmin"]), mzmax = max(fullpeak[,"mzmax"]),
    rtmin = min(fullpeak[,"rtmin"]), rtmax = max(fullpeak[,"rtmax"])))
}


