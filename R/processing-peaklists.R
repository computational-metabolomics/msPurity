create_enriched_peaklist <- function(xset,rsd_i_blank, minfrac_blank,rsd_rt_blank,
                                     ithres_blank, s2b, ref.class, rsd_i_sample,
                                     minfrac_sample, rsd_rt_sample, ithres_sample){
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
  grp_peaklist <- sum_calc_peaklist(xset)

  ##################################################
  # THIS IS FOR THE BLANK
  # Flag peaks valid peaks to be used for the blank
  #################################################
  grp_peaklist <- flag_peaks(peaklist = grp_peaklist,
                             RSD_I_filter = rsd_i_blank,
                             minfrac_filter = minfrac_blank,
                             RSD_RT_filter = rsd_rt_blank,
                             i_thre_filter = ithres_blank,
                             s2b = s2b,
                             fclass = ref.class,
                             xset=xset)

  # Filter peaks for the other classes
  classes <- as.character(unique(xset@phenoData$class))

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
                               xset=NA)

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

sum_calc_peaklist <- function(xset){
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
  grp_peaklist <- xcms::groups(xset)
  rownames(grp_peaklist) <- seq(1, nrow(grp_peaklist))

  # Get the itensity values
  indiv_peaks <- groupval(xset, 'medret', intensity="into", value='into')
  grp_peaklist <- cbind(grp_peaklist, indiv_peaks)


  # Calculate summary stats for each class
  classes <- xset@phenoData
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
    mediani <- apply(subseti, 1, median, na.rm=T)

    med_i_m[,cnt] <- mediani
    colnames(med_i_m)[cnt] <- paste(c, "median_I", sep="_")

    rsd_i_m[,cnt] <- rsdi
    colnames(rsd_i_m)[cnt] <- paste(c, "RSD_I", sep="_")
    cnt = cnt+1
  }

  # Combine the summarys stats with the group peaklist
  grp_peaklist <- cbind(grp_peaklist, med_i_m, rsd_i_m, coverage_m)


  # Get the retention time values
  grouprt <- groupval(xset, method='medret', value = 'rt')

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

flag_peaks <- function(peaklist, RSD_I_filter, minfrac_filter, RSD_RT_filter, i_thre_filter, fclass, s2b, xset){
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
    classes <- as.character(unique(xset@phenoData$class))
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
    print(paste('Warning (ignore this statement if blank peaks have been remove). No peaks for ', fclass))
  }else{
    peaklist[rownames(sub),][,valid_col_nm] = 1
  }

  return(peaklist)

}

rsd <- function(x){ (sd(x,na.rm = T)/mean(x, na.rm=T))*100 }

remove_spectra <- function(xset, peaklist, rclass, rm_peak_out=FALSE){
  # remove the valid peaks of one selected class from an xcms object

  peaklist = get_full_peak_width(peaklist, xset)

  valid_blank_bool <- peaklist[,paste(rclass, 'valid', sep='_')]==1
  invalid_sample_bool <- peaklist[,'all_sample_valid']==0


  removed_peaks <- peaklist[(valid_blank_bool | invalid_sample_bool), ,drop=F]

  grp_ids_rm <- unlist(removed_peaks[,'grpid'])

  peak_ids_rm <- unlist(xset@groupidx[grp_ids_rm])

  xset@peaks <- xset@peaks[-peak_ids_rm,]

  if( rm_peak_out){
    return(list(xset, removed_peaks))
  }else{
    return(xset)
  }


}

get_full_peak_width <- function(peaklist, xsa){
  ###########################################
  # Get full peak width
  ###########################################
  # Args:
  #   peaklist: the peak list generated from either XCMS or CAMERA.
  #              Use the CAMERA peak list for this piplein
  #   xsa: The CAMERA annotation object
  #
  # Returns:
  #   An updated peaklist with the full retention window ranges (and full mz ranges)
  #
  # See also:
  #   full_minmax, getpeaks, ldply (from the plyr library)

  print("Getting full peak widths")
  # Get 'peaks' (xcms features) from the xcmsSet object stored
  # in the camera annotation object

  print("Get 'individual' peaks from camera-xcms object")

  #Modification by M Jones to allow for retrieving full chromatographic peak widths of the excluded peaks before removal from peaklist.
  if(attributes(xsa)$class[1] != "xcmsSet"){
    peaks = data.frame(xsa@xcmsSet@peaks)
    obj = xsa@xcmsSet
  }else{
    peaks = data.frame(xsa@peaks)
    obj = xsa
  }


  rt.min = groupval(obj, method = "medret", value = "rtmin", intensity = "into")
  rt.min = apply(rt.min, 1, min, na.rm = TRUE)

  rt.max = groupval(obj, method = "medret", value = "rtmax", intensity = "into")
  rt.max = apply(rt.max, 1, min, na.rm = TRUE)

  mz.min = groupval(obj, method = "medret", value = "mzmin", intensity = "into")
  mz.min = apply(mz.min, 1, min, na.rm = TRUE)

  mz.max = groupval(obj, method = "medret", value = "mzmax", intensity = "into")
  mz.max = apply(mz.max, 1, min, na.rm = TRUE)

  peaklist_full = cbind(peaklist, "mzmin_full" = mz.min, "mzmax_full" = mz.max, "rtmin_full" = rt.min, "rtmax_full" = rt.max)

  # #peaks <- data.frame(xsa@xcmsSet@peaks)
  # peaks$peakID <- seq(1,nrow(peaks))
  #
  # print("Get associated peak for each group")
  # # Get a list of the 'peaks' assoicated with each group
  # prnt <- lapply(obj@groupidx, function(x){ peaks[x,]})
  # #prnt <- lapply(xsa@xcmsSet@groupidx, function(x){ peaks[x,]})
  #
  # print("update the peak list")
  # peaklist_xcms_full <- full_minmax(grouplist = prnt, peaks = peaks, groupdf = peaklist)

}

full_minmax <- function(grouplist,peaks,groupdf){
  # Use the the plyr package to go through a list() of peaks and output a
  # dataframe of actual rtmin, rtmax, mzmin and mzmax
  full_info <- ldply(grouplist,getpeaks,peaks=peaks)

  # Update users original dataframe
  groupdf$mzmin_full <- full_info$mzmin
  groupdf$mzmax_full <- full_info$mzmax
  groupdf$rtmin_full <- full_info$rtmin
  groupdf$rtmax_full <- full_info$rtmax

  # return the updated dataframe
  groupdf

}

getpeaks <- function(row,peaks) {
  # In case no peaks found. Should be worried if this happens though
  if (is.null(nrow(row)) || nrow(row)==0 ) {
    return(c(mzmin = NA, mzmax = NA, rtmin = NA, rtmax = NA))
  }

  # get the associated peaks for the group
  fullpeak <- peaks[row$peakID,]

  # output a vector of the actual min max values
  c(mzmin = min(fullpeak[,"mzmin"]), mzmax = max(fullpeak[,"mzmax"]),
    rtmin = min(fullpeak[,"rtmin"]), rtmax = max(fullpeak[,"rtmax"]))
}

