###############################################
# Spline interpolation for purity assessment
###############################################
# Only used for testing purposes.
# The spline interpolation is both time conusming and unreliable. This function should only be used
# for comparisons and testing purposes. The default linear interpolation should be used otherwise
splinePurity <- function(row, roi_scns, minOff, maxOff,  ppm, mostIntense, scanids, plotP, plotdir){

  if(mostIntense){
    mztarget <- row$iMz
  }else{
    mztarget <- row$aMz
  }

  # Calculate ppm tolerance to look either side of the original target value
  maxMZ <- mztarget + (mztarget*0.000001)*ppm
  minMZ <- mztarget - (mztarget*0.000001)*ppm

  if(is.na(mztarget)){
    return(NA)
  }

  intensity <- plyr::ldply(roi_scns, function(roi_scn){

    roi_scn <- data.frame(roi_scn)

    sub <- roi_scn[(roi_scn[,1]>=minMZ) & (roi_scn[,1]<=maxMZ),]

    if(nrow(sub)>1){
      # Use the mz value nearest to the
      closeMtch <- sub[which.min(abs(sub[,1] - mztarget)),]
      return(unlist(closeMtch[2]))
    }else if(nrow(sub)==0){
      return(0)
    }else{
      return(sub[,2])
    }
  })

  totalI <- plyr::ldply(roi_scns, function(roi_scn){
    roi_scn <- data.frame(roi_scn)
    sub <- roi_scn[(roi_scn[,1]>=minOff) & (roi_scn[,1]<=maxOff),]
    if(nrow(sub)>1){
      return(sum(sub[,2]))
    }else if(nrow(sub)==0){
      return(0)
    }else if(nrow(sub)==1){
      return(sub[,2])
    }else{
      #print("check")
    }
  })

  idx <- c(scanids$pre, scanids$post)

  df <- cbind(idx, intensity)

  df <- cbind(df, totalI)
  colnames(df) <- c('idx',"i", 'tic')
  df$purity <- df$i/df$tic

  df[is.na(df$purity),] <- 0

  if(nrow(df)==1){
    return(df$purity)
  }


  f <- stats::splinefun(df$idx, df$purity)

  purity <- f(row$seqNum)

  if(is.nan(purity) || (is.na(purity))){
    purity <- NA
  }else if(purity>1){
    purity <- 1
  }else if (purity<0){
    purity <- 0
  }


  # save aplot of the interpolation
  if(plotP){

    name <- file.path(plotdir, paste(row$file, "_", row$id, "_interpolate_plot_", mztarget,".png", sep=""))
    plotPurity(df, row$seqNum, name, f)
  }

  return(c("inPurity"=purity))

}
