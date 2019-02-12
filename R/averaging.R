group1d <- function(v, thr){

  v <- sort(v)
  a <- ''
  cl <-  1
  c <- 1

  for(i in 1:(length(v)-1)){
    a[i] <- (1e6*(v[i+1]-v[i]))/v[i]

    if(as.numeric(a[i])<thr){
      cl[i+1] <- c
    }else{
      c <- c+1
      cl[i+1] <- c
    }

  }
  return(cl)
}

performHc <- function(mzs, ppm){

  if(length(mzs)==1){return(1)}

  mdist = dist(mzs)
  averageMzPair = as.dist(outer(mzs, mzs, "+")/2)
  relativeErrors = averageMzPair * 0.000001
  m_massTolerance = mdist / relativeErrors
  clh <- fastcluster::hclust(m_massTolerance)
  # cut a desired level
  cl <- stats::cutree(clh, h = ppm)
  return(cl)
}


clustering <- function(mz, ppm=1.5, clustType="hc", cores=2){
  #print("Performing clustering/grouping")
  # Sort out the cores to use

  # perform clustering
  if(clustType=="hc"){
    #print("hc clustering")
    if (cores>1) {
      clust<-parallel::makeCluster(cores)
      doSNOW::registerDoSNOW(clust)
      parallelBool = TRUE

    } else {
      parallelBool = FALSE

    }

    mzsplit <- split(mz, ceiling(seq_along(mz)/5000))

    # Perform averaging on multiple single runs (1 file) or
    # multi-core (averaging scans in 1 file)
    clustlist  <- plyr::llply(mzsplit,
                              .parallel = parallelBool,
                              performHc, # FUNCTION
                              ppm=ppm)


    for (i in 1:length(clustlist)){
      if(i>1){
        clustlist[[i]] <- clustlist[[i]]+plus
      }
      plus <- max(clustlist[[i]])
    }

    cl <- unlist(clustlist)

    #     # stop any open connections
    #     if (cores>1) {
    #       parallel::stopCluster(clust)
    #       closeAllConnections()
    #     }

    return(cl)

  }else if (clustType=="simple"){
    # order by mz
    #print("group1d")
    return(group1d(mz, ppm))
  }


}

stde <- function(x) sd(x)/sqrt(length(x))
rsde <- function(x) (sd(x)/mean(x))*100


averageCluster <- function(x, av="median", minnum=1,
                           missingV="ignore", totalScans, normTIC, sumI=FALSE){
  # Filter out any that do not match the minimum
  if(nrow(x)<minnum){
    return(NULL)
  }


  # Get rsd (if normalise flagged then use the normalised intensity)
  if(normTIC){
    rsdRes <- (sd(x$inorm)/mean(x$inorm))*100
  }else{
    rsdRes <- rsde(x$i)

  }

  # Calc average first of the mz. We don't want to mess around with
  # missing values for this stage
  if(av=="median"){
    mz <- median(x$mz)
  }else{
    mz <- mean(x$mz)
  }

  # Zero any missing values for intensity and SNR (this is the thermo approach,
  # but might not be best approach)
  if(missingV=="zero"){
    toadd <- totalScans-nrow(x)
    if(toadd>0){
      for(i in 1:toadd){
        x <- rbind(x, c(mz,0,0,unique(x$cl)))
      }

    }
  }

  # Finally calculate the averaged SNR and intensity
  if(av=="median"){
    i <- median(x$i)
    snr <- median(x$snr)
    inorm <- median(x$inorm)
    purity <- median(x$inPurity)
    ra <- median(x$ra)
  }else{
    i <- mean(x$i)
    snr <- mean(x$snr)
    inorm <- mean(x$inorm)
    purity <- mean(x$inPurity)
    ra <- mean(x$ra)
  }

  if (sumI){
    i <- sum(x$i)
  }



  return(c("mz" = mz, "i" = i, "snr" = snr, "rsd" = rsdRes,
           "inorm" = inorm, 'count' =length(x$i), 'total'=totalScans, 'inPurity'=purity, 'ra'=ra))

}

snrFilter <- function(x, snthr, snMeth){
  # If snMethod is "precalc" then the SNR does not need to
  # be calculated required
  if (snMeth=="median"){
    x$snr <- x$i/median(x$i)
  }else if(snMeth=="mean"){
    x$snr <- x$i/mean(x$i)
  }

  x <- x[x$snr>snthr, ]

  return(x)

}
