#' @include purityPD-class.R purityPD-constructor.R
NULL

#' @title Using purityPD object, calculates to average mz, intensity
#' and signal-to-noise of multiple scans from multiple MS datafiles
#' (mzML or .csv)
#'
#' @description
#' Uses a purityPD object with references to multiple MS files. For each file:
#' Averages multiple scans together,
#' see averageSpectraSingle for more information
#'
#' @param Object object = purityPD object
#' @inheritParams averageSpectraSingle
#
#' @aliases averageSpectra
#' @return  purityPD object with averaged spectra
#' @examples
#'
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityPD(fileList=inDF, cores=1, mzML=TRUE)
#' ppDIMS <- averageSpectra(ppDIMS)
#' @seealso \code{\link{averageSpectraSingle}}
#' @export
setMethod(f="averageSpectra", signature="purityPD", definition =
                  function(Object, rtscn = "all", scanRange=NA, timeRange = NA,
                           clustType="hc", ppm=1.5, snthr = 3, av="median",
                           missingV="zero", minfrac=0.6667, normTIC=FALSE,
                           snMeth="median") {
  nfiles <- nrow(Object@fileList)

  if(Object@mzML==TRUE){
    msfrOpt <- FALSE
  }else{
    msfrOpt <- TRUE
  }

  # Check if multicore
  if (Object@cores>1){
    operator <- foreach::'%dopar%'
    # The cores are split across the files. If enough cores available they
    # are also split across the clustering algorithm as well

    if (Object@cores>nfiles){
      cores4cl = floor(Object@cores/nfiles)
      cores4files <- nfiles
    } else {
      cores4cl = 1
      cores4files <- Object@cores
    }

    cl<-parallel::makeCluster(cores4files, type = "SOCK")
    doSNOW::registerDoSNOW(cl)

  }else{
    cores4files = 1
    cores4cl = 1
    operator <- foreach::'%do%'
  }

  # Perform averaging using "averageSpectraSingle" function
  # on multiple single runs (1 file) multi-core (averaging scans in 1 file)
  Object@avPeaks$orig <- operator(foreach::foreach(i=1:nfiles,
                                            .packages = c("Rcpp", "mzR")),
                                            averageSpectraSingle(filePth = as.character(Object@fileList$filepth[i]),
                                                      rtscn = rtscn,
                                                      scanRange= scanRange,
                                                      timeRange = timeRange,
                                                      clustType=clustType,
                                                      snthr = snthr,
                                                      ppm=ppm,
                                                      av = av,
                                                      missingV = missingV,
                                                      minfrac=minfrac,
                                                      snMeth = snMeth,
                                                      cores = cores4cl,
                                                      MSFileReader = msfrOpt,
                                                      normTIC = normTIC) )

  if(Object@cores>1){
    parallel::stopCluster(cl)
  }

  names(Object@avPeaks$orig) <- Object@fileList[,'name']
  Object@avPeaks$processed <- Object@avPeaks$orig

  Object@avParam$missingV  <- missingV
  Object@avParam$ppm  <- ppm
  Object@avParam$clustType <- clustType
  Object@avParam$snthr <- snthr
  Object@avParam$av <- av
  Object@avParam$minfrac <- minfrac
  Object@avParam$snMeth <- snMeth
  Object@avParam$rtscn <- rtscn
  Object@avParam$scanRange <- scanRange
  Object@avParam$timeRange <- timeRange

  return(Object)

})

#' @title Calculates to average mz, intensity and signal-to-noise of multiple
#'  scans from 1 MS datafile (mzML or .csv)
#'
#' @description
#' Averages multiple scans of mass spectrometry data together.
#' Each scan consisting of a minimum of intensity and mz values.
#'
#' Works for either mzML or a .csv file consisting of mz, i, scanid,
#'  (optional: noise, backgroun, snr)
#'
#' Signal-to-noise (SNR) can be calculated a number of ways. Default is to
#' calculate the SN for every scan as the
#' "Intensity of peak / the median intensity of the scan".
#'
#' Alternatively if using a .CSV file a precalculated snr can be on of the
#' columns and this can be used.
#'
#' The function works for LC-MS or DI-MS datasets.
#'
#' @param filePth character = Path of the file to be processed
#' @param rtscn character = Whether it is scans or retention time to be filtered. Use "all" if all scans to be used. ['rt', 'scns', 'all']
#' @param scanRange vector = Scan range (if rtscn='scns') e.g. c(40, 69)
#' @param timeRange vector = Time range (if rtscn='rt) e.g. c(10.3, 400.8) (only if using mzML file)
#' @param clustType character = Type of clustering used either Hierarchical or just simple 1dgrouping ['hc', 'simple'], default 'hc'
#' @param ppm numeric = the ppm error to cluster mz together default 1.5
#' @param snthr numeric = Signal to noise ratio threshold, default 0
#' @param av character = What type of averaging to do between peaks
#' @param missingV character = What to do with missing values (zero or ignore)
#' @param minfrac numeric = Min fraction of scans with a grouped peak to be an accepted averaged peak
#' @param cores numeric = Number of cores used to perform Hierarchical clustering WARNING: memory intensive, default 2
#' @param MSFileReader boolean = For thermo files a the MSFileReader API can extract peaklist. This can consist of an .csv file with
#'  the following columns c('mz', 'i', 'scanid', 'snr')
#' @param snMeth character = Type of snMethod to use
#
#' @return  dataframe of the median mz, intensity, signal-to-noise ratio.
#' @examples
#' mzmlPth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
#' avP <- averageSpectraSingle(mzmlPth)
#' @export
averageSpectraSingle <- function(filePth,
                                 rtscn = "all",
                                 scanRange=NA,
                                 timeRange = NA,
                                 clustType="hc",
                                 ppm=1.5,
                                 snthr = 3,
                                 cores=1,
                                 av="median",
                                 missingV="ignore",
                                 minfrac=0.6667,
                                 snMeth = "median",
                                 MSFileReader=FALSE,
                                 normTIC = FALSE){

  # Get the peaks from each scan
  # (filtering out any above the signal to noise thres)
  if(MSFileReader){
    # CSV file created from Thermo MSFileReader
    peaklist <- msfrProcess(filePth, scanRange, snthr, snMeth)
  }else{
    # mzML file to be read in by mzR
    peaklist <- mzMLProcess(filePth, rtscn, scanRange, timeRange, snthr, snMeth)
  }

  # get normalised TIC intensity
  peaklist <- ddply( peaklist, .(scanid), function(x){
    x$inorm = x$i/sum(x$i)
    return(x)
    })

  if(nrow(peaklist)<1){
    message("No peaks to process (perhaps SNR filter too high)")
    return(data.frame())
  }

  # calculate the min number of scans required for each peak
  totalScans <- length(unique(peaklist$scanid))
  minnum <- ceiling(totalScans*minfrac)

  # Order before clustering
  peaklist <- peaklist[order(peaklist$mz, decreasing = FALSE),]

  # Perform clustering
  cl <- clustering(peaklist$mz, ppm=ppm, clustType=clustType, cores = cores)

  # add groups onto the peaklist
  peaklist <- cbind(peaklist, cl)

  # re-index the matrix
  rownames(peaklist) <- seq(1, nrow(peaklist))

  if (cores>1){
    clust <- parallel::makeCluster(cores, type = "SOCK")
    doSNOW::registerDoSNOW(clust)
    parallelBool <- TRUE
  }else{
    parallelBool <- FALSE
  }

  averages <- plyr::ddply(peaklist, ~ cl, .parallel = parallelBool,
                    averageCluster, # FUNCTION
                    av=av,
                    minnum=minnum,
                    missingV=missingV,
                    totalScans=minnum,
                    normTIC = normTIC)

  if(cores>1){
    parallel::stopCluster(cl)
  }

  # get rid of the first column (the cluster id)
  if(nrow(averages)>1){
    final <-  averages[order(averages$mz, decreasing = FALSE),]
    # get rid of the first column (the cluster id)
    final$cl <- seq(1, nrow(averages))
    colnames(final)[1] <- "peakID"
  }else{
    message("no peaks! possibly the parameters used are too harsh for the data")
    return(data.frame())
  }

  return(final)

}

mzMLProcess <- function(mzmlPth, rtscn, scanRange, timeRange, snthr, snMeth){
  # Read in mzml file from mzR

  #print("reading in mzML")
  mr <- mzR::openMSfile(mzmlPth)

  # Get the peaks
  scanPeaks <- mzR::peaks(mr)
  h <- mzR::header(mr)

  # get the ms1 scans
  hms1 <- h[h$msLevel==1,]
  scans <- hms1$seqNum

  # header takes up alot of memory
  rm(h)

  # Convert RT to scan time (only available when using mzML)
  if(rtscn=="rt"){
    # Get the scan range (roi) for the rt range supplied
    rtimes <- hms1$retentionTime
    difs <- abs(rtimes-timeRange[1])
    dife <- abs(rtimes-timeRange[2])
    startTime <- rtimes[difs==min(difs)]
    strtIdx <- match(startTime, rtimes)
    endTime <- rtimes[dife==min(dife)]
    endIdx <- match(endTime, rtimes)
    scansROI <- scans[strtIdx:endIdx]
  }else if(rtscn =="all"){
    scansROI <- seq(1, length(scanPeaks))
  }else {
    #print("scans used")
    scansROI <- scans[as.numeric(scanRange[1]):as.numeric(scanRange[2])]
  }

  names(scanPeaks) <- seq(1, length(scanPeaks))

  # Get the scans for the roi
  scanList <- scanPeaks[scansROI]

  peaklist  <- plyr::ldply(scanList)

  colnames(peaklist)[1] <- "scanid"
  colnames(peaklist)[2] <- "mz"
  colnames(peaklist)[3] <- "i"

  # remove any peaks that are have zero intensity (for waters)
  peaklist <- peaklist[peaklist$i>0,]

  peaklist <- plyr::ddply(peaklist, ~ scanid,
                          snrFilter, # FUNCTION
                          snthr=snthr,
                          snMeth=snMeth)

  # get only scans of interest if selected to do so
  if (!is.na(scanRange)){
    peaklist <- peaklist[peaklist$scanid %in% scansROI]
  }

  return(peaklist)

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

msfrProcess <- function(filePth, scanRange, snthr, snMeth){
  MSFileReaderOut <- read.csv(filePth)

  keep <- c('mz', 'i', 'scanid', 'snr')
  MSFileReaderOut <- MSFileReaderOut[,(names(MSFileReaderOut ) %in% keep)]

  # First filter out on scan by scan basis those peaks below a certain threshold
  cat("snmeth:", snMeth)

  peaklist <- plyr::ddply(MSFileReaderOut, ~ scanid,
                          snrFilter, # FUNCTION
                          snthr=snthr,
                          snMeth=snMeth)

  # get only scans of interest if selected to do so
  if (!is.na(scanRange)){
    sr <- seq(scanRange[1], scanRange[2])
    peaklist <- peaklist[peaklist$scanid %in% sr]
  }

  return(peaklist)
}


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
                           missingV="ignore", totalScans, normTIC){
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
  }else{
    i <- mean(x$i)
    snr <- mean(x$snr)
    inorm <- mean(x$inorm)
  }



  return(c("mz" = mz, "i" = i, "snr" = snr, "rsd" = rsdRes, "inorm" = inorm))

}
