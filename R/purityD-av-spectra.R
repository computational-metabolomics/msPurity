#' @include purityD-class.R purityD-constructor.R
NULL


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


#' @rdname msPurity-deprecated
#' @section averageSpectra:
#' If \code{cores} is omitted, the call uses \code{Object@@cores}
#' without changing that stored setting.
#'
#' Uses a purityD object with references to multiple MS files. For each file:
#' Averages multiple scans together,
#' see averageSpectraSingle for more information
#'
#' \strong{Value:} purityD object with averaged spectra
#'
#' @param Object object; purityD object
#' @param cores numeric; Maximum number of workers for this call. If omitted,
#'   uses \code{Object@@cores} without changing that stored setting.
#'
#' @aliases averageSpectra
#' @export
setMethod(f="averageSpectra", signature="purityD", definition =
                  function(Object, rtscn = "all", scanRange=NA, timeRange = NA,
                           clustType="hc", ppm=1.5, snthr = 3, av="median",
                           missingV="zero", minfrac=0.6667, normTIC=FALSE,
                           snMeth="median", mzRback="pwiz", cores=1, blockSize=5000L) {
  .msp_deprecate_dims("averageSpectra")
  if (missing(cores)) cores <- max(1L, Object@cores)
  .clustering_controls(cores, blockSize)
  Object <- .pd_begin(Object)
  nfiles <- nrow(Object@fileList)

  if(Object@mzML==TRUE){
    msfrOpt <- FALSE
  }else{
    msfrOpt <- TRUE
  }

  # Check if multicore
  if (cores>1){
    operator <- foreach::'%dopar%'
    # The cores are split across the files. If enough cores available they
    # are also split across the clustering algorithm as well

    if (cores>nfiles){
      cores4cl = floor(cores/nfiles)
      cores4files <- nfiles
    } else {
      cores4cl = 1
      cores4files <- cores
    }

    cl<-parallel::makeCluster(cores4files, type = "SOCK")
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doSNOW::registerDoSNOW(cl)

  }else{
    cores4files = 1
    cores4cl = 1
    operator <- foreach::'%do%'
  }

  # Perform averaging using "averageSpectraSingle" function
  # on multiple single runs (1 file) multi-core (averaging scans in 1 file)
  Object@avPeaks$orig <- operator(foreach::foreach(i=1:nfiles,
                                            .packages = c("Rcpp", "Spectra")),
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
                                                      blockSize = blockSize,
                                                      csvFile = msfrOpt,
                                                      normTIC = normTIC,
                                                      mzRback=Object@mzRback) )


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
  Object@avParam$cores <- cores
  Object@avParam$blockSize <- blockSize

  return(.pd_end(Object))

})

#' @rdname msPurity-deprecated
#' @section averageSpectraSingle:
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
#' Alternatively if using a .csv file as input (and assigning the csvFile parameter to TRUE), a precalculated SNR can be one of the
#' columns. The precalculated SNR can then be chosen by using the option 'precalc' for the parameter snMethod
#'
#' The function will work for both LC-MS or DI-MS datasets.
#'
#' \strong{Value:} dataframe of the median mz, intensity, signal-to-noise ratio.
#'
#'
#' @param filePth character; Path of the file to be processed
#' @param rtscn character; Whether it is scans or retention time to be filtered. Use "all" if all scans to be used. \['rt', 'scns', 'all'\]
#' @param scanRange vector; Scan range (if rtscn='scns') e.g. c(40, 69)
#' @param timeRange vector; Time range (if rtscn='rt) e.g. c(10.3, 400.8) (only if using mzML file)
#' @param clustType character; Type of clustering used either Hierarchical or just simple 1D grouping \['hc', 'simple'\]
#' @param ppm numeric; The ppm error to cluster mz together
#' @param snthr numeric; Signal to noise ratio threshold
#' @param av character; What type of averaging to do between peaks
#' @param missingV character; What to do with missing values (zero or ignore)
#' @param minfrac numeric; Min fraction of scans with a grouped peak to be an accepted averaged peak
#' @param csvFile boolean; A csv file can be used as input. Useful for thermo files where the MSFileReader API can extract peaklist. This can consist of an .csv file with
#'  the following columns c('mz', 'i', 'scanid', 'snr')
#' @param snMeth character; Type of snMethod to use \['mean', 'median', 'precalc'\]. Precalc only applicable when using the csvFile parameter as TRUE
#' @param normTIC boolean; If TRUE then RSD calculation will use the normalised intensity (intensity divided by TIC) if FALSE will use standard intensity
#' @param mzRback character; deprecated and ignored. Raw data is read through Spectra, which uses the pwiz reader
#' @param MSFileReader boolean; Deprecapted. Use csvFile parameter
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
                                 csvFile=FALSE,
                                 normTIC = FALSE,
                                 mzRback='pwiz',
                                 MSFileReader=FALSE,
                                 blockSize=5000L){
  .msp_deprecate_dims("averageSpectraSingle")
  .clustering_controls(cores, blockSize)

  if (MSFileReader){
    message('MSFileReader option deprecated, please use csvFile parameter for further development')
    csvFile = MSFileReader
  }

  # Get the peaks from each scan
  # (filtering out any above the signal to noise thres)
  if(csvFile){
    # CSV file created from Thermo csvFile
    peaklist <- msfrProcess(filePth, scanRange, snthr, snMeth)
  }else{
    # mzML file to be read in by mzR
    peaklist <- mzMLProcess(filePth, rtscn, scanRange,
                            timeRange, snthr, snMeth, mzRback)
  }

  # get normalised TIC intensity
  peaklist <- plyr::ddply( peaklist, ~ scanid, function(x){
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
  cl <- clustering(peaklist$mz, ppm=ppm, clustType=clustType, cores = cores,
    blockSize = blockSize)

  # add groups onto the peaklist
  peaklist <- cbind(peaklist, cl)

  # re-index the matrix
  rownames(peaklist) <- seq(1, nrow(peaklist))

  if (cores>1){
    clust <- parallel::makeCluster(cores, type = "SOCK")
    on.exit({
      parallel::stopCluster(clust)
      foreach::registerDoSEQ()
    }, add = TRUE)
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

mzMLProcess <- function(mzmlPth, rtscn, scanRange, timeRange, snthr, snMeth, backend){
  # Read in mzml file from mzR

  #print("reading in mzML")
  sp <- .msp_read(mzmlPth)

  # Get the peaks
  scanPeaks <- .msp_peaks(sp)
  h <- .msp_header(sp)

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



msfrProcess <- function(filePth, scanRange, snthr, snMeth){
  csvFileOut <- .msp_msfr_frame(.msp_read_msfr(filePth))

  keep <- c('mz', 'i', 'scanid', 'snr')
  csvFileOut <- csvFileOut[,(names(csvFileOut ) %in% keep)]

  # First filter out on scan by scan basis those peaks below a certain threshold
  cat("snmeth:", snMeth)

  peaklist <- plyr::ddply(csvFileOut, ~ scanid,
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

#' @rdname msPurity-deprecated
#' @section groupPeaks:
#' If \code{cores} is omitted, the call uses \code{Object@@cores}
#' without changing that stored setting.
#'
#' Uses a purityD object to group all the peaklists in the 'avPeaks$processing' slot
#'
#' \strong{Value:} data.frame of peaklists grouped together by mz
#'
#' @param Object object = purityD object
#' @param ppm numeric = The ppm tolerance to group peaklists
#' @param clustType = if 'hc' the hierarchical clustering, if 'simple' the mz values will just be grouped using a simple 1D method
#' @param cores numeric; Maximum number of workers for this call. If omitted,
#'   uses \code{Object@@cores} without changing that stored setting.
#'
#' @aliases groupPeaks
#' @export
setMethod(f="groupPeaks", signature="purityD", definition =
            function(Object, ppm=3, sampleOnly=FALSE, clustType='hc',
                     cores=1, blockSize=5000L) {
              .msp_deprecate_dims("groupPeaks")
              if (missing(cores)) cores <- max(1L, Object@cores)
              .clustering_controls(cores, blockSize)
              Object <- .pd_begin(Object)
              if (sampleOnly){
                idx = Object@sampleIdx
              }else{
                idx = seq(1, nrow(Object@fileList))
              }
              # only show a limited number of columns otherwise the data frame gets too big
              shrt <- sapply(Object@avPeaks$processed[idx], simplify = FALSE,USE.NAMES = TRUE, function(x){
                x <- x[ , which(names(x) %in% c("peakID","mz",'i', 'snr', 'inorm', 'medianPurity'))]
              })

              Object@groupedPeaks <- groupPeaksEx(shrt, cores = cores, clustType = clustType,
                ppm = ppm, blockSize = blockSize)

              return(.pd_end(Object))


})

#' @rdname msPurity-deprecated
#' @section groupPeaksEx:
#' Group a list of dataframes by their m/z values
#'
#' \strong{Value:} data.frame of peaklists grouped together by mz
#'
#' @details
#' With hierarchical clustering, matching peaks stay together. The call stops
#' with an error if more than \code{blockSize} nearby peaks cannot form one group or be
#' separated without splitting matching peaks.
#'
#' @param peak_list list = A list (named) of dataframes consiting of a least the following columns \['peakID', 'mz'\]
#' @param ppm numeric = The ppm tolerance to group peaklists
#' @param clustType = if 'hc' the hierarchical clustering, if 'simple' the mz values will just be grouped using a simple 1D method
#' @param cores numeric; Maximum number of workers, as a positive integer.
#'   More workers can increase memory use.
#' @param blockSize numeric; Maximum input peak rows per hierarchical-clustering
#'   calculation, as a positive integer (default 5,000). Safe blocks can be smaller;
#'   larger groups fitting one tolerance window need no distance matrix. Used only
#'   for hierarchical clustering; distance-matrix memory grows quadratically
#'   with this limit.
#'
#' @aliases groupPeaksEx
#' @export
groupPeaksEx <- function(peak_list, cores = 1, clustType = 'hc',  ppm = 2,
    blockSize = 5000L){
  .msp_deprecate_dims("groupPeaksEx")
  .clustering_controls(cores, blockSize)
  comb <- plyr::ldply(peak_list)

  # Need to be ordered! important for when clustering (especially if when the peaks have
  # to be split for hc
  comb <- comb[order(comb$mz),]

  mz <- comb$mz

  # Cluster the peaks togther
  comb$cl <- clustering(mz, clustType = clustType, cores = cores, ppm = ppm,
    blockSize = blockSize)

  # Create a dataframe with each file is a column
  sampnms <- unique(comb$.id)

  total <-data.frame()

  for(i in 1:length(sampnms)){
    snmi = sampnms[i]

    sub <- comb[comb$.id == snmi, ]

    # if multiple mz values in the group then average
    sub <- plyr::ddply(sub, ~ cl, medGroup)

    addnew <- sub[ , -which(names(sub) %in% c("cl", ".id"))]

    for(j in 1:ncol(addnew)){
      colnames(addnew)[j] <- paste(colnames(addnew)[j], snmi, sep="_")

    }
    addnew$cl <- sub$cl

    if (nrow(total)==0){
      total <- addnew
    }else{
      total <- merge(x = total, y =addnew, by = "cl", all = TRUE)
    }
  }


  mzall <- total[ , which(names(total) %in% grep(".*mz.*", colnames(total), value=TRUE))]

  mzmedian <- apply(mzall, 1, median, na.rm=TRUE)
  total <- total[-1]

  total <- total[,order(colnames(total))]
  total <- cbind(mzmedian, total)

  return(total)



}

medGroup <- function(x){
  if(nrow(x)>1){
    medx <- apply(x[,-which(names(x) %in% c(".id"))], 2, median, na.rm=TRUE)
    x <- data.frame(.id=unique(x$.id), t(medx))
  }

  return(x)
}





