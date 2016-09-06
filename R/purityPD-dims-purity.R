#' @include purityPD-constructor.R
NULL

#' @title Using purityPD object, predict the precursor purity of processed peaks
#'
#' @description
#' Uses a purityPD object with references to multiple MS files. Predicts the
#' purity of the processed sample files
#'
#' @aliases dimsPredictPurity
#'
#' @param Object object = purityPD object
#' @inheritParams dimsPredictPuritySingle
#
#' @return  purityPD object with predicted purity of peaks
#' @examples
#'
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' ppDIMS <- purityPD(fileList=inDF, cores=1, mzML=TRUE)
#' ppDIMS <- averageSpectra(ppDIMS)
#' ppDIMS <- filterp(ppDIMS)
#' ppDIMS <- subtract(ppDIMS)
#' ppDIMS <- dimsPredictPurity(ppDIMS)
#' @return purityPD object
#' @seealso \code{\link{dimsPredictPuritySingle}}
#'
#'
#' @export
setMethod(f="dimsPredictPurity", signature="purityPD",
          definition= function(Object, ppm = 1.5, minOffset=0.5, maxOffset=0.5,
                               iwNorm=FALSE, iwNormFun=NULL, ilim=0.05) {
            requireNamespace('foreach')

            # Get the sample peaks to get the purity of
            sidx <- Object@sampleIdx

            Object@purityParam$minOffset = minOffset
            Object@purityParam$maxOffset = minOffset
            Object@purityParam$ppm = ppm
            Object@purityParam$iwNorm = iwNorm
            Object@purityParam$iwNormFun = iwNormFun
            Object@purityParam$ilim= ilim

            # Check if multicore
            if (Object@cores>1){
              operator <- foreach::'%dopar%'
              cl<-parallel::makeCluster(Object@cores, type = "SOCK")
              doSNOW::registerDoSNOW(cl)

            }else{
              operator <- foreach::'%do%'
            }

            samplePeaksAll <- operator(foreach::foreach(i=1:length(sidx), .packages = "mzR"),
                                       predictPurityExp(Object, sidx[[i]]))

            for (i in 1:length(sidx)){
              Object@avPeaks$processed[[sidx[i]]] <- samplePeaksAll[[i]]
            }

            return(Object)
          })



predictPurityExp <- function(Object, fidx){
  samplePeaks <- Object@avPeaks$processed[[fidx]]

  if(nrow(samplePeaks)==0){
    return(samplePeaks)
  }

  filepth <- as.character(Object@fileList$filepth[fidx])

  # if iwNorm is TRUE and iwNormFun is NULL
  if(is.null(Object@purityParam$iwNormFun)){
    Object@purityParam$iwNormFun <- iwNormRcosine()
  }

  purity <- dimsPredictPuritySingle(mztargets = samplePeaks$mz,
                                    filepth = filepth,
                                    minOffset = Object@purityParam$minOffset,
                                    maxOffset = Object@purityParam$maxOffset,
                                    ppm = Object@purityParam$ppm,
                                    mzML = Object@mzML,
                                    iwNorm=Object@purityParam$iwNorm,
                                    iwNormFun=Object@purityParam$iwNormFun,
                                    ilim=Object@purityParam$ilim)

  samplePeaks <- cbind(samplePeaks, purity)

  return(samplePeaks)
}



#' @title Predict the precursor purity from a DI-MS dataset
#'
#' @description
#' Given a an DI-MS dataset (either mzML or .csv file) calculate the predicted
#' purity for a vector of mz values.
#'
#' Calculated at a given offset e.g. for 0.5 +/- Da the minOffset would be 0.5
#' and the maxOffset of 0.5.
#'
#' A ppm tolerance is used to find the target mz value in each scan.
#'
#' @param mztargets vector = mz targets to get predicted purity for
#' @param filepth character = mzML file path or .csv file path
#' @param minOffset numeric = isolation window minimum offset
#' @param maxOffset numeric = isolation window maximum offset
#' @param ppm numeric = tolerance for target mz value in each scan
#' @param mzML boolean = Whether an mzML file is to be used or .csv file (TRUE == mzML)
#' @param iwNorm boolean = if TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function = A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric = All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5\% (0.05)
#' @examples
#' mzmlPth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
#' predicted <- dimsPredictPuritySingle(c(173.0806, 216.1045), filepth=mzmlPth , minOffset=0.5, maxOffset=0.5, ppm=5, mzML=TRUE)
#' @return a dataframe of the target mz values and the predicted purity score
#' @export
dimsPredictPuritySingle <- function(mztargets,
                                    filepth,
                                    minOffset=0.5,
                                    maxOffset=0.5,
                                    ppm=2.5,
                                    mzML=TRUE,
                                    iwNorm=FALSE,
                                    iwNormFun=NULL,
                                    ilim=0.05){

  # open the file and get the scans
  if(mzML==TRUE){
    # mzML files opened with mzR
    loadNamespace('mzR')
    mr <- mzR::openMSfile(filepth)
    scanPeaks <- mzR::peaks(mr)
    h <- mzR::header(mr)

    # only want the ms1 scans
    hms1 <- h[h$msLevel==1,]
    scans <- hms1$seqNum
    rm(h)
    # get peaks from each scan
    scanPeaks <- mzR::peaks(mr)

    # filter out any that are not ms1
    scanPeaks <- scanPeaks[scans]

  }else{
    # MSFileReader outputs opened with as .csv files
    MSfile <- read.csv(filepth)
    scanPeaks <- plyr::dlply(MSfile, ~ scanid, function(x){x[-1]})
  }

  # if iwNorm is TRUE and iwNormFun is NULL then a gaussian model of the
  # isolation window will be used to normalise intensity
  if(is.null(iwNormFun)){
    # Using a gaussian curve 3 SD either side
    iwNormFun <- iwNormGauss(3)
  }

  # perform the purity prediction on each target mz value
  pureList <- lapply(mztargets, dimsPredictPuritySingleMz,
                     scanPeaks=scanPeaks,
                     minOffset=minOffset,
                     maxOffset=maxOffset,
                     ppm=ppm,
                     iwNorm=iwNorm,
                     iwNormFun=iwNormFun,
                     ilim = ilim)
  puredf <- do.call(rbind.data.frame, pureList)
  colnames(puredf) <- c('medianPurity','meanPurity',
                        'sdPurity', 'cvPurity', 'sdePurity', "medianPeakNum")
  return(puredf)

}

dimsPredictPuritySingleMz <- function(mz, scanPeaks, minOffset, maxOffset, ppm,
                                      plot=FALSE, plotdirpth,
                                      iwNorm=FALSE, iwNormFun=NULL, ilim=0.05){
  # Get isolation window
  minmz <- mz-minOffset
  maxmz <- mz+maxOffset

  purityall <- ""
  pknmall <- ""

  for (i in 1:length(scanPeaks)){
    x <- scanPeaks[[i]]

    pout <- pcalc(x, minmz, maxmz, mz, ppm, iwNorm, iwNormFun, ilim)
    purityi <- pout[1]
    pknm <- pout[2]

    if(plot==TRUE){

      png(file.path(plotdirpth, paste("scan_", i, "_", mz,".png",sep="" )),
          width=10,height=10,units="in",res=1200)

      plot(sub, type="h", xlim=c(minmz, maxmz),  xlab="m/z", ylab="Intensity",
           main=paste("Isolation window surrounding m/z value", mz),
           cex.lab=2, cex.axis=2, cex.main=2, cex.sub=2, lwd=4)

      details <- paste("target I =", round(mtchi,0), "\ntotal I =", round(alli,0),
                       "\nPurity (target/total) =",round(purityi,3))
      points(mtch[1], mtch[2], type="h", col="red", lwd=5)
      text(mtch[1], mtch[2], paste(details),pos = 4, cex=2)
      dev.off()
    }

    purityall <- c(purityall, purityi)
    pknmall <- c(pknmall, pknm)

  }

  purityall <- as.numeric(purityall[-1])
  pknmmpall <- as.numeric(pknmall[-1])

  puritySum <- c(median(purityall), mean(purityall), sd(purityall),
                 covar(purityall), stderror(purityall),
                 median(pknm, na.rm = TRUE))

  return(puritySum)

}

covar <- function(x){ ( 100*sd(x)/mean(x) )} # CV (otherwise known as RSD)
stderror <- function(x){ sd(x)/sqrt(length(x))}
