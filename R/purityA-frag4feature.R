#' @title Assign precursor purity scored fragmentation spectra to XCMS features
#'
#' @description
#'
#' Assign fragmentation spectra (MS/MS) scored via msPurity package to features
#' from an XCMS set object.
#'
#' Allows the user to filter out spectra below a certain threshold for purity.
#'
#' @aliases frag4feature
#'
#' @param pa = purityA object
#' @param xset xcms object = XCMS object derived from the same files as the puritydf
#' @param ppm numeric = ppm tolerance between precursor mz and feature mz
#' @param plim numeric = min purity of precursor to be included
#' @param intense boolean = If the most intense precursor or the centered precursor is used
#' @param convert2RawRT boolean = If retention time correction has been used in XCMS set this to TRUE
#' @return a dataframe of the purity score of the ms/ms spectra
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
#'
#' @export
setMethod(f="frag4feature", signature="purityA",
          definition = function(pa, xset, ppm = 5, plim = 0, intense=TRUE, convert2RawRT=TRUE){

  # Makes sure the same files are being used
  for(i in 1:length(pa@fileList)){
    if(!basename(pa@fileList[i])==basename(xset@filepaths[i])){
      print("xset and pa file paths do not match")
      return(NULL)
    }
  }

  # Get the purity data frame and the xcms peaks data frame
  puritydf <- pa@puritydf
  puritydf$fileid <- as.numeric(puritydf$fileid)
  allpeaks <- data.frame(xset@peaks)
  allpeaks$cid <- seq(1, nrow(allpeaks))
  allpeaks <- plyr::ddply(allpeaks, ~ sample, getname, xset=xset)

  if(convert2RawRT){
    allpeaks$rtminCorrected <- allpeaks$rtmin
    allpeaks$rtmaxCorrected <- allpeaks$rtmax
    allpeaks <- ddply(allpeaks, ~ sample, convert2Raw, xset=xset)
  }


  # Check if is going to be multi-core
  if(pa@cores>1){
    cl <- parallel::makeCluster(pa@cores)
    doSNOW::registerDoSNOW(cl)
    para = TRUE
  }else{
    para = FALSE

    # perform multicore

  }

  para=FALSE

  # Map xcms features to the data frame (takes a while)
  matched <- plyr::ddply(puritydf, ~ fileid, .parallel = para, fsub1,
                                        allpeaks=allpeaks,
                                        ppm = ppm,
                                        intense = intense)

  if(pa@cores>1){
      parallel::stopCluster(cl)
  }

  #Group by the xcms groups
  grpedp <- plyr::llply(xset@groupidx, grpByXCMS, matched=matched)
  names(grpedp) <- seq(1, length(grpedp))
  grpedp <- plyr::ldply(grpedp, .id = TRUE)
  colnames(grpedp)[1] <- "grpid"

  # Add some extra info for filtering purposes
  shrt <- puritydf[,c('fileid', 'seqNum', 'inPurity','pid')]
  grpm <- merge(grpedp, shrt, by = c("seqNum", "fileid"))

  # remove the first two (columns duplicates)
  grpm <- grpm[,!(names(grpm) %in% c("seqNum", "fileid"))]

  # Make sure order is by grpid
  grpm <- grpm[order(grpm$grpid),]

  # Filter out any precursor below purity threshold
  grpm <- grpm[grpm$inPurity>plim,]

  # add to the slots
  pa@grped_df <- grpm
  pa@grped_ms2 <- getMS2scans(grpm, pa@fileList, mzRback = pa@mzRback)

  return(pa)

})

fsub1  <- function(prod, allpeaks, intense, ppm){
  # go through all the MS/MS files from the each file
  allpeakfile <- allpeaks[allpeaks$filename==unique(prod$filename),]

  grpdFile <- plyr::ddply(prod, ~ seqNum,
                          fsub2, # FUNCTION
                          allpeaks = allpeakfile,
                          intense = intense,
                          ppm = ppm)
}

fsub2  <- function(pro, allpeaks, intense, ppm){
  # check for each MS/MS scan if there is an associated feature
   #found in that region for that file
  if(intense){
    mz1 <- pro$iMz
  }else{
    mz1 <- pro$aMz
  }

  if(is.na(mz1) | is.null(mz1)){
    return(NULL)
  }

  prt <- pro$precursorRT

  mtchRT <- allpeaks[prt>=allpeaks$rtmin & prt<=allpeaks$rtmax, ]

  if(nrow(mtchRT)==0){
    return(NULL)
  }

  mtchMZ <- plyr::ddply(mtchRT, ~ cid, mzmatching, mz1=mz1, ppm=ppm, pro=pro)

  return(mtchMZ)

}



check_ppm <- function(mz1, mz2){ return(abs(1e6*(mz1-mz2)/mz2)) }

getMS2scans  <- function(grpm, filepths, mzRback){
  # Get all MS2 scans

  scans <- getscans(filepths, mzRback)

  grpm$fid <- seq(1, nrow(grpm))

  ms2l <- plyr::dlply(grpm, ~ grpid, getScanLoop, scans=scans)

  return(ms2l)
}


mzmatching <- function(mtchRow, mz1=mz1, ppm=ppm, pro=pro){
  mz2 <- mtchRow$mz
  ppmerror <- check_ppm(mz1, mz2)

  if(ppmerror<ppm){
    mtchRow$precurMtchID <- pro$seqNum
    mtchRow$precurMtchRT <- pro$precursorRT
    mtchRow$precurMtchMZ <- mz1
    mtchRow$precurMtchPPM <- ppmerror
    return(mtchRow)
  }else{
    return(NULL)
  }
}

getScanLoop <- function(peaks, scans){
  grpl <-  list()
  for(i in 1:nrow(peaks)){
    x <- peaks[i,]
    grpl[[i]] <- scans[[x$sample]][[x$precurMtchID]]

  }
  return(grpl)
}

getname <- function(x, xset){
  x$filename <- basename(xset@filepaths[x$sample])
  return(x)
}

grpByXCMS <- function(x, matched){
  matched[matched$cid %in% x,]
}

convert2Raw <- function(x, xset){
  sid <- unique(x$sample)
  # for each file get list of peaks
  x$rtmin <- xset@rt$raw[[sid]][match(x$rtmin, xset@rt$corrected[[sid]])]
  x$rtmax <- xset@rt$raw[[sid]][match(x$rtmax, xset@rt$corrected[[sid]])]
  return(x)

}

