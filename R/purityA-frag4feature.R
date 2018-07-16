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
#'
#' @export
setMethod(f="frag4feature", signature="purityA",
          definition = function(pa, xset, ppm=5, plim=0, intense=TRUE, convert2RawRT=TRUE, create_db=FALSE,
                                out_dir='.', db_name=NA, grp_peaklist=NA, use_group=FALSE){

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
  }

  para=FALSE

  if(use_group){
    fullpeakw <- data.frame(get_full_peak_width(xset@groups, xset))
    fullpeakw$grpid <- seq(1, nrow(fullpeakw))

    if(convert2RawRT){
      fullpeakw$rtminCorrected_full <- fullpeakw$rtmin_full
      fullpeakw$rtmaxCorrected_full <- fullpeakw$rtmax_full
      fullpeakw <- ddply(fullpeakw, ~ sample, convert2Raw, xset=xset)
    }
    matched <- ddply(puritydf, ~ pid, fsub2, allpeaks=fullpeakw, intense=intense, ppm=ppm, fullp=TRUE, use_grped=TRUE)

  }else{
    # Map xcms features to the data frame (takes a while)
    matched <- plyr::ddply(puritydf, ~ fileid, .parallel = para, fsub1,
                           allpeaks=allpeaks,
                           ppm = ppm,
                           intense = intense)
  }


  if(pa@cores>1){
      parallel::stopCluster(cl)
  }

  #shrt <- puritydf[,c('fileid', 'seqNum', 'inPurity','pid')]

  if (use_group){
    grpedp <- matched
  }else{
    #Group by the xcms groups
    grpedp <- plyr::llply(xset@groupidx, grpByXCMS, matched=matched)
    names(grpedp) <- seq(1, length(grpedp))
    grpedp <- plyr::ldply(grpedp, .id = TRUE)
    colnames(grpedp)[1] <- "grpid"
  }

  # Add some extra info for filtering purposes
  #grpm <- merge(grpedp, shrt, by = c('pid', 'fileid', 'seqNum'))
  grpm <- grpedp

  # Make sure order is by grpid
  grpm <- grpm[order(grpm$grpid),]

  # Filter out any precursor below purity threshold
  if (!is.na(plim) && plim>0){
    grpm <- grpm[grpm$inPurity>plim,]
  }

  # add to the slots
  pa@grped_df <- grpm
  pa@grped_ms2 <- getMS2scans(grpm, pa@fileList, mzRback = pa@mzRback)


  if (create_db){
    pa@db_path <- create_database(pa=pa, xset=xset, out_dir=out_dir,
                                  db_name=db_name, grp_peaklist=grp_peaklist)
  }


  return(pa)




})


fsub1  <- function(prod, allpeaks, intense, ppm){
  # go through all the MS/MS files from each file
  allpeakfile <- allpeaks[allpeaks$filename==unique(prod$filename),]

  grpdFile <- plyr::ddply(prod, ~ seqNum,
                          fsub2, # FUNCTION
                          allpeaks = allpeakfile,
                          intense = intense,
                          ppm = ppm)
}

fsub2  <- function(pro, allpeaks, intense, ppm, fullp=FALSE, use_grped=FALSE){
  # check for each MS/MS scan if there is an associated feature
   #found in that region for that file
  if(intense){
    mz1 <- pro$iMz
  }else{
    if ('aMz' %in% colnames(pro)){
      mz1 <- pro$aMz
    }else{
      mz1 <- pro$precursorMZ
    }

  }

  if(is.na(mz1) | is.null(mz1)){
    return(NULL)
  }

  prt <- pro$precursorRT
  if (fullp){
    mtchRT <- allpeaks[prt>=allpeaks$rtmin_full & prt<=allpeaks$rtmax_full, ]
  }else{
    mtchRT <- allpeaks[prt>=allpeaks$rtmin & prt<=allpeaks$rtmax, ]
  }


  if(nrow(mtchRT)==0){
    return(NULL)
  }
  if (use_grped){
    # can only use fullp when using the grouped peaklist
    mtchMZ <- plyr::ddply(mtchRT, ~ grpid, mzmatching, mz1=mz1, ppm=ppm, pro=pro)
  }else{
    mtchMZ <- plyr::ddply(mtchRT, ~ cid, mzmatching, mz1=mz1, ppm=ppm, pro=pro)
  }
  return(mtchMZ)

}



check_ppm <- function(mz1, mz2){ return(abs(1e6*(mz1-mz2)/mz2)) }

getMS2scans  <- function(grpm, filepths, mzRback){
  # Get all MS2 scans

  scans <- getscans(filepths, mzRback)

  if(length(filepths)==1){
    scans = list(scans)
  }

  grpm$fid <- seq(1, nrow(grpm))

  ms2l <- plyr::dlply(grpm, ~ grpid, getScanLoop, scans=scans)

  return(ms2l)
}


mzmatching <- function(mtchRow, mz1=mz1, ppm=ppm, pro=pro){
  if ('mzmed' %in% colnames(mtchRow)){
    mz2 <- mtchRow$mzmed
  }else{
    mz2 <- mtchRow$mz
  }

  ppmerror <- check_ppm(mz1, mz2)

  if(ppmerror<ppm){
    mtchRow$pid <- pro$pid
    mtchRow$precurMtchID <- pro$seqNum
    mtchRow$precurMtchScan <- pro$precursorScanNum
    mtchRow$precurMtchRT <- pro$precursorRT
    mtchRow$precurMtchMZ <- mz1
    mtchRow$precurMtchPPM <- ppmerror
    mtchRow$inPurity <- pro$inPurity
    mtchRow$seqNum <- pro$seqNum
    return(mtchRow)
  }else{
    return(NULL)
  }
}

getScanLoop <- function(peaks, scans){
  grpl <-  list()

  if ('sample' %in% colnames(peaks)){
    idx_nm ='sample'
  }else{
    idx_nm = 'fileid'
  }
  for(i in 1:nrow(peaks)){
    x <- peaks[i,]
    idx <- x[,idx_nm]
    grpl[[i]] <- scans[[idx]][[x$precurMtchID]]

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


