#' @title Using a purityA object, link MS/MS data to XCMS features
#'
#' @description
#'
#' **General**
#'
#' Assign fragmentation spectra (MS/MS) stored within a purityA class object to grouped features within an XCMS xset object.
#'
#' XCMS calculates individual chromatographic peaks for each mzML file (saved in xset@@peaks), these are then grouped together
#' (using xcms.group). Ideally the mzML files that contain the MS/MS spectra also contain sufficient MS1 scans for XCMS to detect
#' MS1 chromatographic features. If this is the case, to determine if a MS2 spectra is to be linked to an XCMS grouped feature,
#' the associated acquisition time of the MS/MS event has to be within the retention time window defined for the individual peaks
#' associated for each file. The precursor m/z value also has to be within the user ppm tolerance to XCMS feature.
#'
#' See below for representation of the linking (the \*------\* represent a many-to-many relationship) e.g. 1 or more MS/MS events can be
#' linked to 1 or more individual feature and an individual XCMS feature can be linked to 1 or more grouped XCMS features
#'
#' * \[grouped XCMS feature - across files\] \*------\*  \[individual XCMS feature - per file\] \*------\*  \[MS/MS spectra\]
#'
#' Alternatively, if the "useGroup" argument is set to TRUE, the full width of the grouped peak (determined as the minimum rtmin
#' and maximum rtmax of the all associated individual peaks) will be used. This option should be used if the mzML file with
#' MS/MS has very limited MS1 data and so individual chromatographic peaks might not be detected with the mzML files containing the
#' MS/MS data. However, it should be noted this may lead to potential inaccurate linking.
#'
#' * \[grouped XCMS peaks\] \*------\* \[MS/MS spectra\]
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#' The purityA object can be used for further processing including linking the fragmentation spectra to XCMS features, averaging fragmentation, database creation and spectral matching (from the created database). See below for an example workflow
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
#'  * Fragmentation processing
#'    + (xset, pa) -> **frag4feature** -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' **Additional notes**
#'
#' * If using only a single file, then grouping still needs to be performed within XCMS before frag4feature can be used.
#' * Fragmentation spectra below a certain precursor ion purity can be be removed (see plim argument).
#' * A SQLite database can be created directly here but the functionality has been deprecated and the createDatabase function should now be used
#' * Can experience some problems when using XCMS version < 3 and obiwarp retention time correction.
#'
#'
#' @aliases frag4feature
#'
#' @param pa object; purityA object
#' @param xset object; xcmsSet object derived from the same files as those used to create the purityA object
#' @param ppm numeric; ppm tolerance between precursor mz and XCMS feature mz
#' @param plim numeric; minimum purity of precursor to be included
#' @param intense boolean; If TRUE the most intense precursor will be used. If FALSE the precursor closest to the center of the isolation window will be used
#' @param useGroup boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If retention time correction has been used in XCMS set this to TRUE
#' @param create_db boolean; (Deprecated, to be removed - use createDatabase function) SQLite database will be created of the results
#' @param db_name character; (Deprecated, to be removed - use createDatabase function) If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @param out_dir character; (Deprecated, to be removed - use createDatabase function) Path where database will be created
#' @param grp_peaklist dataframe; (Deprecated, to be removed - use createDatabase function) Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @param use_group boolean; (Deprecated, to be removed - replaced with useGroup argument for style consistency)
#' @return Returns a purityA object (pa) with the following slots populated:
#'
#' * pa@@grped_df: A dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor details is recorded here
#' * pa@@grped_ms2: A list of fragmentation spectra associated with each grouped XCMS feature is recorded here
#' * pa@@f4f_link_type: The linking method is recorded here (e.g. individual peaks or grouped - "useGroup=TRUE")
#'
#'
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- frag4feature(pa, xset)
#' @md
#' @export
setMethod(f="frag4feature", signature="purityA",
          definition = function(pa, xset, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE, useGroup=FALSE, create_db=FALSE,
                                out_dir='.', db_name=NA, grp_peaklist=NA, use_group=NA){
  if(!is.na(use_group)){
    useGroup <- use_group
  }
  #Be sure we have a xset object
  xset <- getxcmsSetObject(xset)

  # Makes sure the same files are being used
  if (!useGroup){
    pa@f4f_link_type = 'individual'
    for(i in 1:length(pa@fileList)){
      if(!basename(pa@fileList[i])==basename(xset@filepaths[i])){
        print("xset and pa file paths do not match")
        return(NULL)
      }
    }
  }else{
    pa@f4f_link_type = 'group'
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



  if(useGroup){
    fullpeakw <- data.frame(get_full_peak_width(xset@groups, xset))
    fullpeakw$grpid <- seq(1, nrow(fullpeakw))

    matched <- plyr::ddply(puritydf, ~ pid, fsub2, allpeaks=fullpeakw, intense=intense, ppm=ppm, fullp=TRUE, use_grped=TRUE)

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

  if (useGroup){
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

  if(intense && !use_grped){
    mz1 <- pro$iMz
  }else{
    if (is.na(pro$aMz)){
      mz1 <- pro$precursorMZ
    }else{
      mz1 <- pro$aMz
    }

  }


  if(is.na(mz1) | is.null(mz1)){
    return(NULL)
  }

  prt <- pro$precursorRT
  if (is.na(prt)){
    prt <- pro$retentionTime
  }

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

    mtchRow$inPurity <- pro$inPurity
    mtchRow$pid <- pro$pid
    mtchRow$precurMtchID <- pro$seqNum
    mtchRow$precurMtchScan <- pro$precursorScanNum
    mtchRow$precurMtchRT <- pro$precursorRT
    mtchRow$precurMtchMZ <- mz1
    mtchRow$precurMtchPPM <- ppmerror
    mtchRow$retentionTime <- pro$retentionTime
    mtchRow$fileid <- pro$fileid

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

# This function retrieve a xset like object
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (!is.null(xset@phenoData$sample_group))
            sampclass(xset) <- xset@phenoData$sample_group
        else
            sampclass(xset) <- "."
        return (xset)
    }
}
