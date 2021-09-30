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



#' @title Using a purityA object, link MS/MS data to XCMS features
#' @description
#'
#' ## General
#'
#' Assign fragmentation spectra (MS/MS) stored within a purityA class object to grouped features within an XCMS xset object.
#'
#' XCMS calculates individual chromatographic peaks for each mzML file (saved in xset@@peaks), these are then grouped together
#' (using xcms.group). Ideally the mzML files that contain the MS/MS spectra also contain sufficient MS1 scans for XCMS to detect
#' MS1 chromatographic features. If this is the case, to determine if a MS2 spectra is to be linked to an XCMS grouped feature,
#' the associated acquisition time of the MS/MS event has to be within the retention time window defined for the individual peaks
#' associated for each file. The precursor m/z value also has to be within the user ppm tolerance to XCMS feature.
#'
#' See below for representation of the linking (the &ast; ------ &ast; represent a many-to-many relationship) e.g. 1 or more MS/MS events can be
#' linked to 1 or more individual feature and an individual XCMS feature can be linked to 1 or more grouped XCMS features
#'
#' * \[grouped XCMS feature - across files\] &ast; ------ &ast;  \[individual XCMS feature - per file\] &ast; ------ &ast;  \[MS/MS spectra\]
#'
#' Alternatively, if the "useGroup" argument is set to TRUE, the full width of the grouped peak (determined as the minimum rtmin
#' and maximum rtmax of the all associated individual peaks) will be used. This option should be used if the mzML file with
#' MS/MS has very limited MS1 data and so individual chromatographic peaks might not be detected with the mzML files containing the
#' MS/MS data. However, it should be noted this may lead to potential inaccurate linking.
#'
#' * \[grouped XCMS peaks\] &ast; ------ &ast; \[MS/MS spectra\]
#'
#'
#' ## Example LC-MS/MS processing workflow
#'
#' The purityA object can be used for further processing including linking the fragmentation spectra to XCMS features, averaging fragmentation, database creation and spectral matching (from the created database). See below for an example workflow
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> xcms.adjustRtime -> xcms.groupChromPeaks -> (xdata)
#'  * Fragmentation processing
#'    + (xdata, pa) -> **frag4feature** -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' ## Additional notes
#'
#' * If using only a single file, then grouping still needs to be performed within XCMS before frag4feature can be used.
#' * Fragmentation spectra below a certain precursor ion purity can be be removed (see plim argument).
#' * A SQLite database can be created directly here but the functionality has been deprecated and the createDatabase function should now be used
#' * Can experience some problems when using XCMS version < 3 and obiwarp retention time correction.
#'
#'
#' @aliases frag4feature
#' @param pa object; purityA object
#' @param xcmsObj object; XCMSnExp, xcmsSet or xsAnnotate object derived from the same files as those used to create the purityA object
#' @param ppm numeric; ppm tolerance between precursor mz and XCMS feature mz
#' @param plim numeric; minimum purity of precursor to be included
#' @param intense boolean; If TRUE the most intense precursor will be used. If FALSE the precursor closest to the center of the isolation window will be used
#' @param useGroup boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If retention time correction has been used in XCMS set this to TRUE
#' @param createDb boolean; if yes, generate a database of MS2 spectra
#' @param outDir string; path where (optionally generated) database file should be saved
#' @param grpPeaklist dataframe; Can use any peak dataframe to add to databse. Still needs to be derived from the "obj" object though
#' @param dbName character; name to assign to (optionally exported) database.
#' @param use_group boolean; (Deprecated, to be removed - replaced with useGroup argument for style consistency)
#' @param out_dir character; (Deprecated, to be removed - use createDatabase function) Path where database will be created
#' @param create_db boolean; (Deprecated, to be removed - use createDatabase function) SQLite database will be created of the results
#' @param grp_peaklist dataframe; (Deprecated, to be removed - use createDatabase function) Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @param db_name character; (Deprecated, to be removed - use createDatabase function) If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @return Returns a purityA object (pa) with the following slots populated:
#'
#' * pa@@grped_df: A dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor details is recorded here
#' * pa@@grped_ms2: A list of fragmentation spectra associated with each grouped XCMS feature is recorded here
#' * pa@@f4f_link_type: The linking method is recorded here (e.g. individual peaks or grouped - "useGroup=TRUE")
#'
#'
#' @examples
#' #read in MS data
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' ms_data = readMSData(msmsPths, mode = 'onDisk', msLevel. = 1)
#'
#' ## For xcms version 3.X
#'
#' #find peaks in each file
#' cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3, 30))
#' xcmsObj <- xcms::findChromPeaks(ms_data, param = cwp)

#'
#' #optionally adjust retention time
#' xcmsObj <- adjustRtime(xcmsObj, param = ObiwarpParam(binSize = 0.6))
#'
#' #group features across samples
#' sg = rep(1, length(xcmsObj$sampleNames))
#' pdp <- PeakDensityParam(sampleGroups = sg, minFraction = 0, bw = 30)
#' xcmsObj <- groupChromPeaks(xcmsObj, param = pdp)
#'
#' ## For xcms versions < 3.X
#' xcmsObj <- xcms::xcmsSet(msmsPths)
#' xcmsObj <- xcms::group(xcmsObj)
#' xcmsObj <- xcms::retcor(xcmsObj)
#' xcmsObj <- xcms::group(xcmsObj)
#'
#' ## generate purityA object and run frag4feature
#' pa  <- purityA(msmsPths)
#' pa <- frag4feature(pa, xcmsObj)
#' @export
setMethod(f="frag4feature", signature="purityA",
          definition = function(pa, xcmsObj, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE, useGroup=FALSE, createDb=FALSE,
                                outDir='.', dbName=NA, grpPeaklist=NA, use_group = NA, out_dir = NA, create_db = NA,
                                grp_peaklist = NA, db_name = NA){


  if(!is.na(use_group)){
    useGroup <- use_group
  }

  if(!is.na(out_dir)){
    outDir <- out_dir
  }

  if(!is.na(create_db)){
    createDb <- create_db
  }

  if(!is.na(grp_peaklist)){
    grpPeaklist <- grp_peaklist
  }

  if(!is.na(db_name)){
    dbName <- db_name
  }

  if('XCMSnExp' == class(xcmsObj)[1]){
    XCMSnExp_bool = TRUE
  }else if('xcmsSet' == class(xcmsObj)[1]){
    XCMSnExp_bool = FALSE
  }else if('xsAnnotate' == class(xcmsObj)){
    XCMSnExp_bool = FALSE
    xcmsObj = xcmsObj@xcmsSet
  }else{
    stop('obj is not of class XCMSnExp, xcmsSet or xsAnnotate')
  }

  # Makes sure the same files are being used
  if (!useGroup){
    pa@f4f_link_type = 'individual'
    for(i in 1:length(pa@fileList)){
      f_nms =
      if(XCMSnExp_bool){
        f_nms = basename(xcmsObj@processingData@files[i])
      }else{
        f_nms = basename(xcmsObj@filepaths[i])
      }

      if(!basename(pa@fileList[i])==f_nms){
        print("XCMSnExp/xset and pa file paths do not match")
        return(NULL)
      }
    }
  }else{
    pa@f4f_link_type = 'group'
  }

  # Get the purity data frame and the xcms peaks data frame
  puritydf <- pa@puritydf
  puritydf$fileid <- as.numeric(as.character(puritydf$fileid))

  if(XCMSnExp_bool){
    allpeaks <- data.frame(xcms::chromPeaks(xcmsObj))
    allpeaks$filename = xcmsObj$sampleName[allpeaks$sample]
  }else{
    allpeaks <- data.frame(xcmsObj@peaks)
    allpeaks$filename <- xcms::sampnames(xcmsObj)[allpeaks$sample]
    #allpeaks <- plyr::ddply(allpeaks, ~ sample, getname, xcmsObj=xcmsObj)
  }

  allpeaks$cid <- seq(1, nrow(allpeaks))

  if(convert2RawRT){

    conv_check = FALSE

    if(XCMSnExp_bool){
      if(hasAdjustedRtime(xcmsObj)){
        conv_check = TRUE
      }
    }else{
      if(any(unlist(lapply(xcmsObj@.processHistory, function(mesg){ "Retention time correction" %in% mesg@type })))){
        conv_check = TRUE
      }
    }

    if(conv_check){
      allpeaks$rtminCorrected <- allpeaks$rtmin
      allpeaks$rtmaxCorrected <- allpeaks$rtmax
      allpeaks <- plyr::ddply(allpeaks, ~ sample, convert2Raw, xcmsObj=xcmsObj, XCMSnExp_bool=XCMSnExp_bool)
    }else{
      message('convert2RawRT == TRUE but retention time alignment not applied to xcmsObj. Using raw retention times for features')
      allpeaks$rtminCorrected <- NA
      allpeaks$rtmaxCorrected <- NA
    }

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

    if(XCMSnExp_bool){
      fullpeakw <- data.frame(get_full_peak_width(xcms::featureDefinitions(xcmsObj), xcmsObj = xcmsObj))
    }else{
      fullpeakw <- data.frame(get_full_peak_width(xcmsObj@groups, xcmsObj = xcmsObj))
    }

    fullpeakw$grpid <- seq(1, nrow(fullpeakw))

    matched <- plyr::ddply(puritydf, ~ pid, fsub2, allpeaks=fullpeakw, intense=intense, ppm=ppm,
                           fullp=TRUE, use_grped=TRUE)

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
    #
    if(XCMSnExp_bool){
      grpedp <- plyr::llply(xcms::featureDefinitions(xcmsObj)$peakidx,grpByXCMS, matched=matched)
    }else{
      grpedp <- plyr::llply(xcmsObj@groupidx, grpByXCMS, matched=matched)
    }

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

  if (createDb){
    if(is.null(pa@filter_frag_params$allfrag)){
      pa@filter_frag_params$allfrag = FALSE
    }
    pa@db_path <- createDatabase(pa, xcmsObj = xcmsObj, xsa=NULL, outDir=outDir,
                                 grpPeaklist=grpPeaklist, dbName=dbName)
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
    rtmin_col <- "rtmin_full"
    rtmax_col <- "rtmax_full"
  }else{
    rtmin_col <- "rtmin"
    rtmax_col <- "rtmax"
  }

  mtchRT <- allpeaks[prt>=allpeaks[,rtmin_col] &
                     prt<=allpeaks[,rtmax_col] &
                     !is.na(allpeaks[,rtmin_col]) &
                     !is.na(allpeaks[,rtmax_col]),]

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

getname <- function(x, xcmsObj){
 x$filename <- basename(xcmsObj@filepaths[x$sample])
 return(x)
}

grpByXCMS <- function(x, matched){
  matched[matched$cid %in% x,]
}

#convert2Raw <- function(x, xset){
#  sid <- unique(x$sample)
#  # for each file get list of peaks
#  x$rtmin <- xset@rt$raw[[sid]][match(x$rtmin, xset@rt$corrected[[sid]])]
#  x$rtmax <- xset@rt$raw[[sid]][match(x$rtmax, xset@rt$corrected[[sid]])]
#  return(x)
#
#}


convert2Raw <- function(all_peaks, xcmsObj, XCMSnExp_bool){
  ## all_peaks = dataframe of chrompeaks
  ## xcmsObj = object of class XCMSnExp, xcmsSet or xsAnnotaiton.
  ## XCMSnExp_bol = boolean, where 1 means xcmsObj class == XCMSnExp, 0 means obj class == xcmsSet
  sid <- unique(all_peaks$sample)
  # for each file get list of peaks
  if(XCMSnExp_bool==1 && (class(xcmsObj) == 'XCMSnExp')){
      all_peaks$rtmin <- rtime(xcmsObj, adjusted=FALSE, bySample=T)[[sid]][match(all_peaks$rtmin, rtime(xcmsObj, adjusted = T, bySample = T)[[sid]])]
      all_peaks$rtmax <- rtime(xcmsObj, adjusted=FALSE, bySample=T)[[sid]][match(all_peaks$rtmax, rtime(xcmsObj, adjusted = T, bySample = T)[[sid]])]
  }else if(XCMSnExp_bool==0 && (class(xcmsObj) == 'xcmsSet')){
      all_peaks$rtmin <- xcmsObj@rt$raw[[sid]][match(all_peaks$rtmin, xcmsObj@rt$corrected[[sid]])]
      all_peaks$rtmax <- xcmsObj@rt$raw[[sid]][match(all_peaks$rtmax, xcmsObj@rt$corrected[[sid]])]
  }
  return(all_peaks)
}


# get_full_peak_widths <- function(peaklist, obj, XCMSnExp_bool=TRUE){
#   ###########################################
#   # Get full peak width
#   ###########################################
#   # Args:
#   #   peaklist: the peak list generated from either XCMS or CAMERA.
#   #              Use the CAMERA peak list for this piplein
#   #   xsa: The CAMERA annotation object
#   #
#   # Returns:
#   #   An updated peaklist with the full retention window ranges (and full mz ranges)
#   #
#   # See also:
#   #   full_minmax, getpeaks, ldply (from the plyr library)
#
#   message("Getting full peak widths")
#   # Get 'peaks' (xcms features) from the XCMSnSet object
#   # in the camera annotation object
#
#   message("Get 'individual' peaks from camera-xcms object")
#
#   if(attributes(obj)$class[1] == 'xsAnnotate'){
#    obj = xsa@xcmsSet
#    XCMSnExp_bool = FALSE
#   }
#
#   if(XCMSnExp_bool && (class(obj) == 'XCMSnExp')){
#     rt.min = xcms::featureValues(obj, method = "medret", value = "rtmin", intensity = "into")
#     rt.max = xcms::featureValues(obj, method = "medret", value = "rtmax", intensity = "into")
#     mz.min = xcms::featureValues(obj, method = "medret", value = "mzmin", intensity = "into")
#     mz.max = xcms::featureValues(obj, method = "medret", value = "mzmax", intensity = "into")
#   }else if (XCMSnExp_bool==FALSE && (class(obj) == 'xcmsSet')){
#     rt.min = xcms::groupval(obj, method = "medret", value = "rtmin", intensity = "into")
#     rt.max = xcms::groupval(obj, method = "medret", value = "rtmax", intensity = "into")
#     mz.min = xcms::groupval(obj, method = "medret", value = "mzmin", intensity = "into")
#     mz.max = xcms::groupval(obj, method = "medret", value = "mzmax", intensity = "into")
#   }
#
#   rt.min = apply(rt.min, 1, min, na.rm = TRUE)
#   rt.max = apply(rt.max, 1, max, na.rm = TRUE)
#   mz.min = apply(mz.min, 1, min, na.rm = TRUE)
#   mz.max = apply(mz.max, 1, max, na.rm = TRUE)
#
#   peaklist_full = cbind(peaklist, "mzmin_full" = mz.min, "mzmax_full" = mz.max, "rtmin_full" = rt.min, "rtmax_full" = rt.max)
#   return(peaklist_full)
#
# }

# This function retrieve a xset like object
getxcmsSetObject <- function(xcmsObj) {
    # XCMS 1.x
    if (class(xcmsObj) == "xcmsSet")
        return (xcmsObj)
    # XCMS 3.x
    if (class(xcmsObj) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xcmsObj, 'xcmsSet'))
        if (!is.null(xcmsObj@phenoData$sample_group))
            sampclass(xset) <- xcmsObj@phenoData$sample_group
        else
            sampclass(xset) <- "."
        return (xset)
    }
}
