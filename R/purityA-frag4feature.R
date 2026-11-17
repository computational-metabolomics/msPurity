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
#' XCMS calculates individual chromatographic peaks for each mzML file (retrieved using xcms::chromPeaks(xcmsObj)), these are then grouped together
#' (using xcms::groupChromPeaks). Ideally the mzML files that contain the MS/MS spectra also contain sufficient MS1 scans for XCMS to detect
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
#' With \code{includeMSn = TRUE}, MS3 and MS4 inherit feature membership through
#' validated, same-file acquisition relationships. Root MS2 associations are
#' established before \code{plim} is applied, so a rejected MS2 precursor does
#' not remove valid descendant links. Higher levels retain their own precursor
#' and acquisition metadata and have no MS1-based precursor-purity value.
#' Missing or invalid ancestry does not establish a descendant feature link.
#' Acquisition times remain raw; \code{convert2RawRT} converts feature bounds
#' using the existing xcms retention-time mapping.
#'
#' ## Example LC-MS/MS processing workflow
#'
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> (optionally) xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'    +  --- *Older versions of XCMS* --- (mzML files) -> xcms.xcmsSet -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> **frag4feature** -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
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
#' @param xcmsObj object; XcmsExperiment, XCMSnExp, xcmsSet or xsAnnotate object derived from the same files as those used to create the purityA object
#' @param ppm numeric; ppm tolerance between precursor mz and XCMS feature mz
#' @param plim numeric; minimum MS2 precursor purity to be included.
#'   This threshold does not apply to MS3/MS4 diagnostic purity.
#' @param intense boolean; If TRUE the most intense precursor will be used. If FALSE the precursor closest to the center of the isolation window will be used
#' @param useGroup boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If retention time correction has been used in XCMS set this to TRUE
#' @param includeMSn logical; also associate MS3 and MS4 spectra through their
#'   validated MS2 acquisition ancestors. The default retains MS2-only links.
#' @param parentReference How to interpret parent references. Numeric
#'   \code{precScanNum} values can refer to \code{"acquisitionNum"}, retained
#'   original file positions (\code{"scanIndex"}), or backend source keys
#'   (\code{"sourceSpectrumId"}). The last option is msPurity-specific and
#'   uses retained \code{spectrum_id_} or \code{source_spectrum_id} keys.
#'   \code{"spectrumId"} matches native \code{spectrumRef} metadata to native
#'   \code{spectrumId}; it requires these fields to have been supplied by the
#'   reader. No new metadata extraction is performed. None of these identifiers
#'   is the current position in a filtered object. All matches use the same run.
#'   \code{"auto"} tries these interpretations and accepts a consistent mapping;
#'   conflicting valid parents or inconsistent conventions raise an error.
#'   Missing references remain unresolved. No parent is inferred from peak
#'   similarity or proximity in retention time.
#' @details Purity is assessed by \code{purityA(msLevel = 2:4)} before feature
#'   association. This function reuses prepared parent relationships and keeps
#'   the recorded assessments when features are relinked; it does not calculate
#'   purity. With \code{includeMSn = TRUE}, MS3/MS4 inherit the feature links of
#'   their MS2 ancestors and retain their own purity values and status.
#'   Both this function's \code{plim} and the purity threshold in
#'   \code{filterFragSpectra()} currently apply only to MS2. For unassessed MS2
#'   roots, use \code{intense = FALSE} to associate recorded precursor m/z.
#'   MSn averaging and processed export are separate operations.
#' @param createDb boolean; if yes, generate a database of MS2 spectra
#' @param outDir string; path where (optionally generated) database file should be saved
#' @param grpPeaklist dataframe; Can use any peak dataframe to add to databse. Still needs to be derived from the "obj" object though
#' @param dbName character; name to assign to (optionally exported) database.
#' @param use_group boolean; (Deprecated, to be removed - replaced with useGroup argument for style consistency)
#' @param out_dir character; (Deprecated, to be removed - use createDatabase function) Path where database will be created
#' @param create_db boolean; (Deprecated, to be removed - use createDatabase function) SQLite database will be created of the results
#' @param grp_peaklist dataframe; (Deprecated, to be removed - use createDatabase function) Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @param db_name character; (Deprecated, to be removed - use createDatabase function) If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @param xset object; (Deprecated, to be removed - use xcmsObj) 'xcmsSet' object derived from the same files as those used to create the purityA object
#'
#' @return Returns a purityA object (pa) with the following slots populated:
#'
#' * pa@@grped_df: A dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor details is recorded here
#' * pa@@fragSpectra: The linked fragmentation spectra as a \code{Spectra} object, each scan once, identified by
#'   its \code{pid}. \code{groupedSpectra(pa)} returns them as a list per XCMS feature
#' * pa@@grped_ms2: The legacy list of fragmentation spectra for each grouped XCMS feature, filled while
#'   \code{options(msPurity.legacySlots = TRUE)} (the default)
#' * pa@@f4f_link_type: The linking method is recorded here (e.g. individual peaks or grouped - "useGroup=TRUE")
#'
#'
#' @examples
#' #====== XCMS =================================
#' ## Read in MS data
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#' #           package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #ms_data = readMSData(msmsPths, mode = 'onDisk', msLevel. = 1)
#'
#' ## Find peaks in each file
#' #cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3, 30))
#' #xcmsObj  <- xcms::findChromPeaks(ms_data, param = cwp)
#'
#' ## Optionally adjust retention time
#' #xcmsObj  <- adjustRtime(xcmsObj , param = ObiwarpParam(binSize = 0.6))
#'
#' ## Group features across samples
#' #pdp <- PeakDensityParam(sampleGroups = c(1, 1), minFraction = 0, bw = 30)
#' #xcmsObj <- groupChromPeaks(xcmsObj , param = pdp)
#'
#' ## Or if using the old XCMS functions
#' #xcmsObj <- xcms::xcmsSet(msmsPths)
#' #xcmsObj <- xcms::group(xcmsObj)
#' #xcmsObj <- xcms::retcor(xcmsObj)
#' #xcmsObj <- xcms::group(xcmsObj)
#'
#' #====== msPurity ============================
#' #pa <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xcmsObj)
#'
#' ## Run from previously generated data
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#'            package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "1_purityA_pa.rds", package="msPurity"))
#' xcmsObj <- readRDS(system.file("extdata", "tests", "xcms",
#'                                "msms_only_xcmsnexp.rds", package="msPurity"))
#' ## The saved objects record the file paths of the machine that made them
#' pa@fileList[1:2] <- msmsPths
#' xcmsObj@processingData@files <- msmsPths
#' pa <- frag4feature(pa, xcmsObj)
#'
#' @export
setMethod(f="frag4feature", signature="purityA",
          definition = function(pa, xcmsObj, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE, useGroup=FALSE, createDb=FALSE,
                                outDir='.', dbName=NA, grpPeaklist=NA, use_group = NA, out_dir = NA, create_db = NA,
                                grp_peaklist = NA, db_name = NA, xset = NA,
                                includeMSn = FALSE, parentReference = "auto") {
  pa <- .pa_update(pa)

  if (!is.logical(includeMSn) || length(includeMSn) != 1L || anyNA(includeMSn))
    stop("includeMSn must be TRUE or FALSE.")
  parentReference <- match.arg(parentReference,
    c("auto", "acquisitionNum", "scanIndex", "spectrumId", "sourceSpectrumId"))
  .pa_msn_validate_preparation(pa)
  processing_started <- Sys.time()
  processing_input <- .pa_processing_input(pa, "frag4feature")

  if(!is.na(xset)){
    message('The param xset is deprecated - please use xcmsObj instead')
    xcmsObj <- xset
  }

  if(!is.na(use_group)){
    message('The param use_group is deprecated - please use useGroup instead')
    useGroup <- use_group
  }

  if(!is.na(out_dir)){
    message('The param out_dir is deprecated - please use outDir instead')
    outDir <- out_dir
  }

  if(!is.na(create_db)){
    message('The param create_db is deprecated - please use createDb instead')
    createDb <- create_db
  }

  if(!is.na(grp_peaklist)){
    message('The param grp_peaklist is deprecated - please use grpPeaklist instead')
    grpPeaklist <- grp_peaklist
  }

  if(!is.na(db_name)){
    message('The param db_name is deprecated - please use dbName instead')
    dbName <- db_name
  }

  if (includeMSn && isTRUE(createDb))
    stop("MSn feature association cannot use the deprecated createDb route. ",
      "Associate spectra first; higher-level database export requires a ",
      "supported MSn export route.")

  if(.xcms_is_modern(xcmsObj)){
    XCMSnExp_bool = TRUE
  }else if(is(xcmsObj, 'xcmsSet')){
    XCMSnExp_bool = FALSE
  }else if(is(xcmsObj, 'xsAnnotate')){
    XCMSnExp_bool = FALSE
    xcmsObj = xcmsObj@xcmsSet
  }else{
    stop('xcmsObj is not of class XcmsExperiment, XCMSnExp, xcmsSet or xsAnnotate')
  }

  # Makes sure the same files are being used
  if (!useGroup){
    pa@f4f_link_type = 'individual'
    for(i in 1:length(pa@fileList)){
      
      if(XCMSnExp_bool){
        f_nms = basename(.xcms_files(xcmsObj)[i])
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
  puritydf <- purityTable(pa)
  puritydf$fileid <- as.numeric(as.character(puritydf$fileid))

  if(XCMSnExp_bool){
    allpeaks <- data.frame(xcms::chromPeaks(xcmsObj))
    allpeaks$filename = basename(.xcms_files(xcmsObj))[allpeaks$sample]
    #allpeaks$filename = basename(xcmsObj$sampleName)[allpeaks$sample]
  }else{
    allpeaks <- data.frame(xcmsObj@peaks)
    allpeaks$filename <- basename(xcmsObj@filepaths)[allpeaks$sample]
    #allpeaks <- plyr::ddply(allpeaks, ~ sample, getname, xcmsObj=xcmsObj)
  }

  allpeaks$cid <- seq_len(nrow(allpeaks))

  if (convert2RawRT && nrow(allpeaks)) {

    conv_check = FALSE

    if(XCMSnExp_bool){
      if(xcms::hasAdjustedRtime(xcmsObj)){
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

  # No MS2 roots or no matching features still yield a typed empty result.
  grpm <- allpeaks[FALSE, , drop = FALSE]
  empty <- list(grpid = character(), pid = integer(), sample = integer(),
    cid = integer(), inPurity = numeric(), fileid = integer(),
    seqNum = integer(), retentionTime = numeric(), precurMtchMZ = numeric())
  for (name in setdiff(names(empty), names(grpm))) grpm[[name]] <- empty[[name]]

  if (nrow(puritydf) && nrow(allpeaks)) {
    para <- pa@cores > 1
    if (para) {
      cl <- parallel::makeCluster(pa@cores)
      doSNOW::registerDoSNOW(cl)
    }
    if (useGroup) {
      definitions <- if (XCMSnExp_bool) xcms::featureDefinitions(xcmsObj) else
        xcmsObj@groups
      fullpeakw <- data.frame(get_full_peak_width(definitions,
        xcmsObj = xcmsObj))
      fullpeakw$grpid <- seq_len(nrow(fullpeakw))
      matched <- plyr::ddply(puritydf, ~ pid, fsub2, allpeaks = fullpeakw,
        intense = intense, ppm = ppm, fullp = TRUE, use_grped = TRUE)
    } else {
      matched <- plyr::ddply(puritydf, ~ fileid, .parallel = para, fsub1,
        allpeaks = allpeaks, ppm = ppm, intense = intense)
    }
    if (para) parallel::stopCluster(cl)
    if (nrow(matched)) {
      if (useGroup) {
        grpm <- matched
      } else {
        membership <- if (XCMSnExp_bool)
          xcms::featureDefinitions(xcmsObj)$peakidx else xcmsObj@groupidx
        grouped <- plyr::llply(membership, grpByXCMS, matched = matched)
        names(grouped) <- seq_along(grouped)
        grouped <- plyr::ldply(grouped, .id = TRUE)
        if (nrow(grouped)) {
          colnames(grouped)[1L] <- "grpid"
          grpm <- grouped
        }
      }
    }
  }

  # Make sure order is by grpid
  grpm <- grpm[order(grpm$grpid), , drop = FALSE]

  # Descendants use the MS2 associations before MS2 precursor-purity exclusion.
  pa <- .pa_processing_relink(pa)
  if (includeMSn) {
    preparation <- .pa_msn_preparation(pa, parentReference)
    grpm <- .pa_msn_feature_links(grpm, preparation, plim)
    fields <- intersect(c("inPurity", "msPurityStagePurity",
        "msPurityPurityStatus"), Spectra::spectraVariables(pa@spectra))
    for (field in fields)
      grpm[[field]] <- pa@spectra[[field]][match(grpm$pid, pa@spectra$pid)]
    pa@params$msnLineage <- preparation
  } else if (!is.na(plim) && plim > 0) {
    grpm <- grpm[which(grpm$inPurity > plim), , drop = FALSE]
  }

  # Record the arguments for provenance.
  prm <- if (methods::.hasSlot(pa, "params")) pa@params else list()
  prm$frag4feature <- list(ppm = ppm, plim = plim, intense = intense,
                           convert2RawRT = convert2RawRT, useGroup = useGroup,
                           xcmsObj_class = class(xcmsObj)[1],
                           includeMSn = includeMSn,
                           parentReference = parentReference)
  pa@params <- prm
  if (includeMSn) {
    pa@params$frag4feature$resolvedParentReference <-
      preparation$parentReference
  }

  # add to the slots
  pa@grped_df <- grpm
  # Each linked scan is stored once; grped_ms2 is rebuilt from the links.
  pa@fragSpectra <- if (nrow(grpm)) .pa_raw_scans(pa, unique(grpm$pid)) else
    pa@spectra[integer()]
  pa <- .pa_sync_legacy(pa, grped_ms2 = .pa_grouped_legacy(pa))

  if (createDb){
    if(is.null(pa@filter_frag_params$allfrag)){
      pa@filter_frag_params$allfrag = FALSE
    }
    pa@db_path <- .createDatabase_sqlite(pa, xcmsObj = xcmsObj, xsa=NULL, outDir=outDir,
                                         grpPeaklist=grpPeaklist, dbName=dbName)
    .msp_deprecate_sqlite("frag4feature(createDb = TRUE)",
                          paste("Call createDatabase(format = \"parquet\")",
                                "after frag4feature() instead."))
  }

        effective_parameters <- c(pa@params$frag4feature,
            list(createDb = createDb, outDir = outDir, dbName = dbName,
                grpPeaklist = grpPeaklist))
        pa <- .pa_processing_record(pa, "frag4feature", effective_parameters,
            processing_started, processing_input)
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

getname <- function(x, xcmsObj){
 x$filename <- basename(xcmsObj@filepaths[x$sample])
 return(x)
}

grpByXCMS <- function(x, matched){
  matched[matched$cid %in% x,]
}

convert2Raw <- function(all_peaks, xcmsObj, XCMSnExp_bool){
  ## all_peaks = dataframe of chrompeaks
  ## xcmsObj = object of class XCMSnExp, xcmsSet or xsAnnotaiton.
  ## XCMSnExp_bol = boolean, where 1 means xcmsObj class == XCMSnExp, 0 means obj class == xcmsSet
  sid <- unique(all_peaks$sample)
  # for each file get list of peaks
  if(XCMSnExp_bool==1 && (is(xcmsObj,  'XCMSnExp'))){
      all_peaks$rtmin <- xcms::rtime(xcmsObj, adjusted=FALSE, bySample=TRUE)[[sid]][match(all_peaks$rtmin, xcms::rtime(xcmsObj, adjusted = TRUE, bySample = TRUE)[[sid]])]
      all_peaks$rtmax <- xcms::rtime(xcmsObj, adjusted=FALSE, bySample=TRUE)[[sid]][match(all_peaks$rtmax, xcms::rtime(xcmsObj, adjusted = TRUE, bySample = TRUE)[[sid]])]
  }else if(XCMSnExp_bool==1 && (is(xcmsObj, 'XcmsExperiment'))){
      raw <- .xcms_rtime_by_sample(xcmsObj, adjusted = FALSE)[[sid]]
      adj <- .xcms_rtime_by_sample(xcmsObj, adjusted = TRUE)[[sid]]
      all_peaks$rtmin <- raw[match(all_peaks$rtmin, adj)]
      all_peaks$rtmax <- raw[match(all_peaks$rtmax, adj)]
  }else if(XCMSnExp_bool==0 && (is(xcmsObj, 'xcmsSet'))){
      all_peaks$rtmin <- xcmsObj@rt$raw[[sid]][match(all_peaks$rtmin, xcmsObj@rt$corrected[[sid]])]
      all_peaks$rtmax <- xcmsObj@rt$raw[[sid]][match(all_peaks$rtmax, xcmsObj@rt$corrected[[sid]])]
  }
  return(all_peaks)
}
