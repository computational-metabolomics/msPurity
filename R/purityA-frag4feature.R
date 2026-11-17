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
#' and maximum rtmax of the representative peak of each sample, see \code{groupWindow}) will be used. This option should be used if the mzML file with
#' MS/MS has very limited MS1 data and so individual chromatographic peaks might not be detected with the mzML files containing the
#' MS/MS data. However, it should be noted this may lead to potential inaccurate linking.
#'
#' * \[grouped XCMS peaks\] &ast; ------ &ast; \[MS/MS spectra\]
#'
#' With \code{msLevel = 2:4}, MS3 and MS4 inherit feature membership through
#' validated, same-file acquisition relationships. Root MS2 associations are
#' established before \code{plim} is applied, so a rejected MS2 precursor does
#' not remove valid descendant links. Higher levels retain their own precursor
#' and acquisition metadata and their own recorded purity assessments.
#' Missing or invalid ancestry does not establish a descendant feature link.
#' Only MS2 precursors are compared with peak and feature windows; the
#' precursor time of an MS3/MS4 row is mapped to the adjusted scale of its file
#' for reference.
#' Selecting only higher levels, such as \code{msLevel = 4L}, retains their
#' required MS2/MS3 source context while returning only the selected feature
#' links.
#'
#' ## Retention time correction
#'
#' MS/MS times are raw times read from each mzML file, while XCMS peak and feature windows are on the adjusted scale once
#' retention time correction has been applied. With \code{convert2RawRT = TRUE}, the time of each MS/MS scan is mapped to the
#' adjusted scale of its file, by linear interpolation between the raw and adjusted scan times XCMS holds for that file, and is
#' compared there, in both linking modes. A time equal to an XCMS scan time maps exactly to the adjusted time XCMS gave that scan.
#'
#' Files are matched to XCMS samples by name, so the MS/MS files can be any subset of the XCMS samples in any order. Each file gets a
#' mapping status, stored with the other diagnostics in \code{pa@@params$frag4feature$rt_diagnostics}:
#'
#' * aligned: the file is in the XCMS object, which holds its adjusted times;
#' * unaligned: no retention time correction was applied, times are used unchanged;
#' * external: the mapping comes from \code{ms2RtMap};
#' * estimated: the mapping was estimated from the MS/MS precursors (\code{estimateRtMap = TRUE});
#' * missing: the features are aligned but nothing maps this file, so its raw times are compared with adjusted windows (with a warning);
#' * not_converted: \code{convert2RawRT = FALSE}.
#'
#' When the MS/MS files were not part of the retention time correction (for example DDA files processed apart from the full-scan
#' samples), supply their mapping with \code{ms2RtMap}. The most reliable route is landmark alignment in xcms: run
#' \code{findChromPeaks} on the DDA files, then
#' \code{adjustRtime(ddaObj, LamaParama(lamas = featureDefinitions(ms1Obj)[, c("mzmed", "rtmed")]))}, and pass \code{ddaObj}.
#' When XCMS fits the correction, it is best fitted on the full-scan samples only, with
#' \code{PeakGroupsParam(subset = <MS1 samples>, subsetAdjust = "average")}, as DDA files have sparse MS1 data.
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
#' * Peaks added by xcms::fillChromPeaks are not used for linking unless \code{includeFilled = TRUE}.
#' * One MS/MS scan can be linked to several peaks, so it can appear more than once for a feature when the feature holds two peaks
#'   of the same file. \code{matchRank} orders the matches of each scan.
#' * A SQLite database can be created directly here but the functionality has been deprecated and the createDatabase function should now be used
#'
#'
#' @aliases frag4feature
#' @param pa object; purityA object
#' @param xcmsObj object; XcmsExperiment, XCMSnExp, xcmsSet or xsAnnotate object derived from the same files as those used to create the purityA object
#' @param ppm numeric; ppm tolerance between precursor mz and XCMS feature mz
#' @param plim numeric; minimum MS2 precursor purity to be included.
#'   This threshold does not apply to MS3/MS4 diagnostic purity.
#' @param intense boolean; If TRUE the most intense precursor (iMz) will be used. If FALSE the precursor closest to the center of the isolation
#'        window (aMz) will be used. Missing values fall back to aMz and then to the precursor m/z of the mzML file. Applies to both linking modes.
#' @param useGroup boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If TRUE, MS/MS times are mapped to the adjusted retention time scale before being compared with XCMS windows.
#'        If FALSE, raw MS/MS times are compared with the XCMS windows as they are.
#' @param createDb boolean; if yes, generate a database of MS2 spectra
#' @param outDir string; path where (optionally generated) database file should be saved
#' @param grpPeaklist dataframe; Can use any peak dataframe to add to databse. Still needs to be derived from the "obj" object though
#' @param dbName character; name to assign to (optionally exported) database.
#' @param ms2RtMap object; mapping to the adjusted scale for MS/MS files the XCMS object does not align. Either an XCMS object holding those files
#'        with adjusted retention times, a list of functions from raw to adjusted time, or a numeric vector of offsets (adjusted = raw + offset).
#'        Lists and vectors are named by file. Files aligned in \code{xcmsObj} keep that alignment.
#' @param rtTol numeric; tolerance in seconds added to both sides of each retention time window, to absorb error left after mapping
#' @param rtSource character; MS/MS time to compare: "precursor" uses the time of the MS1 survey scan (falling back to the MS/MS scan time),
#'        "ms2" uses the time of the MS/MS scan
#' @param groupWindow character; window of a feature when \code{useGroup = TRUE}: "full" uses the minimum rtmin and maximum rtmax of the
#'        representative peak of each sample, "median" their median rtmin and rtmax
#' @param includeFilled boolean; If TRUE, peaks added by xcms::fillChromPeaks are used for linking and flagged in the is_filled column.
#'        A link through a filled peak means the MS/MS fell where the feature was expected in that file, not inside a detected peak.
#'        DDA files have sparse MS1 data, so many of their features exist only as filled peaks.
#' @param estimateRtMap boolean; If TRUE and \code{useGroup = TRUE}, files with no mapping get one estimated from their precursors. Scans matching exactly one
#'        feature are landmarks: a first pass looks 60 seconds (or rtTol if wider) either side of each feature rtmed and takes the median
#'        offset, and two more passes look 15 and 8 seconds around the current estimate. A robust loess of the offset against time is used
#'        when it fits the landmarks clearly better than a constant offset. The method and the number of landmarks are reported in the
#'        diagnostics. An estimate is weaker than a mapping from xcms, so it is off by default
#' @param use_group boolean; (Deprecated, to be removed - replaced with useGroup argument for style consistency)
#' @param out_dir character; (Deprecated, to be removed - use createDatabase function) Path where database will be created
#' @param create_db boolean; (Deprecated, to be removed - use createDatabase function) SQLite database will be created of the results
#' @param grp_peaklist dataframe; (Deprecated, to be removed - use createDatabase function) Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @param db_name character; (Deprecated, to be removed - use createDatabase function) If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @param xset object; (Deprecated, to be removed - use xcmsObj) 'xcmsSet' object derived from the same files as those used to create the purityA object
#' @param msLevel Nonempty numeric vector of output levels selected from
#'   2, 3 and 4. The default retains MS2-only links. MS3/MS4 inherit membership
#'   through their validated MS2 acquisition ancestors.
#' @param param One complete \code{\link{Frag4FeatureParam}} object, or
#'   \code{NULL} to use scalar settings. Explicit scalar association arguments
#'   supplied alongside the object must agree, including \code{use_group}.
#'   \code{ms2RtMap} is an input rather than a setting and is always supplied
#'   as an argument.
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
#'   purity. When requested by \code{msLevel}, MS3/MS4 inherit their MS2
#'   ancestors' feature links and retain their own purity values and status.
#'   This function's \code{plim} applies only to MS2. Use
#'   \code{filterFragSpectra()} for optional purity thresholds at each level.
#'   For unassessed MS2
#'   roots, use \code{intense = FALSE} to associate recorded precursor m/z.
#'   MSn averaging and processed export are separate operations.
#'
#'   Processing history records requested and effective workers, whether the
#'   retention time mapping was requested and applied, and a fingerprint of the
#'   prepared feature bounds and memberships used for association. The mapping
#'   is used by both association modes. No XCMS object or peak arrays are
#'   copied into history. The fingerprint identifies the prepared association
#'   inputs, not the complete original data or a replayable analysis.
#'   Association runs as one serial join, so effective workers is one when
#'   matching runs and zero when no matching runs because MS2 roots or
#'   chromatographic peaks are unavailable.
#'
#' @return Returns a purityA object (pa) with the following slots populated:
#'
#' * pa@@grped_df: A dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor details is recorded here.
#'   Besides the XCMS columns and the precursor details, it holds
#'   + precurMtchRTadj: the MS/MS time on the adjusted scale, as compared with the window
#'   + rtMapStatus: the mapping status of the file (see above)
#'   + precurMzSource: the precursor m/z used (iMz, aMz or precursorMZ)
#'   + is_filled: whether the peak was added by fillChromPeaks (individual mode)
#'   + rtminWin, rtmaxWin: the feature window used (grouped mode)
#'   + rtApexDist: distance from the peak or feature apex, divided by the window width
#'   + matchRank: rank of the match among the matches of the scan, by rtApexDist then ppm error
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
                                outDir='.', dbName=NA, grpPeaklist=NA, ms2RtMap=NULL, rtTol=0,
                                rtSource=c("precursor", "ms2"), groupWindow=c("full", "median"),
                                includeFilled=FALSE, estimateRtMap=FALSE, use_group = NA, out_dir = NA, create_db = NA,
                                grp_peaklist = NA, db_name = NA, xset = NA,
                                msLevel = 2L, parentReference = "auto",
                                param = NULL) {
  supplied <- c(ppm = !missing(ppm), plim = !missing(plim),
    intense = !missing(intense), convert2RawRT = !missing(convert2RawRT),
    useGroup = !missing(useGroup), parentReference = !missing(parentReference),
    rtTol = !missing(rtTol), rtSource = !missing(rtSource),
    groupWindow = !missing(groupWindow), includeFilled = !missing(includeFilled),
    estimateRtMap = !missing(estimateRtMap))
  pa <- .pa_update(pa)

  if (!is.numeric(msLevel) || !length(msLevel) || anyNA(msLevel) ||
      any(!msLevel %in% 2:4))
    stop("msLevel must be a nonempty numeric vector of levels 2, 3 or 4.")
  msLevel <- unique(as.integer(msLevel))
  parentReference <- match.arg(parentReference,
    c("auto", "acquisitionNum", "scanIndex", "spectrumId", "sourceSpectrumId"))
  if (!is.null(param)) {
    settings <- list(ppm = ppm, plim = plim, intense = intense,
      convert2RawRT = convert2RawRT, useGroup = useGroup,
      parentReference = parentReference, rtTol = rtTol, rtSource = rtSource,
      groupWindow = groupWindow, includeFilled = includeFilled,
      estimateRtMap = estimateRtMap)[supplied]
    if (!is.na(use_group)) settings$use_group <- use_group
    settings <- .pa_feature_resolve(param, settings)
    ppm <- settings$ppm
    plim <- settings$plim
    intense <- settings$intense
    convert2RawRT <- settings$convert2RawRT
    useGroup <- settings$useGroup
    parentReference <- settings$parentReference
    rtTol <- settings$rtTol
    rtSource <- settings$rtSource
    groupWindow <- settings$groupWindow
    includeFilled <- settings$includeFilled
    estimateRtMap <- settings$estimateRtMap
  }
  .pa_msn_validate_preparation(pa)
  processing_started <- Sys.time()
  processing_input <- .pa_processing_input(pa, "frag4feature")

  if(!missing(xset)){
    message('The param xset is deprecated - please use xcmsObj instead')
    xcmsObj <- xset
  }

  if(!missing(use_group)){
    message('The param use_group is deprecated - please use useGroup instead')
    useGroup <- use_group
  }

  if(!missing(out_dir)){
    message('The param out_dir is deprecated - please use outDir instead')
    outDir <- out_dir
  }

  if(!missing(create_db)){
    message('The param create_db is deprecated - please use createDb instead')
    createDb <- create_db
  }

  if(!missing(grp_peaklist)){
    message('The param grp_peaklist is deprecated - please use grpPeaklist instead')
    grpPeaklist <- grp_peaklist
  }

  if(!missing(db_name)){
    message('The param db_name is deprecated - please use dbName instead')
    dbName <- db_name
  }

  if (any(msLevel > 2L) && isTRUE(createDb))
    stop("MSn feature association cannot use the deprecated createDb route. ",
      "Associate spectra first; higher-level database export requires a ",
      "supported MSn export route.")

  rtSource <- match.arg(rtSource)
  groupWindow <- match.arg(groupWindow)
  if (!is.numeric(rtTol) || length(rtTol) != 1L || is.na(rtTol) || rtTol < 0)
    stop("rtTol must be a single non-negative number of seconds.")

  if(is(xcmsObj, 'xsAnnotate')){
    xcmsObj <- xcmsObj@xcmsSet
  }
  if(!.xcms_is_modern(xcmsObj) && !is(xcmsObj, 'xcmsSet')){
    stop('xcmsObj is not of class XcmsExperiment, XCMSnExp, xcmsSet or xsAnnotate')
  }

  pa@f4f_link_type <- if (useGroup) 'group' else 'individual'

  # Files are matched to XCMS samples by name
  pa_files <- unname(as.character(pa@fileList))
  source_files <- .xcms_files(xcmsObj)
  sample_idx <- .f4f_file_index(pa_files, source_files)
  if (!useGroup) {
    if (all(is.na(sample_idx)))
      stop("None of the purityA files are in the XCMS object (",
           paste(.f4f_basename(pa_files), collapse = ", "),
           "). Use useGroup = TRUE to link MS/MS from files that XCMS did not process.")
    if (anyNA(sample_idx))
      warning("The XCMS object does not hold ",
              paste(.f4f_basename(pa_files[is.na(sample_idx)]), collapse = ", "),
              ", so their MS/MS cannot be linked to individual peaks. ",
              "Use useGroup = TRUE to link them to features.", call. = FALSE)
  }

  # Retention time map of each file
  aligned <- .f4f_rt_pairs(xcmsObj)
  if (convert2RawRT && is.null(aligned))
    message('convert2RawRT == TRUE but retention time alignment not applied to xcmsObj. Using raw retention times for features')
  if (!convert2RawRT && !is.null(aligned))
    message('convert2RawRT == FALSE: raw MS/MS retention times are compared with adjusted XCMS retention times')
  maps <- .f4f_rt_maps(pa_files, sample_idx, if (convert2RawRT) aligned,
                       convert2RawRT, ms2RtMap)

  # No MS2 roots or no chromatographic peaks still yield a typed empty result.
  puritydf <- purityTable(pa)
  n_peaks <- nrow(if (.xcms_is_modern(xcmsObj)) xcms::chromPeaks(xcmsObj)
                  else xcmsObj@peaks)
  run <- nrow(puritydf) > 0L && n_peaks > 0L
  grpm <- .f4f_empty_links()
  scans <- data.frame(fileid = numeric())
  features <- .f4f_features(xcmsObj)
  estimated <- list()
  effective_cores <- 0L
  association_input <- "not_used"
  used_bounds <- membership <- NULL

  if (run) {
    scans <- .f4f_scans(puritydf, intense, rtSource)

    if (estimateRtMap && useGroup) {
      # In individual mode a file the XCMS object lacks has no peaks to link to
      for (i in which(maps$status == "missing")) {
        est <- .f4f_estimate_map(scans[scans$fileid == i, , drop = FALSE],
                                 features, ppm, max(60, rtTol),
                                 if (is.na(plim)) 0 else plim)
        estimated[[as.character(i)]] <- est[c("method", "landmarks")]
        if (est$method != "none") {
          maps$fun[[i]] <- est$fun
          maps$status[i] <- "estimated"
        } else {
          warning("Too few landmarks (", est$landmarks, ") to estimate a retention time map for ",
                  .f4f_basename(pa_files[i]), call. = FALSE)
        }
      }
    }

    missing_map <- maps$status == "missing" & (useGroup | !is.na(sample_idx)) &
      seq_along(pa_files) %in% scans$fileid
    if (any(missing_map))
      warning("No retention time map for ",
              paste(.f4f_basename(pa_files[missing_map]), collapse = ", "),
              ": their raw MS/MS times are compared with adjusted XCMS windows. ",
              "Supply ms2RtMap, or set estimateRtMap = TRUE.", call. = FALSE)

    scans$rtAdj <- NA_real_
    scans$rtMapStatus <- NA_character_
    for (i in unique(scans$fileid)) {
      r <- which(scans$fileid == i)
      scans$rtAdj[r] <- maps$fun[[i]](scans$rtRaw[r])
      scans$rtMapStatus[r] <- maps$status[i]
    }

    if (useGroup) {
      linked <- .f4f_link_group(scans, xcmsObj, ppm, rtTol, groupWindow)
      association_input <- "group_bounds"
    } else {
      scans$sample <- sample_idx[scans$fileid]
      linked <- .f4f_link_individual(scans, xcmsObj, ppm, rtTol, includeFilled,
                                     convert2RawRT, aligned)
      association_input <- "individual_peaks"
    }
    effective_cores <- 1L
    used_bounds <- linked$bounds
    membership <- linked$membership
    if (nrow(linked$links))
      grpm <- linked$links
  }

  mapped <- any(maps$status %in% c("aligned", "external", "estimated"))
  rt_conversion <- list(
    preparation = if (!convert2RawRT) "not_requested"
                  else if (!n_peaks) "no_chromatographic_peaks"
                  else if (mapped) "applied" else "no_adjustment",
    used_for_association = convert2RawRT && run && mapped)

  processing_input$xcms <- .pa_feature_input_reference(
    class(xcmsObj)[1L], source_files, association_input,
    used_bounds, membership)

  # Descendants use the MS2 associations before MS2 precursor-purity exclusion.
  pa <- .pa_processing_relink(pa)
  if (any(msLevel > 2L)) {
    preparation <- .pa_msn_preparation(pa, parentReference)
    grpm <- .pa_msn_feature_links(grpm, preparation, plim)
    fields <- intersect(c("inPurity", "msPurityStagePurity",
        "msPurityPurityStatus"), Spectra::spectraVariables(pa@spectra))
    for (field in fields)
      grpm[[field]] <- pa@spectra[[field]][match(grpm$pid, pa@spectra$pid)]
    pa@params$msnLineage <- preparation
    grpm <- grpm[grpm$msLevel %in% msLevel, , drop = FALSE]
    grpm <- .f4f_descendant_times(grpm, maps)
  } else if (!is.na(plim) && plim > 0) {
    grpm <- grpm[which(grpm$inPurity > plim), , drop = FALSE]
  }
  rownames(grpm) <- NULL

  roots <- if ("msLevel" %in% names(grpm))
    grpm[grpm$msLevel == 2L, , drop = FALSE] else grpm
  diagnostics <- .f4f_diagnostics(roots, scans, features, pa_files,
                                  sample_idx, maps$status, estimated)

  # Record the arguments for provenance.
  prm <- if (methods::.hasSlot(pa, "params")) pa@params else list()
  prm$frag4feature <- list(ppm = ppm, plim = plim, intense = intense,
                           convert2RawRT = convert2RawRT, useGroup = useGroup,
                           xcmsObj_class = class(xcmsObj)[1],
                           ms2RtMap = if (is.null(ms2RtMap)) NA_character_
                                      else class(ms2RtMap)[1],
                           rtTol = rtTol, rtSource = rtSource,
                           groupWindow = groupWindow,
                           includeFilled = includeFilled,
                           estimateRtMap = estimateRtMap,
                           rt_diagnostics = diagnostics,
                           cores = pa@cores, effective_cores = effective_cores,
                           rt_conversion = rt_conversion,
                           msLevel = msLevel,
                           effective_msLevel = unique(Spectra::msLevel(
                             pa@spectra)[match(grpm$pid, pa@spectra$pid)]),
                           parentReference = parentReference,
                           param = param)
  pa@params <- prm
  if (any(msLevel > 2L)) {
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


.pa_feature_input_reference <- function(
  class, files, mode, bounds, membership
) {
  reference <- list(class = class, files = files, association_input = mode,
    fingerprint = NULL)
  if (identical(mode, "not_used")) return(reference)
  schema <- "frag4feature-input-v1"
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path))
  saveRDS(list(schema = schema, association_input = mode, bounds = bounds,
      memberships = membership), path, compress = FALSE, version = 2L)
  reference$fingerprint <- list(schema = schema, serialization_version = 2L,
    algorithm = "sha256", value = .mzs_sha256(path))
  reference
}

check_ppm <- function(mz1, mz2){ return(abs(1e6*(mz1-mz2)/mz2)) }

#' grped_df with no links, with the column types of a linked result.
#'
#' @noRd
.f4f_empty_links <- function() {
  data.frame(grpid = character(), pid = integer(), sample = integer(),
             cid = integer(), inPurity = numeric(), fileid = integer(),
             seqNum = integer(), retentionTime = numeric(),
             precurMtchMZ = numeric(), precurMtchRTadj = numeric(),
             rtMapStatus = character(), precurMzSource = character(),
             rtApexDist = numeric(), matchRank = integer(),
             stringsAsFactors = FALSE)
}

#' MS3/MS4 rows inherit the link of their MS2 root. Their precursor time and
#' m/z are their own, so the adjusted time is mapped again from the precursor
#' time and the m/z source is the recorded precursor m/z.
#'
#' @noRd
.f4f_descendant_times <- function(grpm, maps) {
  child <- which(grpm$msLevel > 2L)
  if (!length(child))
    return(grpm)
  if ("precurMtchRTadj" %in% names(grpm)) {
    file <- as.integer(as.character(grpm$fileid[child]))
    rt <- as.numeric(grpm$precurMtchRT[child])
    adj <- rep(NA_real_, length(child))
    for (i in unique(file[!is.na(file)])) {
      r <- which(file == i)
      adj[r] <- maps$fun[[i]](rt[r])
    }
    grpm$precurMtchRTadj[child] <- adj
  }
  if ("precurMzSource" %in% names(grpm))
    grpm$precurMzSource[child] <- "precursorMZ"
  grpm
}

#' The MS/MS scans of the purity table, with the precursor m/z and the raw
#' time to compare.
#'
#' With `intense`, the precursor m/z is iMz, then aMz, then the precursor m/z
#' of the mzML file; otherwise aMz, then the mzML value.
#'
#' @noRd
.f4f_scans <- function(puritydf, intense, rtSource) {
  # Purity tables built from Spectra can lack columns of the mzML route.
  col <- function(name) {
    v <- puritydf[[name]]
    if (is.null(v)) rep(NA_real_, nrow(puritydf)) else v
  }
  s <- data.frame(pid = puritydf$pid,
                  fileid = as.numeric(as.character(puritydf$fileid)),
                  seqNum = puritydf$seqNum,
                  inPurity = col("inPurity"),
                  precursorScanNum = col("precursorScanNum"),
                  precursorRT = col("precursorRT"),
                  retentionTime = col("retentionTime"))
  candidates <- if (intense) c("iMz", "aMz", "precursorMZ") else c("aMz", "precursorMZ")
  s$mz1 <- NA_real_
  s$precurMzSource <- NA_character_
  for (cn in candidates) {
    v <- as.numeric(puritydf[[cn]])
    take <- is.na(s$mz1) & !is.na(v)
    s$mz1[take] <- v[take]
    s$precurMzSource[take] <- cn
  }
  s$rtRaw <- if (rtSource == "ms2") s$retentionTime
             else ifelse(is.na(s$precursorRT), s$retentionTime, s$precursorRT)
  s
}

#' Feature id, m/z, retention time and member peaks of each XCMS feature.
#'
#' @noRd
.f4f_features <- function(xcmsObj) {
  if (.xcms_is_modern(xcmsObj)) {
    fd <- xcms::featureDefinitions(xcmsObj)
    idx <- fd$peakidx
  } else {
    fd <- xcmsObj@groups
    idx <- xcmsObj@groupidx
  }
  data.frame(grpid = seq_len(nrow(fd)),
             mzmed = as.numeric(fd[, "mzmed"]),
             rtmed = as.numeric(fd[, "rtmed"]),
             npeaks = lengths(idx))[, c("grpid", "mzmed", "rtmed")]
}

#' Peaks added by fillChromPeaks.
#'
#' @noRd
.f4f_is_filled <- function(xcmsObj, n) {
  out <- logical(n)
  if (.xcms_is_modern(xcmsObj)) {
    f <- xcms::chromPeakData(xcmsObj)$is_filled
    if (!is.null(f))
      out <- as.logical(f) %in% TRUE
  } else if (length(xcmsObj@filled)) {
    out[xcmsObj@filled] <- TRUE
  }
  out
}

#' Interval join of MS/MS scans and windows.
#'
#' `windows` has the columns .lo, .hi (retention time window, adjusted scale)
#' and .mz (reference m/z), and optionally sample. Matches meet
#' check_ppm() < ppm, as before.
#'
#' @return The scan row (.si) and window row (.wi) of each match.
#' @noRd
.f4f_join <- function(scans, windows, ppm, by_sample) {
  s <- data.frame(.si = seq_len(nrow(scans)), rt = scans$rtAdj, mz = scans$mz1)
  # The bounds only preselect; the exact ppm test below decides, so they are
  # widened slightly to keep rounding from excluding a boundary match.
  tol <- ppm * 1e-6 * (1 + 1e-6)
  w <- data.frame(.wi = seq_len(nrow(windows)), .lo = windows$.lo,
                  .hi = windows$.hi, .mzlo = windows$.mz * (1 - tol),
                  .mzhi = windows$.mz * (1 + tol))
  if (by_sample) {
    s$sample <- scans$sample
    w$sample <- windows$sample
  }
  s <- s[!is.na(s$rt) & !is.na(s$mz) & (if (by_sample) !is.na(s$sample) else TRUE), , drop = FALSE]
  w <- w[!is.na(w$.lo) & !is.na(w$.hi) & !is.na(w$.mzlo), , drop = FALSE]
  by <- if (by_sample)
    dplyr::join_by(sample, between(rt, .lo, .hi), between(mz, .mzlo, .mzhi))
  else
    dplyr::join_by(between(rt, .lo, .hi), between(mz, .mzlo, .mzhi))
  m <- dplyr::inner_join(s, w, by = by, relationship = "many-to-many")
  m <- m[check_ppm(m$mz, windows$.mz[m$.wi]) < ppm, c(".si", ".wi"), drop = FALSE]
  as.data.frame(m)
}

#' The precursor details of each matched scan, in the column order of
#' grped_df.
#'
#' @noRd
.f4f_match_columns <- function(out, scans, m, ppm_ref) {
  sc <- scans[m$.si, , drop = FALSE]
  out$inPurity <- sc$inPurity
  out$pid <- sc$pid
  out$precurMtchID <- sc$seqNum
  out$precurMtchScan <- sc$precursorScanNum
  out$precurMtchRT <- sc$precursorRT
  out$precurMtchMZ <- sc$mz1
  out$precurMtchPPM <- check_ppm(sc$mz1, ppm_ref)
  out$retentionTime <- sc$retentionTime
  out$fileid <- sc$fileid
  out$seqNum <- sc$seqNum
  out$precurMtchRTadj <- sc$rtAdj
  out$rtMapStatus <- sc$rtMapStatus
  out$precurMzSource <- sc$precurMzSource
  out
}

#' Rank the matches of each scan by distance from the apex, relative to the
#' window width, and then by ppm error.
#'
#' @noRd
.f4f_rank <- function(out, apex, lo, hi) {
  width <- hi - lo
  out$rtApexDist <- ifelse(width > 0, abs(out$precurMtchRTadj - apex) / width,
                           NA_real_)
  if (!nrow(out)) {
    out$matchRank <- integer()
    return(out)
  }
  o <- order(out$pid, out$rtApexDist, out$precurMtchPPM, na.last = TRUE)
  rank <- integer(nrow(out))
  rank[o] <- stats::ave(seq_along(o), out$pid[o], FUN = seq_along)
  out$matchRank <- as.integer(rank)
  out
}

#' Link MS/MS scans to the individual peaks of their own file, then to the
#' features holding those peaks.
#'
#' @return The links, and the peak windows and feature memberships used.
#' @noRd
.f4f_link_individual <- function(scans, xcmsObj, ppm, rtTol, includeFilled,
                                 convert2RawRT, aligned) {
  modern <- .xcms_is_modern(xcmsObj)
  peaks <- if (modern) data.frame(xcms::chromPeaks(xcmsObj))
           else data.frame(xcmsObj@peaks)
  files <- .xcms_files(xcmsObj)
  peaks$filename <- .f4f_basename(files)[peaks$sample]
  peaks$cid <- seq_len(nrow(peaks))
  filled <- .f4f_is_filled(xcmsObj, nrow(peaks))

  # Windows are compared on the adjusted scale; rtmin and rtmax are reported
  # on the raw scale of the file, with the adjusted values alongside.
  adj_min <- peaks$rtmin
  adj_max <- peaks$rtmax
  if (convert2RawRT) {
    if (!is.null(aligned)) {
      peaks$rtminCorrected <- adj_min
      peaks$rtmaxCorrected <- adj_max
      for (s in unique(peaks$sample)) {
        r <- which(peaks$sample == s)
        p <- aligned[[s]]
        peaks$rtmin[r] <- .f4f_interp(adj_min[r], p$adjusted, p$raw)
        peaks$rtmax[r] <- .f4f_interp(adj_max[r], p$adjusted, p$raw)
      }
    } else {
      peaks$rtminCorrected <- NA
      peaks$rtmaxCorrected <- NA
    }
  }

  use <- if (includeFilled) rep(TRUE, nrow(peaks)) else !filled
  windows <- data.frame(.lo = adj_min - rtTol, .hi = adj_max + rtTol,
                        .mz = peaks$mz, sample = peaks$sample)
  windows$.lo[!use] <- NA
  m <- .f4f_join(scans, windows, ppm, by_sample = TRUE)

  # Peak to feature, once
  idx <- if (modern) xcms::featureDefinitions(xcmsObj)$peakidx
         else xcmsObj@groupidx
  p2g <- data.frame(.wi = as.integer(unlist(idx)),
                    grpid = rep(seq_along(idx), lengths(idx)))
  m <- merge(m, p2g, by = ".wi", sort = FALSE)
  m <- m[order(m$grpid, scans$fileid[m$.si], scans$seqNum[m$.si], m$.wi), ,
         drop = FALSE]

  out <- peaks[m$.wi, , drop = FALSE]
  out <- cbind(grpid = factor(m$grpid, levels = seq_along(idx)), out)
  out <- .f4f_match_columns(out, scans, m, peaks$mz[m$.wi])
  out$is_filled <- filled[m$.wi]
  out <- .f4f_rank(out, peaks$rt[m$.wi], adj_min[m$.wi], adj_max[m$.wi])
  rownames(out) <- NULL
  bounds <- peaks
  bounds$is_filled <- filled
  bounds$rtminWin <- windows$.lo
  bounds$rtmaxWin <- windows$.hi
  list(links = out, bounds = bounds, membership = idx)
}

#' Link MS/MS scans to feature windows.
#'
#' @return The links, and the feature windows used.
#' @noRd
.f4f_link_group <- function(scans, xcmsObj, ppm, rtTol, groupWindow) {
  modern <- .xcms_is_modern(xcmsObj)
  fd <- if (modern) data.frame(get_full_peak_width(
                       xcms::featureDefinitions(xcmsObj), xcmsObj = xcmsObj))
        else data.frame(get_full_peak_width(xcmsObj@groups, xcmsObj = xcmsObj))
  fd$grpid <- seq_len(nrow(fd))
  if (groupWindow == "median") {
    f <- if (modern) xcms::featureValues else xcms::groupval
    lo <- f(xcmsObj, method = "medret", value = "rtmin", intensity = "into")
    hi <- f(xcmsObj, method = "medret", value = "rtmax", intensity = "into")
    lo <- apply(lo, 1, stats::median, na.rm = TRUE)
    hi <- apply(hi, 1, stats::median, na.rm = TRUE)
  } else {
    lo <- fd$rtmin_full
    hi <- fd$rtmax_full
  }
  windows <- data.frame(.lo = lo - rtTol, .hi = hi + rtTol, .mz = fd$mzmed)
  m <- .f4f_join(scans, windows, ppm, by_sample = FALSE)
  m <- m[order(m$.wi, scans$pid[m$.si]), , drop = FALSE]

  out <- fd[m$.wi, , drop = FALSE]
  out <- .f4f_match_columns(out, scans, m, fd$mzmed[m$.wi])
  out$rtminWin <- lo[m$.wi]
  out$rtmaxWin <- hi[m$.wi]
  out <- .f4f_rank(out, fd$rtmed[m$.wi], lo[m$.wi], hi[m$.wi])
  rownames(out) <- NULL
  bounds <- fd
  bounds$rtminWin <- windows$.lo
  bounds$rtmaxWin <- windows$.hi
  list(links = out, bounds = bounds, membership = NULL)
}

#' Per-file linking diagnostics, as a JSON-safe list.
#'
#' For each file: its mapping status, the number of MS/MS scans, the share
#' linked and the median difference between the MS/MS time and the rtmed of
#' its best-ranked feature. A file whose median difference is above 10
#' seconds has a poor mapping and raises a warning.
#'
#' @noRd
.f4f_diagnostics <- function(grpm, scans, features, pa_files, sample_idx,
                             status, estimated) {
  best <- grpm[grpm$matchRank == 1L, , drop = FALSE]
  gid <- as.integer(as.character(best$grpid))
  best$.offset <- best$precurMtchRTadj - features$rtmed[match(gid, features$grpid)]
  out <- lapply(seq_along(pa_files), function(i) {
    n <- sum(scans$fileid == i)
    linked <- length(unique(grpm$pid[grpm$fileid == i]))
    off <- best$.offset[best$fileid == i]
    est <- estimated[[as.character(i)]]
    list(file = .f4f_basename(pa_files[i]), fileid = i,
         sample = if (is.na(sample_idx[i])) NA_integer_ else sample_idx[i],
         status = status[i], scans = n, linked = linked,
         share_linked = if (n) linked / n else NA_real_,
         median_offset = if (length(off)) stats::median(off) else NA_real_,
         estimate_method = if (is.null(est)) NA_character_ else est$method,
         landmarks = if (is.null(est)) NA_integer_ else as.integer(est$landmarks))
  })
  offsets <- vapply(out, `[[`, numeric(1), "median_offset")
  far <- !is.na(offsets) & abs(offsets) > 10
  if (any(far))
    warning("MS/MS times differ from feature rtmed by more than 10 seconds (median) in ",
            paste(vapply(out[far], `[[`, character(1), "file"), collapse = ", "),
            ". Check the retention time mapping of these files.", call. = FALSE)
  out
}
