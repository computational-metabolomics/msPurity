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



#' @title Assess the acquired precursor ion purity of MS/MS spectra (constructor)
#'
#' @description
#' ## General
#'
#' Given a vector of LC-MS/MS or DI-MS/MS mzML file paths calculate the precursor ion purity of
#' each MS/MS scan.
#'
#' The precursor ion purity represents the measure of the contribution of a selected precursor
#' peak in an isolation window used for fragmentation and can be used as away of assessing the
#' spectral quality and level of "contamination" of fragmentation spectra.
#'
#' The calculation involves dividing the intensity of the selected precursor peak by the total
#' intensity of the isolation window and is performed before and after the MS/MS scan of
#' interest and interpolated at the recorded time of the MS/MS acquisition.
#'
#' Additionally, isotopic peaks are annotated and omitted from the calculation,
#' low abundance peaks are removed that are thought to have minor contribution to the
#' resulting MS/MS spectra and the isolation efficiency of the mass spectrometer can be
#' used to normalise the intensities used for the calculation.
#'
#' The output is a purityA S4 class object (referred to as pa for convenience throughout
#' the manual). The object contains a slot (pa@@puritydf) where the details of the purity
#' assessments for each MS/MS scan. The purityA object can then be used for further processing
#' including linking the fragmentation spectra to XCMS features, averaging fragmentation,
#' database creation and spectral matching (from the created database).
#'
#'
#' ## Example LC-MS/MS processing workflow
#'
#' The purityA object can be used for further processing including linking the fragmentation spectra to XCMS features, averaging fragmentation, database creation and spectral matching (from the created database).
#' See below for an example workflow:
#'
#'  * Purity assessments
#'    +  (mzML files) -> **purityA** -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> (optionally) xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'    +  --- *Older versions of XCMS* --- (mzML files) -> xcms.xcmsSet -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> frag4feature -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' ## Isolation efficiency
#'
#' When the isolation efficiency of an MS instrument is known the peak intensities within an isolation window can be normalised for the precursor purity calculation. The isolation efficiency can be estimated by measuring a single precursor across a sliding window. See figure 3 from the original msPurity paper (Lawson et al 2017). This has been experimentally measured  for a Thermo Fisher Q-Exactive Mass spectrometer using 0.5 Da windows and can be set within msPurity by using msPurity::iwNormQE.5() as the input to the iwNormFunc argument.
#'
#' Other options to model the isolation efficiency the  gaussian isolation window msPurity::iwNormGauss(minOff=-0.5, maxOff = 0.5) or a R-Cosine window msPurity::iwNormRCosine(minOff=-0.5, maxOff=0.5). Where the minOff and maxOff can be altered depending on the isolation window size.
#'
#' A user can also define their own normalisation function. The only requirement of the function is that given a value between the minOff and maxOff a normalisation value between 0-1 is returned.
#'
#' ## Isolation windows
#'
#' Each MS2 event uses its recorded isolation bounds. If bounds are missing or
#' invalid, its purity is NA and a warning identifies the affected count.
#' Supply explicit lower and upper offsets when the acquisition settings are
#' known; no instrument-specific width is assumed. For example,
#' \code{offsets = c(0.65, 0.65)} is appropriate only for a known 1.3 Da window.
#' Resolved bounds and their sources are retained in processing history.
#'
#' @param fileList vector of mzML file paths, or a \code{Spectra},
#'        \code{MsExperiment} or \code{XcmsExperiment} object holding the
#'        spectra of one or more files. Objects are used as they are, so any
#'        filtering applied to them also applies to the purity calculation
#' @param mostIntense boolean; True if the most intense peak is used for calculation. Set to FALSE if the peak closest to mz value detailed in mzML meta data.
#' @param nearest logical; TRUE uses the nearest usable MS1 measurement on either
#' side. FALSE requires the reported preceding MS1 reference. When references are
#' absent throughout a file, the preceding MS1 is used. A later measurement never
#' substitutes for an unavailable preceding reference. Missing permitted
#' measurements produce NA purity; one usable measurement supports a single estimate.
#' @param offsets numeric; Two nonnegative lower and upper offsets
#' in Da, with positive total width. Explicit offsets set a window
#' around the selected precursor m/z. NA uses each event's recorded
#' absolute bounds. A selected precursor outside those bounds has
#' unavailable purity.
#' @param plotP boolean; If TRUE a plot of the purity is to be saved
#' @param plotdir vector; If plotP is TRUE plots will be saved to this directory
#' @param interpol character; type of interolation to be performed "linear" or "spline" (Spline option is only included for testing purposes,
#'                            linear should be used for all standard cases, isotope removal is also not available for the spline option)
#' @param iwNorm boolean; If TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5% (0.05)
#' @param isotopes boolean; TRUE if isotopes are to be removed
#' @param im matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#' @param mzRback character; deprecated and ignored. Raw data is read through Spectra, which uses the pwiz reader
#' @param ppmInterp numeric; Set the ppm tolerance for the precursor ion purity interpolation. i.e. the ppm tolerence between
#'                           the precursor ion found in the neighbouring scans.
#' @param cores numeric; Number of cores to use
#' @details MS3/MS4 spectra are retained with their acquisition metadata, but
#'   MS1-based precursor purity is calculated only for MS2. Higher-level
#'   feature association is enabled with `frag4feature(includeMSn = TRUE)`.
#'   That operation also estimates diagnostic MS3/MS4 purity from their parent
#'   fragmentation spectra; purity thresholds still apply only to MS2.
#'
#' @return Returns a purityA object (pa). The MS/MS scans and their purity are
#' stored as a \code{Spectra} object in the pa@@spectra slot, with the columns
#' below as spectra variables. \code{purityTable(pa)} returns them as a data
#' frame, which is also kept in the legacy pa@@puritydf slot while
#' \code{options(msPurity.legacySlots = TRUE)} (the default).
#'
#' The purity dataframe (**purityTable(pa)**) consists of the following columns:
#' * pid: unique id for MS/MS scan
#' * fileid: unique id for mzML file
#' * seqNum: scan number
#' * precursorIntensity: precursor intensity value as defined in the mzML file
#' * precursorMZ: precursor m/z value as defined in the mzML file
#' * precursorRT: precursor RT value as defined in the mzML file
#' * precursorScanNum: precursor scan number value as defined in mzML file
#' * id: unique id (redundant)
#' * filename: mzML filename
#' * precursorNearest: MS1 scan nearest to the MS/MS scan
#' * aMz: The m/z value in the "precursorNearest" MS1 scan which most closely matches the precursorMZ value provided from the mzML file
#' * aPurity: The purity score for aMz
#' * apkNm: The number of peaks in the isolation window for aMz
#' * iMz: The m/z value in the precursorNearest MS1 scan that is the most intense within the isolation window.
#' * iPurity: The purity score for iMz
#' * ipkNm: The number of peaks in the isolation window for iMz
#' * inPurity: The interpolated purity score (the purity score is calculated at neighbouring MS1 scans and interpolated at the point of the MS/MS acquisition)
#' * inpkNm: The interpolated number of peaks in the isolation window
#'
#' The remaining slots for purityA class include
#' * pa@@cores: The number of CPUs to be used for any further processing with this purityA object
#' * pa@@fileList: list of the mzML files that have been processed
#' * pa@@mzRback: The backend library used by mzR to extract information from the mzML file (e.g. pwiz)
#' * pa@@grped_df: If frag4feature has been performed, a dataframe of the grouped XCMS features linked to the associated fragmentation spectra precursor details is recorded here
#' * pa@@grped_ms2:  If frag4feature has been performed, a list of fragmentation spectra associated with each grouped XCMS feature is recorded here
#' * pa@@f4f_link_type: If frag4feature has been performed, the 'linking method' is recorded here, e.g. 'group' or 'individual'. Default is 'individual', see frag4feature documentation for more details
#' * pa@@av_spectra: if averageIntraFragSpectra, averageInterFragSpectra,  or averageAllFragSpectra have been performed, the average spectra is recorded here
#' * pa@@av_intra_params:  If averageIntraFragSpectra has been performed, the parameters are recorded here
#' * pa@@av_inter_params: if averageInterFragSpectra has been performed, the  parameters are recorded here]
#' * pa@@av_all_params: If averageAllFragSpectra has been performed, the parameters are recorded here
#' * pa@@db_path: If create_database has been performed, the resulting path to the database is recorded here
#'
#' @examples
#' filepths <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")
#' pa <- purityA(filepths)
#' @seealso \code{\link{assessPuritySingle}}
#' @md
#' @export
purityA <- function(fileList,
                    cores=1,
                    mostIntense=FALSE,
                    nearest=TRUE,
                    offsets=NA,
                    plotP=FALSE,
                    plotdir=NULL,
                    interpol="linear",
                    iwNorm=FALSE,
                    iwNormFun=NULL,
                    ilim=0.05,
                    mzRback='pwiz',
                    isotopes=TRUE,
                    im=NULL,
                    ppmInterp=7){

    processing_started <- Sys.time()

  # Spectra, MsExperiment and XcmsExperiment inputs are used as they are;
  # file paths are read once here.
  if (is(fileList, "MsExperiment")){
    fileList <- MsExperiment::spectra(fileList)
  }
  if (is(fileList, "Spectra")){
    raw <- fileList
    fileList <- unique(Spectra::dataOrigin(raw))
  }else{
    if((is.null(fileList)) || (all(fileList == "" ))){
      message("no file list")
      return(NULL)
    }
    raw <- .msp_read(fileList)
  }
  acquisition <- .pa_filter_metadata(raw)
  for (field in c("activationType", "scanFilterAnalyzer"))
    raw[[field]] <- acquisition[[field]]
  names(fileList) <- basename(fileList)
  .msp_deprecate_mzRback(mzRback)
  sps <- .msp_split_files(raw)

  requireNamespace('foreach')
  pa <- new("purityA", fileList = fileList, cores = cores, mzRback=mzRback)
  pa@params$processing <- .pa_processing_empty()
  pa@params$purityA <- list(cores = cores, mostIntense = mostIntense, nearest = nearest,
                            offsets = offsets, plotP = plotP, plotdir = plotdir,
                            interpol = interpol, iwNorm = iwNorm, iwNormFun = iwNormFun,
                            ilim = ilim, mzRback = mzRback, isotopes = isotopes, im = im,
                            ppmInterp = ppmInterp)

  # Check cores and choose if parallel or not (do or dopar)
  if(pa@cores<=1){
    operator <- foreach::'%do%'
  } else{
    cl1 <- parallel::makeCluster(pa@cores, type="SOCK")
    doSNOW::registerDoSNOW(cl1)
    operator <- foreach::'%dopar%'
  }



  # run parallel (or not) using foreach
  purityL <- operator(
    foreach::foreach(
      i = 1:length(pa@fileList),
      .packages = "Spectra",
      .export = ".pa_assess_purity_single"
    ),
    .pa_assess_purity_single(
      filepth = sps[[i]],
      mostIntense = mostIntense,
      nearest = nearest,
      offsets = offsets,
      plotP = plotP,
      plotdir = plotdir,
      interpol = interpol,
      iwNorm = iwNorm,
      iwNormFun = iwNormFun,
      ilim = ilim,
      mzRback = mzRback,
      isotopes = isotopes,
      im = im,
      ppmInterp = ppmInterp
    )
  )

  if(pa@cores>1){
      parallel::stopCluster(cl1)
  }

  names(purityL) <- seq(1, length(purityL))

  resolved <- lapply(seq_along(purityL), function(i) {
    settings <- purityL[[i]]$resolved
    settings$fileid <- as.integer(i)
    settings
  })
  names(resolved) <- names(purityL)
  puritydf <- plyr::ldply(lapply(purityL, `[[`, "purity"), .id = TRUE)

  if(nrow(puritydf)>0){
    colnames(puritydf)[1] <- "fileid"
    pid <- seq(1, nrow(puritydf))
    puritydf <- cbind(pid, puritydf)
  }

  o <- Spectra::dataOrigin(raw)
  pa@spectra <- .pa_spectra_new(raw, match(o, unique(o)), puritydf)
  pa <- .pa_sync_legacy(pa, puritydf = puritydf)
  pa@cores = cores
  effective_parameters <- pa@params$purityA
  effective_parameters$effective_cores <- max(1L, as.integer(pa@cores))
  effective_parameters$resolved <- resolved
  pa <- .pa_processing_record(pa, "purityA", effective_parameters,
    processing_started,
    list(previous_steps = character(), spectra = .pa_processing_refs(raw)))

  return(pa)
}


#' @title Assess the purity of a single LC-MS/MS or DI-MS/MS file
#'
#' @description
#'
#' Given a filepath to an mzML file the precursor purity for any MS/MS scans
#' will be outputed into a dataframe
#'
#' @param filepth character; mzML file path for MS/MS spectra
#' @param fileid numeric; adds a fileid column (primarily for internal use for msPurity)
#' @param mostIntense boolean; True if the most intense peak is used for calculation. False if the centered peak is used
#' @param nearest logical; TRUE uses the nearest usable MS1 on either side.
#' FALSE requires the reported preceding reference. Missing permitted measurements
#' produce NA purity; a later measurement is not substituted.
#' @param offsets numeric; Two nonnegative lower and upper offsets
#' in Da, with positive total width. Explicit offsets set a window
#' around the selected precursor m/z. NA uses each event's recorded
#' absolute bounds. A selected precursor outside those bounds has
#' unavailable purity.
#' @param cores numeric; Number of cores to use
#' @param plotP boolean; If TRUE a plot of the purity is to be saved
#' @param plotdir vector; If plotP is TRUE plots will be saved to this directory
#' @param interpol character; Type of interolation to be performed "linear", "spline" or "none"
#' @param iwNorm boolean; If TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5% (0.05)
#' @param mzRback character; deprecated and ignored. Raw data is read through Spectra, which uses the pwiz reader
#' @param isotopes boolean; TRUE if isotopes are to be removed
#' @param im matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#' @param ppmInterp numeric; Set the ppm tolerance for the precursor ion purity interpolation. i.e. the ppm tolerence between
#'                           the precursor ion found in the neighbouring scans.
#' @return a dataframe of the purity score of the ms/ms spectra
#'
#' @examples
#' filepth <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")
#'
#' puritydf <- assessPuritySingle(filepth)
#' @seealso \code{\link{purityA}}
#' @export
assessPuritySingle <- function(filepth,
                               fileid=NA,
                               mostIntense=FALSE,
                               nearest=TRUE,
                               offsets=NA,
                               cores=1,
                               plotP=FALSE,
                               plotdir=NULL,
                               interpol="linear",
                               iwNorm=FALSE,
                               iwNormFun=NULL,
                               ilim=0,
                               mzRback='pwiz',
                               isotopes=TRUE,
                               im=NULL,
                               ppmInterp=7){
  .pa_assess_purity_single(
    filepth = filepth, fileid = fileid, mostIntense = mostIntense,
    nearest = nearest, offsets = offsets, cores = cores, plotP = plotP,
    plotdir = plotdir, interpol = interpol, iwNorm = iwNorm,
    iwNormFun = iwNormFun, ilim = ilim, mzRback = mzRback,
    isotopes = isotopes, im = im, ppmInterp = ppmInterp
  )$purity
}

# The public function returns only the historical purity table. The constructor
# also retains settings resolved by this same assessment, without reading twice.
.pa_assess_purity_single <- function(
  filepth,
  fileid = NA,
  mostIntense = FALSE,
  nearest = TRUE,
  offsets = NA,
  cores = 1,
  plotP = FALSE,
  plotdir = NULL,
  interpol = "linear",
  iwNorm = FALSE,
  iwNormFun = NULL,
  ilim = 0,
  mzRback = "pwiz",
  isotopes = TRUE,
  im = NULL,
  ppmInterp = 7
) {
  #=================================
  # Load in files and initial setup
  #=================================
  # Get the mzR dataframes
  automatic_offsets <- .pa_automatic_offsets(offsets)
  sp <- .msp_as_spectra(filepth)
  if (is(filepth, "Spectra")){
    filepth <- Spectra::dataOrigin(sp)[1]
  }
  resolved <- list(
    dataOrigin = as.character(filepth),
    status = "no_ms2",
    offsets = c(NA_real_, NA_real_),
    offsets_source = "not_used",
    offsets_method = NA_character_,
    iwNorm = iwNorm,
    iwNormFun = NULL,
    iwNormFun_source = "not_used"
  )
  mrdf <- getmrdf(sp, mzRback)
  if (is.null(mrdf) || !any(mrdf$msLevel == 2L)) {
    message(paste("No MS/MS spectra for file: ", filepth))
    return(list(purity = NULL, resolved = resolved))
  }

  # having id column makes things easier to track
  mrdf$id <- seq(1, nrow(mrdf))

  # add the fileid info
  if(!is.na(fileid)){
    mrdf$fileid <- rep(fileid, nrow(mrdf))
  }

  # add filename
  mrdf$filename <- basename(filepth)

  # Get a shortened mzR dataframe
  mrdfshrt <- mrdf[mrdf$msLevel==2,][,c("seqNum","acquisitionNum","precursorIntensity",
                                        "precursorMZ", "precursorRT",
                                        "precursorScanNum", "id", "filename", "retentionTime")]

  if (!any(mrdf$msLevel == 1L)) {
    message("No MS1 measurements; only MS2 data cannot supply precursor purity.")
    mrdfshrt[ , c("precursorNearest", "aMz", "aPurity", "apkNm", "iMz", "iPurity", "ipkNm", "inPkNm", "inPurity")] <- NA
    resolved$status <- "no_ms1"
    resolved$observations <- .pa_ms1_observations(mrdfshrt, mrdf,
        rep(list(integer()), nrow(mrdfshrt)), "no_ms1")
    return(list(purity = mrdfshrt, resolved = resolved))
  }

  # get scans (list of mz and i) from mzR
  scans <- getscans(sp, mzRback)

  selected <- sp[Spectra::msLevel(sp) == 2L]
  if (automatic_offsets) {
    windows <- .msp_isolation_settings(selected)
  } else {
    windows <- .pa_isolation_windows(Spectra::precursorMz(selected),
        Spectra::isolationWindowLowerMz(selected),
        Spectra::isolationWindowUpperMz(selected),
        Spectra::isolationWindowTargetMz(selected), offsets)
  }
  windows$seqNum <- mrdfshrt$seqNum
  windows$acquisitionNum <- mrdfshrt$acquisitionNum
  invalid <- windows$status != "valid"
  if (any(invalid))
    warning(sum(invalid), " MS2 spectra have missing or invalid isolation ",
        "bounds; their purity is unavailable. Supply offsets when known.",
        call. = FALSE)
  isolation <- .pa_isolation_history(windows, iwNorm, iwNormFun)
  resolved[names(isolation)] <- isolation
  resolved$status <- "assessed"
  minoff <- maxoff <- NA_real_




  # For n MS1 scans before and after the MS2 scan. For linear interpolation
  # only two points needed. More needed for spline
  if((interpol=="linear") || (interpol=="none")){
    # if linear interpolation we only need the purity before and after.
    # If no interpolation required we still need to find the nearest ms1 scan
    # so the nump is 1 as well
    nump <- 1
  }else if (interpol=="spline"){
    # if performing spline we need multiple data points
    nump <- 10
  }

  # Get the precursor (ms1) scans
  usable <- vapply(scans, function(peaks) {
    nrow(peaks) > 0L && all(is.finite(peaks)) &&
        all(peaks[, 1L] > 0 & peaks[, 2L] >= 0)
  }, logical(1L))
  prec_scans <- get_prec_scans(mrdf, nump, usable)
  nearest_position <- vapply(prec_scans, `[[`, integer(1L), "nearest")
  mrdfshrt$precursorNearest <- mrdf$seqNum[nearest_position]
  reference_position <- match(mrdfshrt$precursorScanNum, mrdf$acquisitionNum)
  permitted <- !is.na(reference_position) &
      mrdf$msLevel[reference_position] == 1L & usable[reference_position] &
      reference_position < which(mrdf$msLevel == 2L)
  permitted[is.na(permitted)] <- FALSE
  reference_position[!permitted] <- NA_integer_
  mrdfshrt$.purity_position <- if (nearest) nearest_position else reference_position
  unavailable <- is.na(mrdfshrt$.purity_position)
  if (any(unavailable))
    warning(sum(unavailable), " MS2 spectra have no usable permitted MS1 ",
        "measurement; their purity is unavailable.", call. = FALSE)


  #=====================================
  # Initial purity calculation
  #=====================================
  # Get the initial purity score and target mz associated with every ms2 scan.
  # Calculate the purity on the nearest ms1 scan. This function also determines
  # which is the target peak used for interpolations
  initial <- plyr::ddply(mrdfshrt, ~ seqNum, get_init_purity,
                           scans=scans,
                           minoff=minoff,
                           maxoff=maxoff,
                           nearest=nearest,
                           mostIntense=mostIntense,
                           iwNorm=iwNorm,
                           iwNormFun=iwNormFun,
                           ilim=ilim,
                           isotopes=isotopes,
                           im = im, windows = windows)

    # Keep the input key once; ddply also returns its grouping column.
    if (nrow(initial) != nrow(mrdfshrt) ||
            !identical(as.numeric(initial$seqNum),
                as.numeric(mrdfshrt$seqNum)))
        stop("Initial purity results do not match the spectrum positions.")
    mrdfshrt <- cbind(mrdfshrt,
        initial[, names(initial) != "seqNum", drop = FALSE])

  # If no interpolation required then all purity calculations are completed
  used <- lapply(seq_len(nrow(mrdfshrt)), function(i) {
    if (windows$status[i] != "valid" || is.na(mrdfshrt$.purity_position[i]))
        integer() else mrdfshrt$.purity_position[i]
  })
  if(interpol=="none"){
    resolved$observations <- .pa_ms1_observations(mrdfshrt, mrdf, used,
        ifelse(windows$status == "valid", "initial", windows$status))
    mrdfshrt$.purity_position <- NULL
    return(list(purity = mrdfshrt, resolved = resolved))
  }

  #=====================================
  # Interpolation
  #=====================================
  # The interpolation can be time consuming if using spline interpolation,
  # so can be run in parallel
  if(cores>1){
    requireNamespace(parallel)
    clust<-parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(clust)
    pBool = TRUE
  }else{
    pBool = FALSE
  }

  interDF <- plyr::ddply(mrdfshrt, ~ seqNum, .parallel = pBool,
                                   get_interp_purity, # FUNCTION
                                   scan_peaks=scans,
                                   ppm=ppmInterp,
                                   ms2=mrdf[mrdf$msLevel==2,]$seqNum,
                                   prec_scans=prec_scans,
                                   minoff=minoff,
                                   maxoff=maxoff,
                                   mostIntense=mostIntense,
                                   plotP=plotP,
                                   plotdir=plotdir,
                                   interpol=interpol,
                                   nearest=nearest,
                                   iwNorm=iwNorm,
                                   iwNormFun=iwNormFun,
                                   ilim=ilim,
                                   isotopes=isotopes,
                                   im = im, windows = windows)

  if(cores>1){
    requireNamespace(parallel)
    parallel::stopCluster(clust)
  }

  if(interpol=="spline"){
    mrdfshrt$inPurity <- interDF$inPurity
  }else{
    mrdfshrt$inPkNm <- interDF$inPkNm
    mrdfshrt$inPurity <- interDF$inPurity
  }


  resolved$observations <- .pa_ms1_observations(mrdfshrt, mrdf,
      interDF$ms1_used, ifelse(windows$status == "valid", interpol, windows$status))
  mrdfshrt$.purity_position <- NULL
  return(list(purity = mrdfshrt, resolved = resolved))

}


get_init_purity <- function(ms2h, scans, minoff, maxoff, nearest,
                            mostIntense, iwNorm, iwNormFun, ilim, isotopes, im,
                            windows = NULL) {

  if (!is.null(windows)) {
    window <- .pa_purity_window(ms2h, windows, iwNormFun)
    if (!window$valid)
      return(c(aMz = NA_real_, aPurity = NA_real_, apkNm = NA_real_,
          iMz = NA_real_, iPurity = NA_real_, ipkNm = NA_real_))
    minoff <- window$minoff
    maxoff <- window$maxoff
    iwNormFun <- window$iwNormFun
  }

  #========================================
  # Get the scan to perform calculations on
  #========================================
  if (".purity_position" %in% names(ms2h)) {
    precScn <- ms2h$.purity_position
  } else if (nearest) {
    precScn <- ms2h$precursorNearest
  }else{
    precScn <- ms2h$precursorScanNum
  }

  # What the mzML says the mz and intensity of the target is
  precMZ <- ms2h$precursorMZ
  precI <- ms2h$precursorIntensity

  mzmin <-  precMZ - minoff
  mzmax <-  precMZ + maxoff

  if (length(precScn) != 1L || is.na(precScn) || precScn < 1L ||
          precScn > length(scans) || !nrow(scans[[precScn]]))
    return(c(aMz = NA_real_, aPurity = NA_real_, apkNm = NA_real_,
        iMz = NA_real_, iPurity = NA_real_, ipkNm = NA_real_))
  allscans <- scans[[precScn]]

  #=======================================
  # Get centered target peak
  #=======================================
  # Get target based on matching to be within isolation window of precursorMZ
  # and 1 decimal of intensity
  target_approx <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax)
                            & (round(allscans[,2],1) == round(precI,1)),]

  # no precusor within isolation window and itensity
  if(length(target_approx)==0){
    # get nearest match for just mz
    target_approx <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax),]
    if(length(target_approx)==0){
      return(c("aMz"= NA, "aPurity" = 0, "apkNm" = NA,
               "iMz" = NA, "iPurity" = 0, "ipkNm" = NA ))
    }
  }

  if (length(target_approx)==2) {
    target <- target_approx
  }else{
    # Get the mz that is closest to the target
    target <- target_approx[which.min(abs(target_approx[,1]- precMZ)),]
  }

  # The centered peak
  aMz <- unname(target[1])

  #==================================================
  # Get most intense peak within isolation window
  #==================================================
  subp <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax),  ]

  if(is.vector(subp)){
    iTarget <- subp
  }else if (nrow(subp)==0){
    # nothing within in the mz tolerance range
    return(c("aMz"=aMz, "aPurity" = 0, "iMz" = NA, "iPurity" = 0 ))
  }else{
    iTarget <- subp[order(subp[,2], decreasing = TRUE),][1,]
  }

  # the most intense peak
  iMz <- unname(iTarget[1])

  #==================================================
  # Calculate initial purity
  #==================================================
  if(is.vector(subp)){
    aPurity = iPurity = apkNm = ipkNm =  1
  } else if (nrow(subp)==1){
    aPurity = iPurity = apkNm = ipkNm =  1

  } else {
    pouta <- pcalc(peaks=subp, mzmin=mzmin, mzmax=mzmax, mztarget=aMz, ppm=NA,
                   iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim, isotopes=isotopes,
                   im=im)
    aPurity <- unname(pouta[1])
    apkNm <- unname(pouta[2])
    pouti <- pcalc(peaks=subp, mzmin=mzmin, mzmax=mzmax, mztarget=iMz, ppm=NA,
                   iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim, isotopes=isotopes,
                   im=im)
    iPurity <- unname(pouti[1])
    ipkNm <- unname(pouti[2])

  }
  fileinfo <- c("aMz"=aMz, "aPurity" = aPurity, "apkNm" = apkNm,
                "iMz" = iMz, "iPurity" = iPurity, "ipkNm" = ipkNm )
  return(fileinfo)

}

get_interp_purity <- function(rowi, scan_peaks, prec_scans, ms2, ppm,
                              minoff, maxoff, mostIntense=FALSE, plotP=FALSE,
                              plotdir=NULL, interpol, nearest=TRUE,
                              iwNorm, iwNormFun, ilim, isotopes, im,
                              windows = NULL) {

  if (!is.null(windows)) {
    window <- .pa_purity_window(rowi, windows, iwNormFun)
    if (!window$valid || is.na(rowi$.purity_position))
      return(data.frame(inPkNm = NA_real_, inPurity = NA_real_,
          ms1_used = I(list(integer()))))
    minoff <- window$minoff
    maxoff <- window$maxoff
    iwNormFun <- window$iwNormFun
  }

  # get the region of interest for scans
  scanids <- prec_scans[which(ms2==rowi$seqNum)][[1]]
  roi_scns <- scan_peaks[c(scanids$pre, scanids$post)]

  if(interpol=="linear"){
    purity <- linearPurity(rowi, scan_peaks, minoff, maxoff, ppm, scanids,
                           nearest, mostIntense, iwNorm, iwNormFun,
                           ilim, plotP, plotdir, isotopes, im)
  }else if (interpol=="spline"){
    purity <- splinePurity(rowi, roi_scns, rowi$precursorMZ - minoff,
                           rowi$precursorMZ + maxoff, ppm,
                           mostIntense, scanids, plotP, plotdir)
    attr(purity, "ms1_used") <- if (is.finite(purity))
        c(scanids$pre, scanids$post) else integer()
  }
  if (!is.null(windows)) {
    used <- attr(purity, "ms1_used")
    purity <- if (interpol == "linear") as.data.frame(as.list(purity)) else
        data.frame(inPkNm = NA_real_, inPurity = purity)
    purity$ms1_used <- I(list(used))
  }
  return(purity)

}

linearPurity <- function(rowi, scan_peaks, minoff, maxoff, ppm, scanids,
                         nearest, mostIntense, iwNorm, iwNormFun, ilim,
                         plotP, plotdir, isotopes, im){
  finish <- function(value, used) {
    if (".purity_position" %in% names(rowi)) attr(value, "ms1_used") <- used
    value
  }
  if (".purity_position" %in% names(rowi)) {
    scn1 <- rowi$.purity_position
  } else if (nearest) {
    scn1 <- rowi$precursorNearest
  }else{
    scn1 <- rowi$precursorScanNum
  }

  if(mostIntense){
    purity1 <- rowi$iPurity
    mztarget1 <- rowi$iMz
    pknm1 <- rowi$ipkNm
  }else{
    purity1 <- rowi$aPurity
    mztarget1 <- rowi$aMz
    pknm1 <- rowi$apkNm
  }

  precMZ <- rowi$precursorMZ

  mzmin <-  precMZ - minoff
  mzmax <-  precMZ + maxoff


  # Calculate ppm tolerance to look either side of the original target value
  targetMaxMZ <- mztarget1 + (mztarget1*0.000001)*ppm
  targetMinMZ <- mztarget1 - (mztarget1*0.000001)*ppm


  # Then calculate purity for the other nearest MS1 scan (purity2)
  if (scn1 %in% scanids$pre) {
    scn2 <- scanids$post
  }else{
    scn2 <- scanids$pre
  }

  # The last scan (nothing to interpolate from)
  if((length(scn2)==0) || (is.na(mztarget1))){
    return(finish(c("inPkNm" = pknm1, "inPurity" = purity1), scn1))
  }

  # Calculate purity for the other nearest MS1 scan
  allp <- scan_peaks[[scn2]]

  pout <- pcalc(peaks=allp, mzmin=mzmin, mzmax=mzmax, mztarget=mztarget1,
                iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim,
                targetMinMZ = targetMinMZ, targetMaxMZ = targetMaxMZ,
                isotopes = isotopes, im = im)

  purity2 <- pout[1]
  pknm2 <- pout[2]



  ###########################################
  # Interpolate for purity
  ###########################################
  # Cant interpolate if only have 1 NA value
  if(is.na(purity2)){
    inPurity <- purity1
    fpurity <- NULL
  }else if (is.na(purity1)){
    inPurity <- purity2
    fpurity <- NULL
  }else{
    # Create linear functions for poth purity and the number of peaks
    coordinates <- if (is.null(scanids$coordinates)) c(scn1, scn2) else
        scanids$coordinates[c(scn1, scn2)]
    fpurity <- approxfun(coordinates, c(purity1, purity2))
    inPurity <- fpurity(rowi$seqNum)
  }

  ###########################################
  # Interpolate for peak number
  ###########################################
  if(is.na(pknm2)){
    inPkNm <- pknm1
  }else if (is.na(pknm1)){
    inPkNm <- pknm2
  }else{
    # Create linear functions for poth purity and the number of peaks
    coordinates <- if (is.null(scanids$coordinates)) c(scn1, scn2) else
        scanids$coordinates[c(scn1, scn2)]
    fpk <- approxfun(coordinates, c(pknm2, pknm2))
    inPkNm <- fpk(rowi$seqNum)
  }

  # save a plot of the interpolation
  if(plotP){
    df <- data.frame("idx"=c(scn1, scn2),"purity"=c(purity1, purity2))
    intername <- paste(rowi$file, "_", rowi$id, "_interpolate_plot_",
                       mztarget1,".png", sep="")
    name <- file.path(plotdir, intername)
    plotPurity(df, rowi$seqNum, name, fpurity)
  }

  return(finish(c("inPkNm" = inPkNm, "inPurity" = inPurity), c(scn1, scn2)))
}

get_prec_scans <- function(mrdf, num, usable = rep(TRUE, nrow(mrdf))) {

  ms1 <- which(mrdf$msLevel == 1L & usable)
  ms2 <- which(mrdf$msLevel == 2L)

  prec_scans <- plyr::alply(ms2, 1, function(scan2){
    post <- ms1[ms1>scan2]
    pre <- ms1[ms1<scan2]

    if(length(pre)>num){
      pre <- pre[(length(pre)-(num-1)):length(pre)]
    }

    if(length(post)>num){
      post <- post[1:num]
    }

    nearest <- if (length(ms1))
        ms1[which.min(abs(mrdf$seqNum[ms1] - mrdf$seqNum[scan2]))] else NA_integer_
    return(list(pre = pre, post = post, scan = scan2, nearest = nearest,
        coordinates = mrdf$seqNum))
  })

  return(prec_scans)

}

# Get the Data
getmrdf <- function(files, backend='pwiz'){
  # files can be file paths or a Spectra object; either way one header per file
  if (is(files, "Spectra")){
    sps <- .msp_split_files(files)
    files <- vapply(sps, function(s) Spectra::dataOrigin(s)[1], "")
  }else{
    sps <- NULL
  }
  mrdf <- NULL

  for(i in 1:length(files)){
    #message(paste("processing file:" ,i))
    mrdfn <- .msp_header(if (is.null(sps)) .msp_read(files[i]) else sps[[i]])
    if(length(unique(mrdfn$msLevel))<2){
      if (unique(mrdfn$msLevel)==1){
        message("only MS1 data")
        next
      }
      #else{
      #  message("only fragmentation data")
      #  next
      #}
    }else{
      if(length(unique(mrdfn$precursorScanNum))<2){
        # Note: will be of length 1 even if no scans associated because
        # the mrdf will be zero for not assigned
        message("MS2 data has no associated scan data, will use most recent full scan for information")
        mrdfn  <- missing_prec_scan(mrdfn)
      }

    }


    #mrdfn$fileid <- rep(i,nrow(mrdfn))
    mrdfn$filename <- rep(basename(files[i]),nrow(mrdfn))
    mrdfn$precursorRT <- NA
    # precursorScanNum matches to the acuisitionNum, get row matching row number and relevant retntion time
    mrdfn[mrdfn$msLevel==2,]$precursorRT <- mrdfn[match(mrdfn[mrdfn$msLevel==2,]$precursorScanNum, mrdfn$acquisitionNum),]$retentionTime

    if(!is.data.frame(mrdf)){
      mrdf <- mrdfn
    }else{
      mrdf <- rbind(mrdf,mrdfn)
    }
  }
  return(mrdf)
}

missing_prec_scan <- function(mrdfn){
  for(i in 1:nrow(mrdfn)){
    if(mrdfn[i,]$msLevel>1){
      # Find most recent ms1 level scan
      for (x in i:1){
        if(mrdfn[x,]$msLevel==1){
          mrdfn[i,]$precursorScanNum = mrdfn[x,]$acquisitionNum
          break
        }
      }
    }
  }
  return(mrdfn)
}

getscans <- function(files, backend='pwiz'){
  # files can be file paths or a Spectra object. One file gives a list of
  # peak matrices; several give one such list per file.
  sps <- if (is(files, "Spectra")) .msp_split_files(files)
         else lapply(files, .msp_read)
  scan_peaks <- lapply(sps, .msp_peaks)
  if(length(scan_peaks)==1){
    return(scan_peaks[[1]])
  }
  return(scan_peaks)
}

# MSMSperMS <- function(filepths){
#   l <- lapply(filepths, MSMSperMSsingle)
#   return(plyr::ldply(l, summary))
# }
#
#
# MSMSperMSsingle <- function(filepth){
#   mr <- mzR::openMSfile(filepth)
#   mrdf <- mzR::header(mr)
#
#   m <- mrdf$msLevel
#
#   # First remove repeating ms1 elements
#   c <- 1
#   idx <- c()
#   for (i in 1:length(m)){
#
#     if(i==length(m)){
#       break
#     }
#     if(m[i]==1 && m[i+1]==1) {
#       idx[c] <- i
#       c = c+1
#     }
#   }
#   if(!is.null(idx)){
#     m <- m[-idx]
#   }
#
#
#   # then count how many ms2 are present at each ms1---ms2ms2...ms2---ms1 interval
#   ms2all <- c()
#   c <- 0
#   cidx <- 1
#   for(i in 1:length(m)){
#
#     if(m[i]==2) {
#       c = c+1
#     }
#
#     if(m[i]==2 && (m[i+1]==1 || is.na(m[i+1]))) {
#       ms2all[cidx] <- c
#       cidx = cidx+1
#       c = 0
#     }
#
#   }
#
#   return(ms2all)
#
# }

plotPurity <- function(df, inter, name, f){

  dfm <- reshape2::melt(df, id.vars = c("purity", "idx"))
  png(filename = name)
  p1 <- ggplot2::ggplot(dfm, ggplot2::aes(x=idx, y=purity)) +
    ggplot2::geom_point(ggplot2::aes(colour=variable),shape=20, colour="red") +
    ggplot2::geom_vline(xintercept = inter, colour="red", linetype = "longdash")+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::theme(text = ggplot2::element_text(size=23))+
    ggplot2::labs(x="scan" , y="precursor purity score") +
    ggplot2::stat_function(fun=f, colour="blue")+
    ggplot2::ylim(0, 1)+
    ggplot2::scale_color_manual(
            "Legend Title\n",labels = c("EIC feature", "+/- 0.5 Contamination"),
             values = c("red", "black") )
  print(p1)
  dev.off()
  dfm$pick <- inter
  saveRDS(f, paste(name, ".rds", sep=""))
  write.csv(dfm, paste(name, ".csv", sep=""), row.names = FALSE)
}
