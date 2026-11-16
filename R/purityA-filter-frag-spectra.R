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


#' @title Filter fragmentation spectra associated with an XCMS feature
#' @aliases filterFragSpectra
#' @description
#' **General**
#'
#' Flag and filter features based on signal-to-noise ratio, relative abundance, intensity threshold and purity of the precursor ion.
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> (optionally) xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'    +  --- *Older versions of XCMS* --- (mzML files) -> xcms.xcmsSet -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> frag4feature -> **filterFragSpectra** -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' @param pa object; purityA object
#' @param ilim numeric; min intensity of a peak
#' @param plim numeric; minimum MS2 precursor ion purity. MS3/MS4 spectra
#'   use fragment-quality thresholds without another precursor-purity
#'   requirement. Diagnostic MS3/MS4 purity estimates are retained for
#'   inspection and never determine whether spectra or peaks pass this filter.
#' @param ra numeric; minimum relative abundance of a peak
#' @param snr numeric; minimum signal-to-noise of a peak within each file
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged.
#' @param snmeth character; Method to calculate signal to noise ration (either median or mean)
#' @param cores numeric; Accepted for compatibility.
#'   Fragment filtering is serial.
#' @param param NULL, a complete \code{FilterFragSpectraParam}, or a named list
#'   of complete parameter objects. Parameter-set filtering requires one core.
#' @param allfrag logical; consider all retained fragmentation spectra when
#'   TRUE, or only feature-associated spectra when FALSE. In either case,
#'   \code{msLevel} selects the levels included in the processed result.
#' @param msLevel nonempty numeric vector of integer MS levels from 2 to 4.
#'   Defaults to MS2. Use \code{msLevel = 2:4} for all supported levels.
#'   Requested levels with no eligible spectra produce an empty result.
#'
#' Level selection narrows processed spectra and feature associations,
#' retaining source measurements and ancestry. Start alternatives from the
#' same linked object: filtering cannot restore excluded associations.
#' Parameter \code{select} assigns settings within the selected levels.
#' History records requested \code{msLevel} and actual \code{effective_msLevel}.
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
#' #====== msPurity ============================
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xcmsObj)
#' #pa <- filterFragSpectra(pa)
#'
#' ## Run from previously generated data
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "2_frag4feature_pa.rds", package="msPurity"))
#' pa <- filterFragSpectra(pa)
#'
#' @return Returns a purityA object. The flags are stored in pa@@fragSpectra: purity_pass_flag as a spectra
#' variable and the others as peak variables (stored as 1 and 0). \code{groupedSpectra(pa)}, and the legacy
#' pa@@grped_ms2 slot while \code{options(msPurity.legacySlots = TRUE)}, give the spectra matrices with the
#' following columns
#'
#' * snr: Signal to noise ratio (calculated at scan level)
#' * ra: Relative abundance (calculated at scan level)
#' * purity_pass_flag: Precursor ion purity flag (1 pass, 0 fail)
#' * intensity_pass_flag: Intensity flag (1 pass, 0 fail)
#' * snr_pass_flag: Signal-to-noise pass flag (1 pass, 0 fail)
#' * ra_pass_flag: Relative abundance pass flag (1 pass, 0 fail)
#' * pass_flag: Overall pass flag, all flags must pass for this to pass (1 pass, 0 fail)
#'
#' For MS3/MS4, \code{purity_pass_flag} is unknown (\code{NA});
#' \code{pass_flag} describes product-quality criteria only. Filtering an MS2
#' spectrum does not automatically exclude its linked MS3/MS4 descendants.
#' @md
#' @export
setMethod(f = "filterFragSpectra", signature = "purityA",
    definition = function(pa, ilim = 0, plim = 0.8, ra = 0, snr = 3,
        cores = 1, rmp = FALSE, snmeth = "median", allfrag = FALSE,
        param = NULL, msLevel = 2L) {
        pa <- .pa_update(pa)
        if (length(pa@filter_frag_params) && pa@filter_frag_params$rmp) {
            message("Fragmentation peaks have been previously filtered and ",
                "removed - function can't be performed")
            return(pa)
        }
        .pa_msn_validate_preparation(pa)
        processing_started <- Sys.time()
        selection <- .pa_filter_selection(pa, allfrag, msLevel)
        pa@grped_df <- selection$links
        if (length(pa@fragSpectra))
            pa@fragSpectra <- pa@fragSpectra[
                pa@fragSpectra$pid %in% selection$pids]

        settings <- list(ilim = ilim, plim = plim, ra = ra, snr = snr,
            snmeth = snmeth, rmp = rmp, allfrag = allfrag)
        if (!is.null(param)) {
            fallback <- c(settings, list(cores = cores, msLevel = msLevel,
                    effective_msLevel = selection$levels))
            return(.pa_filter_parameters(pa, param, fallback,
                    selection, processing_started))
        }

        processing_input <- .pa_processing_input(pa, "filterFragSpectra",
            allfrag = allfrag, source = selection$source)

        # Reset linked peak arrays using the previous filtering state.
        if (!allfrag)
            grped_ms2 <- lapply(.pa_grouped_legacy(pa), lapply, resetGrpedMS2)
        filter_frag_params <- c(settings, list(msLevel = msLevel,
                effective_msLevel = selection$levels))
        pa@filter_frag_params <- filter_frag_params
        if (nrow(pa@grped_df)) {
            rows <- match(pa@grped_df$pid, pa@spectra$pid)
            if (anyNA(rows))
                stop("Feature links refer to unavailable spectra.")
            assessment <- .pa_filter_purity(Spectra::msLevel(pa@spectra)[rows],
                pa@grped_df$inPurity, plim, source = pa@spectra[rows])
            pa@grped_df$purity_pass_flag <- assessment$flag
        }

        if (allfrag) {
            puritydf <- purityTable(pa)
            selected <- puritydf$pid %in% selection$pids
            puritydf$purity_pass_flag <- rep(NA, nrow(puritydf))
            puritydf$purity_pass_flag[selected] <-
                puritydf$inPurity[selected] > plim
            pa@spectra <- .pa_set_purity_vars(pa@spectra, puritydf)
            pa@fragSpectra <- .pa_frag_all(pa, filter_frag_params, puritydf,
                selection$pids)
            all_frag_scans <- if (.pa_legacy()) .pa_allfrag_legacy(pa)
            pa <- .pa_sync_legacy(pa, puritydf = puritydf,
                grped_ms2 = .pa_grouped_legacy(pa),
                all_frag_scans = all_frag_scans)
        } else {
            if (nrow(pa@grped_df)) {
                internal <- pa@grped_df
                internal$purity_pass_flag <- assessment$gate
                grped_ms2 <- plyr::dlply(internal, ~grpid,
                    setPuritySpectraGrp, list(grped_ms2 = grped_ms2))
            } else {
                grped_ms2 <- list()
            }

            grped_ms2 <- lapply(grped_ms2, function(group) {
                lapply(group, function(m) {
                    if (is.null(m) || !nrow(m)) return(NULL)
                    setFlagMatrix(m, filter_frag_params)
                })
            })

            pa@fragSpectra <- .pa_frag_from_grouped(pa, grped_ms2)
            if (length(pa@fragSpectra)) {
                rows <- match(pa@fragSpectra$pid, pa@grped_df$pid)
                pa@fragSpectra$purity_pass_flag <- assessment$flag[rows]
                pa@fragSpectra$msPurityStagePurity <- assessment$purity[rows]
                pa@fragSpectra$msPurityPurityStatus <- assessment$status[rows]
            }
            pa <- .pa_sync_legacy(pa,
                grped_ms2 = .pa_grouped_legacy(pa),
                all_frag_scans = .pa_empty_frame())
        }

        pa <- .pa_processing_invalidate(pa)
        effective <- c(filter_frag_params,
            list(cores = cores, effective_cores = 1L))
        .pa_processing_record(pa, "filterFragSpectra", effective,
            processing_started, processing_input)

    })

# Select metadata and association occurrences once for both filtering routes.
.pa_filter_selection <- function(pa, allfrag, msLevel) {
    if (!is.logical(allfrag) || length(allfrag) != 1L || anyNA(allfrag))
        stop("allfrag must be TRUE or FALSE.")
    if (!is.numeric(msLevel) || !length(msLevel) || anyNA(msLevel) ||
            any(!msLevel %in% 2:4))
        stop("msLevel must be a nonempty numeric vector of levels 2, 3 or 4.")
    links <- pa@grped_df
    if (allfrag) {
        source <- pa@spectra
    } else if (!nrow(links)) {
        source <- pa@fragSpectra[integer()]
    } else {
        rows <- match(unique(links$pid), pa@fragSpectra$pid)
        if (anyNA(rows)) stop("Feature links refer to unavailable spectra.")
        source <- pa@fragSpectra[rows]
    }
    levels <- Spectra::msLevel(source)
    if (anyNA(levels) || any(!levels %in% 2:4))
        stop("Fragment filtering requires known MS2-MS4 spectra.")
    source <- source[which(levels %in% msLevel)]
    pids <- if (length(source)) source$pid else integer()
    if (length(pids) != length(source))
        stop("Filtering source spectra require pid identifiers.")
    list(source = source, pids = pids,
        links = links[links$pid %in% pids, , drop = FALSE],
        levels = unique(Spectra::msLevel(source)))
}

.pa_filter_purity <- function(level, purity, threshold, typed = FALSE,
    source = NULL) {
    if (anyNA(level) || any(!level %in% 2:4))
        stop("Fragment filtering requires known MS2-MS4 spectra.")
    if (!length(purity)) purity <- rep(NA_real_, length(level))
    if (length(purity) != length(level))
        stop("Purity measurements do not match the selected spectra.")
    threshold <- rep(threshold, length.out = length(level))
    ms2 <- level == 2L
    flag <- rep(NA, length(level))
    flag[ms2] <- purity[ms2] > threshold[ms2]
    if (typed) {
        flag[ms2 & !is.na(threshold) & is.na(purity)] <- FALSE
        flag[ms2 & is.na(threshold)] <- TRUE
    }
    status <- rep("not_assessed", length(level))
    status[ms2] <- "assessed"
    status[ms2 & is.na(purity)] <- "missing"
    status[ms2 & is.na(threshold)] <- if (typed) "skipped" else
        "unavailable_threshold"
    purity[!ms2] <- NA_real_
    if (!is.null(source) && "msPurityStagePurity" %in%
            Spectra::spectraVariables(source)) {
        purity[!ms2] <- source$msPurityStagePurity[!ms2]
        status[!ms2] <- source$msPurityPurityStatus[!ms2]
    }
    gate <- flag
    gate[!ms2] <- TRUE
    list(flag = flag, gate = gate, status = status, purity = purity)
}


.pa_filter_parameters <- function(pa, param, fallback, selection, started) {
    source <- selection$source
    source_ids <- selection$pids
    resolved <- .pa_filter_resolve(source, param, fallback)
    levels <- resolved$metadata$msLevel
    input <- .pa_processing_input(pa, "filterFragSpectra",
        allfrag = fallback$allfrag, source = source)
    if (fallback$allfrag && length(source))
        source <- .pa_raw_scans(pa, source_ids)
    for (field in c("activationType", "scanFilterAnalyzer"))
        source[[field]] <- resolved$metadata[[field]]
    selected <- .pa_filter_source_peaks(source)
    pk <- selected$peaks
    peak_ids <- if (fallback$allfrag) .pa_source_peak_ids(source, selected)
    else lapply(pk, function(m) rep(NA_real_, nrow(m)))
    tab <- purityTable(pa)
    purity <- if (length(source_ids))
        tab$inPurity[match(source_ids, tab$pid)] else numeric()
    thresholds <- vapply(resolved$parameters, `[[`, numeric(1), "plim")
    assessment <- .pa_filter_purity(levels, purity, thresholds, typed = TRUE,
        source = source)
    mats <- lapply(seq_along(source), function(i) {
        if (!nrow(pk[[i]])) return(NULL)
        settings <- resolved$parameters[[i]]
        settings$rmp <- FALSE
        value <- setFlagMatrix(cbind(pk[[i]], assessment$gate[i]), settings)
        value <- value[, setdiff(.PA_FRAG_LEGACY_COLS,
                "purity_pass_flag"), drop = FALSE]
        colnames(value)[2] <- "intensity"
        value <- cbind(value, sid = peak_ids[[i]],
            original_peak_index = selected$original_positions[[i]])
        if (fallback$rmp)
            value <- value[value[, "pass_flag"] == 1, , drop = FALSE]
        value
    })
    columns <- setdiff(Spectra::spectraVariables(source),
        Spectra::peaksVariables(source))
    sv <- as.data.frame(Spectra::spectraData(source, columns = columns),
        optional = TRUE)
    sv$purity_pass_flag <- assessment$flag
    sv$single_peak <- vapply(pk, nrow, integer(1)) == 1L
    sv$msPurityEffectiveParameterId <- resolved$assignments$parameter_id
    sv$msPurityParameterSetNames <- resolved$assignments$parameter_set_names
    sv$msPurityParameterFallback <- resolved$assignments$fallback_used
    sv$msPurityOriginalPeakIndexBasis <- selected$basis
    sv$msPurityFilterInputPeakCount <- vapply(pk, nrow, integer(1))
    sv$msPurityStagePurity <- assessment$purity
    sv$msPurityPurityStatus <- assessment$status
    frag <- .msp_memory(sv, mats, c(.PA_FRAG_PEAK_VARS, "original_peak_index"))
    .pa_spectra_metadata(frag) <- .pa_spectra_metadata(source)
    pa@fragSpectra <- frag
    g <- pa@grped_df
    idx <- if (nrow(g)) match(g$pid, source_ids) else integer()
    g$purity_pass_flag <- assessment$flag[idx]
    pa@grped_df <- g
    settings <- fallback
    settings$configuration <- resolved[c("assignments", "parameter_sets",
            "parameter_objects", "unused_parameter_sets")]
    pa@filter_frag_params <- settings
    if (fallback$allfrag) {
        tab$purity_pass_flag <- if (nrow(tab))
            assessment$flag[match(tab$pid, source_ids)] else logical()
        pa@spectra <- .pa_set_purity_vars(pa@spectra, tab)
        pa <- .pa_sync_legacy(pa, puritydf = tab,
            grped_ms2 = .pa_grouped_legacy(pa),
            all_frag_scans = .pa_allfrag_legacy(pa))
    } else {
        pa@all_frag_scans <- .pa_empty_frame()
        pa <- .pa_sync_legacy(pa, grped_ms2 = .pa_grouped_legacy(pa))
    }
    pa <- .pa_processing_invalidate(pa)
    settings$effective_cores <- 1L
    .pa_processing_record(pa, "filterFragSpectra", settings, started, input)
}




# Empty raw-peak schema, including results with every peak removed.
.pa_empty_scan_peaks <- function(pid = integer()) {
    data.frame(sid = integer(),
        pid = if (is.null(pid)) integer() else pid[0],
        fileid = integer(), scan = factor(), mz = numeric(),
        i = numeric(), type = character())
}

getScanPeaks <- function(pa) {
    scans <- purityTable(pa)
    key <- paste(as.character(scans$fileid), scans$seqNum)
    if (anyNA(scans$fileid) || anyNA(scans$seqNum) ||
            anyNA(scans$pid) || anyDuplicated(key) || anyDuplicated(scans$pid))
        stop("Spectrum identities must be unique by pid and file position.")
    parts <- lapply(seq_along(pa@fileList), function(fileid) {
        sp <- .msp_read(as.character(pa@fileList[[fileid]]))
        header <- .msp_header(sp)
        peaks <- .msp_peaks(sp)
        index <- which(header$msLevel > 1L)
        n <- vapply(peaks[index], nrow, integer(1L))
        values <- do.call(rbind, peaks[index])
        if (is.null(values) || !nrow(values)) return(NULL)
        data.frame(fileid = fileid,
            scan = as.character(rep(header$seqNum[index], n)),
            mz = values[, 1L], i = values[, 2L])
    })
    peaks <- do.call(rbind, parts)
    if (is.null(peaks))
        return(.pa_empty_scan_peaks(scans$pid))
    # Preserve peak-row ids before excluding unassessed higher-level spectra.
    peaks$sid <- seq_len(nrow(peaks))
    peaks$pid <- scans$pid[match(paste(peaks$fileid, peaks$scan), key)]
    peaks <- peaks[!is.na(peaks$pid), , drop = FALSE]
    peaks$scan <- factor(peaks$scan, levels = unique(as.character(
        scans$seqNum)))
    peaks$type <- rep("scan", nrow(peaks))
    rownames(peaks) <- NULL
    peaks[, c("sid", "pid", "fileid", "scan", "mz", "i", "type")]
}

resetGrpedMS2 <- function(m){
  m[,1:2, drop=FALSE]
}

setPuritySpectraGrp <- function(x, pa){
  grpid <- as.character(unique(x$grpid))

  # pa is a purityA object or a list holding grped_ms2
  msms_l <- if (is.list(pa)) pa$grped_ms2[as.character(grpid)][[1]]
            else pa@grped_ms2[as.character(grpid)][[1]]

  purity_pass_flag <- x$purity_pass_flag

  x$subsetid <- 1:nrow(x)
  result <- plyr::dlply(x, ~subsetid, setPurityFlagMatrix, msms_l=msms_l)
  result <- unname(result)
  attr(result, "split_type") <- NULL
  attr(result, "split_labels") <- NULL
  return(result)
}



setPurityFlagMatrix <- function(grpinfo, msms_l){
  m <- msms_l[[grpinfo$subsetid]]
  if (is.null(m) || !nrow(m)) return(NULL)
  m <- cbind(m, grpinfo$purity_pass_flag)
  return(m)
}

setFlagMatrix <- function(x, filter_frag_params){

  snmeth <- filter_frag_params$snmeth
  i_thre <- filter_frag_params$ilim
  ra_thre <- filter_frag_params$ra
  snr_thre <- filter_frag_params$snr
  rmp <- filter_frag_params$rmp

  # check if columns defined for intensity
  if('i' %in% colnames(x)){
    intensityIdx <- which(colnames(x)=='i')
  }else{
    intensityIdx <- 2
  }

  if (snmeth=="median"){
    snr <- x[,intensityIdx]/median(x[,intensityIdx])
  }else if(snmeth=="mean"){
    snr <- x[,intensityIdx]/mean(x[,intensityIdx])
  }

  ra <- x[,intensityIdx]/max(x[,intensityIdx])*100

  intensity_pass_flag <- x[,intensityIdx]>i_thre
  ra_pass_flag <- ra>ra_thre
  snr_pass_flag <- snr>snr_thre

  x <- cbind(x, snr, ra, intensity_pass_flag, ra_pass_flag, snr_pass_flag)


  # check if columns available for sid and pid
  if('pid' %in% colnames(x)){
    col_order <-   c('sid', 'pid', 'fileid', 'scan', 'mz', 'i', 'snr', 'ra','type', 'purity_pass_flag', 'intensity_pass_flag','ra_pass_flag', 'snr_pass_flag', 'pass_flag')

  }else{
    colnames(x) <- c('mz', 'i','purity_pass_flag', 'snr', 'ra', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')
    col_order <-   c('mz', 'i', 'snr', 'ra', 'purity_pass_flag', 'intensity_pass_flag','ra_pass_flag', 'snr_pass_flag', 'pass_flag')
  }




  if (nrow(x)==1){
    pass_flag <- sum(unlist(x[,c('purity_pass_flag', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')]))==4
  }else{
    pass_flag <- rowSums(x[,c('purity_pass_flag', 'intensity_pass_flag', 'ra_pass_flag', 'snr_pass_flag')])==4
  }

  x <- cbind(x, pass_flag)
  # reoder so it is easier to read


  if (rmp){
    x <- x[x[,'pass_flag']==1,,drop=FALSE]
    if (nrow(x)==0){
      return(NULL)
    }
  }

  return(x[,col_order, drop=FALSE])
}

