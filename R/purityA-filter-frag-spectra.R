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
#' @param plim numeric; min precursor ion purity of the associated precursor for fragmentation spectra scan
#' @param ra numeric; minimum relative abundance of a peak
#' @param snr numeric; minimum signal-to-noise of a peak within each file
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged.
#' @param snmeth character; Method to calculate signal to noise ration (either median or mean)
#' @param cores numeric; Accepted for compatibility.
#'   Fragment filtering is serial.
#' @param allfrag boolean; Whether to filter on all fragmentation spectra or just the fragmentation spectra grouped to XCMS feature
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
#' @md
#' @export
setMethod(f="filterFragSpectra", signature="purityA",
    definition = function(pa, ilim = 0, plim = 0.8, ra = 0, snr = 3,
        cores = 1, rmp = FALSE, snmeth = "median", allfrag = FALSE) {
            pa <- .pa_update(pa)

            if ((!length(pa@filter_frag_params)==0) && (pa@filter_frag_params$rmp)){
              message('Fragmentation peaks have been previously filtered and removed - function can\'t be performed')
              return(pa)
            }


        processing_started <- Sys.time()
        processing_input <- .pa_processing_input(pa, "filterFragSpectra",
            allfrag = allfrag)

            filter_frag_params = list()
            filter_frag_params$ilim = ilim
            filter_frag_params$plim = plim
            filter_frag_params$ra = ra
            filter_frag_params$snr = snr
            filter_frag_params$snmeth = snmeth
            filter_frag_params$rmp = rmp
            filter_frag_params$allfrag = allfrag

            # The grouped spectra as they were before any filtering (in case
            # filterFragSpectra has already been run)
            grped_ms2 <- lapply(.pa_grouped_legacy(pa), lapply, resetGrpedMS2)

            pa@filter_frag_params <- filter_frag_params

            # Calculate and add flags to matrix
            # Add the purity flag
            if (nrow(pa@grped_df)) {
                pa@grped_df$purity_pass_flag <- pa@grped_df$inPurity > plim
                grped_ms2 <- plyr::dlply(pa@grped_df, ~grpid,
                    setPuritySpectraGrp, list(grped_ms2 = grped_ms2))
            } else {
                grped_ms2 <- list()
            }

            # calculate snr, ra and combine all flags
            grped_ms2 <- lapply(grped_ms2, lapply, setFlagMatrix, filter_frag_params=filter_frag_params)

            # if allfrag is TRUE calculate the flags for all fragmentation spectra; these then hold
            # the flags of the grouped spectra as well
            if (allfrag) {
                puritydf <- purityTable(pa)
                puritydf$purity_pass_flag <- puritydf$inPurity > plim
                pa@spectra <- .pa_set_purity_vars(pa@spectra, puritydf)
                pa@fragSpectra <- .pa_frag_all(pa, filter_frag_params, puritydf)
                all_frag_scans <- if (.pa_legacy()) {
                    .pa_allfrag_legacy(pa)
                }
                pa <- .pa_sync_legacy(pa, puritydf = puritydf,
                    grped_ms2 = grped_ms2, all_frag_scans = all_frag_scans)
            } else {
                pa@fragSpectra <- .pa_frag_from_grouped(pa, grped_ms2)
                pa <- .pa_sync_legacy(pa, grped_ms2 = grped_ms2)
            }

        pa <- .pa_processing_invalidate(pa)
        effective <- c(filter_frag_params,
            list(cores = cores, effective_cores = 1L))
        .pa_processing_record(pa, "filterFragSpectra", effective,
            processing_started, processing_input)

          }
)


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
  m <- cbind(msms_l[[grpinfo$subsetid]], grpinfo$purity_pass_flag)
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

