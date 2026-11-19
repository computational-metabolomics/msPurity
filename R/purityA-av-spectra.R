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



#' @title Using a purityA object, average and filter fragmentation spectra for each XCMS feature within a MS data file
#' @aliases averageIntraFragSpectra
#' @description
#'
#' **General**
#'
#' Average and filter fragmentation spectra for each XCMS feature within a MS data file.
#'
#' The averaging is performed using hierarchical clustering of the m/z values of each peaks, where m/z values within a set ppm tolerance will be clustered. The clustered peaks are then averaged (or summed).
#'
#' The fragmentation can be filtered on the averaged spectra (with the arguments snr, rsd, minfrac and ra)
#'
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * XCMS processing (version >= 3)
#'    +  (mzML files) -> MSnBase.readMSdata -> xcms.findChromPeaks -> xcms.groupChromPeaks-> xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> frag4feature -> filterFragSpectra -> **averageIntraFragSpectra** -> averageIntraFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#'
#' **Averaging with parameter objects**
#'
#' With \code{param}, each parameter set averages the feature-linked MS2, MS3
#' and MS4 spectra it selects within each file, keeping incompatible
#' acquisition conditions and precursor paths apart, as described for
#' \code{\link{averageAllFragSpectra}} and
#' \code{\link{AverageFragSpectraParam}}. Spectra of different files are
#' never averaged together; \code{msPurityAverageGroup} is numbered within
#' each feature and file. The result's \code{sample} identifies the acquisition
#' file in \code{pa@@fileList}; contributor records retain the XCMS sample
#' index where available. Grouped feature associations do not need an
#' individual chromatographic peak. The call replaces earlier intra- and
#' inter-file averages and keeps all-spectrum averages. Without \code{param}, the
#' function replaces earlier intra- and inter-file averages too, and stops
#' with an error while all-spectrum averages computed with \code{param} are
#' present.
#'
#' @inherit groupPeaksEx details
#'
#' @param pa object; purityA object
#' @param ppm numeric; ppm threshold to average within each file
#' @param minnum numeric; minimum number of times peak is present across fragmentation spectra within each file
#' @param minfrac numeric; minimum ratio of the peak fraction (peak count / total peaks) within each file
#' @param ra numeric; minimum relative abundance of the peak within each file
#' @param snr numeric; minimum signal-to-noise of the peak within each file
#' @param av character; type of averaging to use (median or mean)
#' @param sumi boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#' @param param NULL, one \code{\link{AverageFragSpectraParam}} object or a
#'   named list of them, averaged independently.
#' @param cores Maximum number of workers; defaults to \code{pa@@cores}.
#'   Overrides apply to the current call. With \code{param}, workers process
#'   independent clustering blocks; a single block runs directly.
#' @inheritParams groupPeaksEx
#'
#' @return Returns a purityA object (pa) with the following slots now with data
#'
#' * pa@@avSpectra: the averaged spectra as a \code{Spectra} object, one spectrum per grouped feature, averaging
#'   level and file, with the averaging statistics as peak variables.
#' * pa@@av_spectra: the legacy list of averaged spectra, filled while \code{options(msPurity.legacySlots = TRUE)}.
#'   e.g. "averagedSpectra(pa)$`1`$av_intra$`1`" gives the average spectra for grouped feature 1 and for file 1.
#' * pa@@av_intra_params: The parameters used are recorded here
#'
#' Each averaged spectrum in the legacy list contains the following columns (peak variables in pa@@avSpectra):
#'
#' * cl: id of clustered (averaged) peak
#' * mz: average m/z
#' * i: average intensity
#' * snr: average signal to noise ratio
#' * rsd: relative standard deviation
#' * count: number of clustered peaks
#' * total: total number of potential scans to be used for averaging
#' * inPurity: average precursor ion purity
#' * ra: average relative abundance
#' * frac: the fraction of clustered peaks (e.g. the count/total)
#' * snr_pass_flag: TRUE if snr threshold criteria met
#' * minfrac_pass_flag: TRUE if minfrac threshold criteria
#' * ra_pass_flag: TRUE if ra threshold criteria met
#' * pass_flag: TRUE if all threshold criteria met
#'
#' @examples
#'
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
#' #pa <- averageIntraFragSpectra(pa)
#'
#' # Run from previously generated data (where class is 'XCMSnExp'):
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'               "2_frag4feature_pa.rds", package="msPurity"))
#' pa <- averageIntraFragSpectra(pa)
#' @md
#' @export
setMethod(f="averageIntraFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                               av='median', sumi=TRUE, rmp=FALSE, cores=1,
                               blockSize=5000L, param=NULL
                                ){
            if (missing(cores)) cores <- max(1L, pa@cores)
            .clustering_controls(cores, blockSize)
            pa <- .pa_update(pa)
            if (!is.null(param)) {
              supplied <- list(minfrac = minfrac, minnum = minnum, ppm = ppm,
                snr = snr, ra = ra, av = av, sumi = sumi)[c(
                  !missing(minfrac), !missing(minnum), !missing(ppm),
                  !missing(snr), !missing(ra), !missing(av), !missing(sumi))]
              return(.pa_average_configured(pa, "intra", param, supplied, rmp,
                cores, blockSize))
            }
            pa <- .pa_average_scalar_ready(pa, "intra")

            pa@av_intra_params$minfrac = minfrac
            pa@av_intra_params$minnum = minnum
            pa@av_intra_params$ppm = ppm
            pa@av_intra_params$snr = snr
            pa@av_intra_params$ra = ra

            pa@av_intra_params$av_type = av
            pa@av_intra_params$sumi = sumi


            pa@av_intra_params$cores = cores
            pa@av_intra_params$blockSize = blockSize
            pa@av_intra_params$rmp = rmp

            return(average_xcms_grouped_msms(pa, "intra", cores, blockSize))

          }
)


#' @title Using a purityA object, average and filter fragmentation spectra for each XCMS feature across multiple MS data files
#' @aliases averageInterFragSpectra
#' @description
#'
#' **General**
#'
#' Average and filter fragmentation spectra for each XCMS feature across MS data files. This can only be run after averageIntraFragSpectra has been used.
#'
#' The averaging is performed using hierarchical clustering of the m/z values of each peaks, where m/z values within a set ppm tolerance will be clustered. The clustered peaks are then averaged (or summed).
#'
#' The fragmentation can be filtered on the averaged spectra (with the arguments snr, rsd, minfrac and ra)
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
#'    + (xcmsObj, pa) -> frag4feature -> filterFragSpectra -> averageIntraFragSpectra -> **averageInterFragSpectra** -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#'
#'
#' **Averaging with parameter objects**
#'
#' With \code{param}, inter-file averaging combines intra-file averages
#' computed with \code{param}. Each parameter set must have an intra-file
#' parameter set of the same name (a single object is named
#' \code{"default"}) with the same \code{select} and \code{combine}; its
#' thresholds and tolerances may differ. Intra-file averages of the same
#' feature are combined when their conditions match as for
#' \code{\link{averageAllFragSpectra}} and every MS3/MS4 precursor of the
#' contributors to passing peaks agrees within the inter-file tolerance. The
#' call stops if these precursors are incompatible within an intra-file input.
#' Each file contributes at most one input, and each input counts once.
#'
#' Only intra-file averages with passing peaks are inputs. Compatibility and
#' numerical averaging use those same peaks, whether failed intra-file peaks
#' were flagged or removed. Contributors are traced through the intra-file
#' clusters: \code{msPurityContributors} lists links contributing to retained
#' final peaks. Use the result's \code{spectrumId} with
#' \code{\link{averageMembers}} for its immediate intra-file inputs or with
#' \code{\link{averagePeakMembers}} to trace a final fragment through
#' intermediate peaks to the original scans. Without
#' \code{param},
#' the function stops with an error while intra-file or all-spectrum averages
#' computed with \code{param} are present.
#'
#' @inherit groupPeaksEx details
#'
#' @param pa object; purityA object
#' @param cores Maximum number of workers; defaults to \code{pa@@cores}.
#'   Overrides apply to the current call. With \code{param}, workers process
#'   independent clustering blocks; a single block runs directly.
#' @inheritParams groupPeaksEx
#' @param ppm numeric; ppm threshold to average across files
#' @param minnum numeric; minimum number of times peak is present across fragmentation spectra across files
#' @param minfrac numeric; minimum ratio of the peak fraction (peak count / total peaks) across files
#' @param ra numeric; minimum relative abundance of the peak across files
#' @param snr numeric; minimum signal-to-noise of the peak across files
#'
#' @param av character; type of averaging to use (median or mean)
#' @param sumi boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#' @param param NULL, one \code{\link{AverageFragSpectraParam}} object or a
#'   named list of them, averaged independently.
#'
#' @return Returns a purityA object (pa) with the following slots now with data
#'
#' * pa@@avSpectra: the averaged spectra as a \code{Spectra} object, one spectrum per grouped feature, averaging
#'   level and file, with the averaging statistics as peak variables.
#' * pa@@av_spectra: the legacy list of averaged spectra, filled while \code{options(msPurity.legacySlots = TRUE)}.
#'   e.g. "averagedSpectra(pa)$`1`$av_inter" gives the average spectra for grouped feature 1
#' * pa@@av_intra_params: The parameters used are recorded here
#'
#' Each averaged spectrum in the legacy list contains the following columns (peak variables in pa@@avSpectra):
#' *
#' * cl: id of clustered (averaged) peak
#' * mz: average m/z
#' * i: average intensity
#' * snr: average signal to noise ratio
#' * rsd: relative standard deviation
#' * count: number of clustered peaks
#' * total: total number of potential scans to be used for averaging
#' * inPurity: average precursor ion purity
#' * ra: average relative abundance
#' * frac: the fraction of clustered peaks (e.g. the count/total)
#' * snr_pass_flag: TRUE if snr threshold criteria met
#' * minfrac_pass_flag: TRUE if minfrac threshold criteria
#' * ra_pass_flag: TRUE if ra threshold criteria met
#' * pass_flag: TRUE if all threshold criteria met
#'
#' @examples
#'
#' #====== XCMS =====================================
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
#' #====== msPurity =================================
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xcmsObj)
#' #pa <- averageIntraFragSpectra(pa)
#' #pa <- averageInterFragSpectra(pa)
#'
#' ## Run from previously generated data
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "4_averageIntraFragSpectra_no_filter_pa.rds",
#'                           package="msPurity"))
#' pa <- averageInterFragSpectra(pa)
#' @md
#' @export
setMethod(f="averageInterFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                                av='median', sumi=TRUE,  rmp=FALSE, cores=1,
                               blockSize=5000L, param=NULL
          ){
            if (missing(cores)) cores <- max(1L, pa@cores)
            .clustering_controls(cores, blockSize)
            pa <- .pa_update(pa)
            if (!is.null(param)) {
              supplied <- list(minfrac = minfrac, minnum = minnum, ppm = ppm,
                snr = snr, ra = ra, av = av, sumi = sumi)[c(
                  !missing(minfrac), !missing(minnum), !missing(ppm),
                  !missing(snr), !missing(ra), !missing(av), !missing(sumi))]
              return(.pa_average_configured(pa, "inter", param, supplied, rmp,
                cores, blockSize))
            }
            pa <- .pa_average_scalar_ready(pa, "inter")

            pa@av_inter_params$minfrac = minfrac
            pa@av_inter_params$minnum = minnum
            pa@av_inter_params$ppm = ppm
            pa@av_inter_params$snr = snr
            pa@av_inter_params$ra = ra

            pa@av_inter_params$av_type = av
            pa@av_inter_params$sumi = sumi


            pa@av_inter_params$cores = cores
            pa@av_inter_params$blockSize = blockSize
            pa@av_inter_params$rmp = rmp

            first_grp <- unique(as.character(pa@grped_df$grpid))[1]
            if (is.null(.pa_av_legacy(pa@avSpectra, pa@grped_df)[[first_grp]][["av_intra"]])){
              stop("Apply averageIntraFragSpectra first")
            }

            return(average_xcms_grouped_msms(pa, "inter", cores, blockSize))

          }
)


#' @title Using a purityA object, average and filter MS/MS spectra for each XCMS feature within
#' and across MS data files (ignoring intra and inter relationships)
#' @aliases averageAllFragSpectra
#' @description
#'
#' **General**
#'
#' Average and filter fragmentation spectra for each XCMS feature within and across MS data files (ignoring intra and inter relationships).
#'
#' The averaging is performed using hierarchical clustering of the m/z values of each peaks, where m/z values within a set ppm tolerance will be clustered. The clustered peaks are then averaged (or summed).
#' With \code{param}, the tolerance can also include an absolute part (see
#' \code{\link{AverageFragSpectraParam}}).
#'
#' The fragmentation can be filtered on the averaged spectra (with the arguments snr, rsd, minfrac, ra)
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> (optionally) xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'    +  --- *Older versions of XCMS* --- (mzML files) -> xcms.xcmsSet -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> frag4feature -> filterFragSpectra -> **averageAllFragSpectra** -> createDatabase -> spectralMatching -> (sqlite spectral database)
#'
#' **Averaging with parameter objects**
#'
#' With \code{param}, each parameter set averages the feature-linked MS2, MS3
#' and MS4 spectra it selects and keeps incompatible acquisition conditions and
#' precursor paths apart (see \code{\link{AverageFragSpectraParam}}), so a
#' feature can have several averaged spectra. The settings come from the
#' parameter objects: an explicitly supplied \code{minfrac}, \code{minnum},
#' \code{ppm}, \code{snr}, \code{ra}, \code{av} or \code{sumi} must agree
#' with them. \code{rmp} and \code{cores} apply to the whole call.
#'
#' These averages are available from \code{averagedSpectra(pa, legacy = FALSE)}.
#' Each records its parameter set (\code{msPurityParameterSet}; a single
#' object is named \code{"default"}) and its group within the feature
#' (\code{msPurityAverageGroup}, \code{"<set>/<number>"}). With
#' \code{rmp = FALSE}, output peaks include a \code{pass_flag} column:
#' \code{1} for peaks passing the averaging thresholds and \code{0} for
#' failed peaks. Retrieve it with \code{Spectra::peaksData(s,
#' columns = c("mz", "intensity", "pass_flag"))}, where \code{s} is the
#' averaged Spectra object. With \code{rmp = TRUE}, failed peaks are removed.
#' Empty peak arrays indicate that no output peaks are stored.
#' Each averaged spectrum has an object-local character \code{spectrumId}. Use it with
#' \code{\link{averageMembers}} to inspect selected inputs and
#' \code{\link{averagePeakMembers}} to trace stored peaks. Measured spectra
#' keep their existing spectrum identifiers. Use \code{msPurityProcessingId}
#' to retrieve the result's processing record with
#' \code{\link{processingHistory}}. In the returned Spectra,
#' \code{msPurityContributors} lists original occurrences supporting stored
#' peaks, and \code{msPurityClusterContributors} maps clusters to those
#' occurrences. The
#' acquisition summaries, precursor m/z, retention time and
#' \code{msPurityPrecursorPath} (MS3 and MS4 precursor m/z) describe the
#' contributing spectra; \code{msPurityAncestorConditions} gives the
#' acquisition conditions of the earlier stages of MS3 and MS4 averages.
#' A feature linked to the same spectrum more than once contributes it once
#' per link, as without \code{param}.
#'
#' The legacy list returned by \code{averagedSpectra(pa)}, \code{createMSP()}
#' with averaged spectra and \code{createDatabase()} hold one average per
#' feature, level and file; they stop with an error while averages computed
#' with \code{param} are present. For MS2-only feature links, averaging again
#' without \code{param} replaces them; averaging without \code{param} does not
#' support MS3/MS4 feature links. Averaging without \code{param} stops with an
#' error while averages computed with \code{param} remain at another level.
#' \code{averageIntraFragSpectra()} and \code{averageInterFragSpectra()} accept
#' \code{param} in the same way.
#'
#' @inherit groupPeaksEx details
#'
#' @param pa object; purityA object
#' @param cores Maximum number of workers; defaults to \code{pa@@cores}.
#'   Overrides apply to the current call. With \code{param}, workers process
#'   independent clustering blocks; a single block runs directly.
#' @inheritParams groupPeaksEx
#' @param ppm numeric; ppm threshold to average across all scans (ignoring intra and inter relationships)
#' @param minnum numeric; minimum number of times peak is present across all fragmentation spectra (ignoring intra and inter relationships)
#' @param minfrac numeric;minimum ratio of the peak fraction (peak count / total peaks) across all (ignoring intra and inter relationships)
#' @param ra numeric; minimum relative abundance of the peak fraction across all (ignoring intra and inter relationships)
#' @param snr numeric;  minimum signal-to-noise of the peak across all (ignoring intra and inter relationships)
#' @param av character; type of averaging to use (median or mean)
#' @param sumi boolean; TRUE if the intensity for each peak is summed across averaged spectra
#' @param rmp boolean; TRUE if peaks are to be removed that do not meet the threshold criteria. Otherwise they will just be flagged
#' @param param NULL, one \code{\link{AverageFragSpectraParam}} object or a
#'   named list of them, averaged independently.
#'
#' @return Returns a purityA object (pa) with the following slots now with data
#'
#' * pa@@avSpectra: the averaged spectra as a \code{Spectra} object, one spectrum per grouped feature, averaging
#'   level and file, with the averaging statistics as peak variables. With \code{param}, a feature can have
#'   several all-spectrum averages, one per parameter set and group.
#' * pa@@av_spectra: the legacy list of averaged spectra, filled while \code{options(msPurity.legacySlots = TRUE)}.
#'   E.g. averagedSpectra(pa)$`1`$av_all gives the average spectra for grouped feature 1. Averages computed with
#'   \code{param} are stored only in pa@@avSpectra.
#' * pa@@av_all_params: The parameters used are recorded here
#'
#' Each averaged spectrum in the legacy list contains the following columns (peak variables in pa@@avSpectra):
#'
#' * cl: id of clustered (averaged) peak
#' * mz: average m/z
#' * i: average intensity
#' * snr: average signal to noise ratio
#' * rsd: relative standard deviation
#' * count: number of clustered peaks
#' * total: total number of potential scans to be used for averaging
#' * inPurity: average precursor ion purity
#' * ra: average relative abundance
#' * frac: the fraction of clustered peaks (e.g. the count/total)
#' * snr_pass_flag: TRUE if snr threshold criteria met
#' * minfrac_pass_flag: TRUE if minfrac threshold criteria
#' * ra_pass_flag: TRUE if ra threshold criteria met
#' * pass_flag: TRUE if all threshold criteria met
#'
#'
#' @examples
#'
#' #====== XCMS ===================================
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
#' #====== msPurity ================================
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xcmsObj)
#' #pa <- filterFragSpectra(pa)
#' #pa <- averageAllFragSpectra(pa)
#'
#' ## Run from previously generated data
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "3_filterFragSpectra_pa.rds", package="msPurity"))
#' pa <- averageAllFragSpectra(pa)
#'
#' ## Keep activation types apart while combining collision energies
#' pa_param <- averageAllFragSpectra(pa,
#'     param = AverageFragSpectraParam(combine = "collisionEnergy"))
#' sp <- averagedSpectra(pa_param, legacy = FALSE)
#' i <- 1L
#' averageMembers(pa_param, sp$spectrumId[i])
#' sp$msPurityContributors[[i]]
#' sp$msPurityClusterContributors[[i]]
#'
#' @md
#' @export
setMethod(f="averageAllFragSpectra", signature="purityA",
          definition = function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                                 av='median', sumi=TRUE, rmp=FALSE, cores=1,
                               blockSize=5000L, param=NULL
          ){
            if (missing(cores)) cores <- max(1L, pa@cores)
            .clustering_controls(cores, blockSize)
            pa <- .pa_update(pa)
            if (!is.null(param)) {
              supplied <- list(minfrac = minfrac, minnum = minnum, ppm = ppm,
                snr = snr, ra = ra, av = av, sumi = sumi)[c(
                  !missing(minfrac), !missing(minnum), !missing(ppm),
                  !missing(snr), !missing(ra), !missing(av), !missing(sumi))]
              return(.pa_average_configured(pa, "all", param, supplied, rmp,
                cores, blockSize))
            }
            pa <- .pa_average_scalar_ready(pa, "all")

            pa@av_all_params$minfrac = minfrac
            pa@av_all_params$minnum = minnum
            pa@av_all_params$ppm = ppm
            pa@av_all_params$snr = snr
            pa@av_all_params$ra = ra

            pa@av_all_params$av_type = av
            pa@av_all_params$sumi = sumi
            pa@av_all_params$cores = cores
            pa@av_all_params$blockSize = blockSize
            pa@av_all_params$rmp = rmp

            return(average_xcms_grouped_msms(pa, "all", cores, blockSize))

          }
)



average_xcms_grouped_msms <- function(pa, av_level, cores = pa@cores,
    blockSize = 5000L){

    if (nrow(pa@grped_df)) {
        rows <- match(unique(pa@grped_df$pid), pa@spectra$pid)
        if (any(Spectra::msLevel(pa@spectra)[rows] > 2L, na.rm = TRUE))
            stop("Averaging MS3/MS4 feature links is not supported.")
    }

    processing_started <- Sys.time()
    processing_step <- switch(av_level,
        intra = "averageIntraFragSpectra", inter = "averageInterFragSpectra",
        all = "averageAllFragSpectra")
    if (av_level == "inter") pa <- .pa_graph_import(pa)
    processing_input <- .pa_processing_input(pa, processing_step)
    previous <- pa@avSpectra
    intra <- if (av_level == "inter")
        .pa_graph_view(pa, previous[previous$av_level == "av_intra"])
    if (av_level == "intra")
        pa <- .pa_processing_invalidate(pa, "inter")


  # Inputs come from the Spectra slots. Levels not being averaged now are
  # carried over unchanged, from the legacy slot when it is filled.
  ms2 <- if (av_level != "inter") .pa_grouped_legacy(pa, peak_identity = TRUE)
  av_in <- .pa_av_legacy(pa@avSpectra, pa@grped_df)
  av_keep <- if (.pa_legacy() && length(pa@av_spectra)) pa@av_spectra else av_in
  groups <- unique(as.character(pa@grped_df$grpid))
  workers <- min(cores, max(1L, length(groups)))
  para <- workers > 1L
  if (para) {
    cl <- parallel::makeCluster(workers)
    on.exit({
      parallel::stopCluster(cl)
      foreach::registerDoSEQ()
    }, add = TRUE)
    doSNOW::registerDoSNOW(cl)
  }


    results <- plyr::alply(groups, 1, average_xcms_grouped_msms_indiv,
        pa = pa, av_level = av_level, ms2 = ms2, av_in = av_in,
        av_keep = av_keep, record_inputs = TRUE, blockSize = blockSize,
        .parallel = para)
    names(results) <- groups
    used_inputs <- lapply(results, `[[`, "used")
    av_spectra <- results
    for (i in seq_along(results))
        av_spectra[i] <- list(results[[i]]$spectra)


  pa@avSpectra <- .pa_av_to_spectra(av_spectra, pa@grped_df)
    pa <- .pa_av_metadata(pa, av_level, previous, used_inputs, intra)
    sp <- pa@avSpectra
    rows <- if (length(sp)) which(sp$av_level == paste0("av_", av_level)) else integer()
    graph_groups <- lapply(rows, function(i) {
        graph <- results[[as.character(sp$grpid[i])]]$graph
        if (av_level != "intra") return(graph)
        if (is.na(sp$sample[i]))
            return(list(links = integer(), edges = .pa_average_peak_edges(NULL, NULL)))
        graph[[as.character(sp$sample[i])]]
    })
    pa <- .pa_graph_store(pa, rows, graph_groups)
  pa <- .pa_sync_legacy(pa, av_spectra = av_spectra)

    effective_parameters <- methods::slot(
        pa, paste0("av_", av_level, "_params"))
    names(effective_parameters)[
        names(effective_parameters) == "av_type"] <- "av"
    effective_parameters$effective_cores <- as.integer(workers)
    pa <- .pa_processing_record(pa, processing_step, effective_parameters,
        processing_started, processing_input)
  return(pa)

}



average_xcms_grouped_msms_indiv <- function(grp_idx, pa, av_level, ms2, av_in,
    av_keep, record_inputs = FALSE, blockSize = 5000L) {

  grp_idx_char <- as.character(grp_idx)
    used_inputs <- list()

  ##############################################################################
  # Get the appropiate details for the xcms grouped feature from purityA object
  ##############################################################################
  graph <- list()
  if (av_level != "inter") {
    link_rows <- which(as.character(pa@grped_df$grpid) == grp_idx_char)
    grped_info <- pa@grped_df[link_rows, , drop = FALSE]
    grped_info$index <- seq_along(link_rows)
    spectra_to_average <- .pa_average_peak_table(grped_info, ms2[[grp_idx_char]],
        filtered = length(pa@filter_frag_params) > 0)
    spectra_to_average <- spectra_to_average[order(spectra_to_average$mz), ]
  }

  # Set return variable to empty list or already existing results
  if (!is.null(av_keep[[grp_idx_char]][["av_intra"]])){
    av_intra = av_keep[[grp_idx_char]][["av_intra"]]
  } else {
    av_intra = NULL
  }

  # Set return variable to empty list or already existing results
  if (!is.null(av_keep[[grp_idx_char]][["av_inter"]])){
    av_inter = av_keep[[grp_idx_char]][["av_inter"]]
  } else {
    av_inter = NULL
  }

  # Set return variable to empty list or already existing results
  if (!is.null(av_keep[[grp_idx_char]][["av_all"]])){
    av_all = av_keep[[grp_idx_char]][["av_all"]]
  } else {
    av_all = NULL
  }

  # filter out peaks below precursor ion purity thres
  if (!(av_level %in% c("intra", "inter", "all"))){
    stop("Incorrect av_level for averaging fragmentation spectra; use intra, inter or all")
  }


  ##############################################################################
  # Performing averaging
  ##############################################################################

  if (av_level=="intra"){
    # Average by sample (file)

    av_intra <- plyr::dlply(spectra_to_average, ~sample, average_spectra,
                                            cores=1,
                                            ppm=pa@av_intra_params$ppm,
                                            minnum=pa@av_intra_params$minnum,
                                            sumi=pa@av_intra_params$sumi,
                                            minfrac=pa@av_intra_params$minfrac,
                                            snthr=pa@av_intra_params$snr,
                                            rathr=pa@av_intra_params$ra,
                                            rathr_pre= pa@av_intra_params$ra_pre,
                                            snrthr_pre= pa@av_intra_params$snr_pre,
                                            av_type = pa@av_intra_params$av_type,
                                            record_inputs = record_inputs, blockSize = blockSize)

        if (record_inputs) {
            results <- av_intra
            members <- lapply(results, `[[`, "members")
            for (i in seq_along(av_intra))
                av_intra[i] <- list(av_intra[[i]]$peaks)
        }

    if (pa@av_intra_params$rmp){
      av_intra  <- plyr::llply(av_intra , function(x){x[x$pass_flag,]})
    }
        if (record_inputs) {
            used_inputs$intra <- Map(.pa_average_used_members, members, av_intra)
            samples <- if ("sample" %in% names(grped_info))
                grped_info$sample else grped_info$fileid
            graph <- lapply(names(results), function(sample) {
                ranks <- which(as.character(samples) == sample)
                list(links = link_rows[ranks], edges = .pa_average_peak_edges(
                    results[[sample]], av_intra[[sample]], ranks))
            })
            names(graph) <- names(results)
        }

  } else if (av_level=="inter") {

    av_intra_df <- plyr::ldply(av_in[[grp_idx_char]][["av_intra"]],
        .id = "sample", function(x) {
            if (is.null(x)) return(NULL)
            x$memberPeak <- seq_len(nrow(x))
            x$basis <- rep("derived", nrow(x))
            x[x$pass_flag, , drop = FALSE]
        })
    samples <- unique(as.character(av_intra_df$sample))
    av_intra_df$memberRank <- match(as.character(av_intra_df$sample), samples)

    # Average the averaged spectra across files
    av_inter <- average_spectra(av_intra_df,
                                 indx='sample',
                                 cores=1,
                                 ppm=pa@av_inter_params$ppm,
                                 minnum=pa@av_inter_params$minnum,
                                 sumi=pa@av_inter_params$sumi,
                                 minfrac=pa@av_inter_params$minfrac,
                                 snthr=pa@av_inter_params$snr,
                                 rathr=pa@av_inter_params$ra,
                                 av_type = pa@av_inter_params$av_type,
                                 record_inputs = record_inputs, blockSize = blockSize
                               )

        if (record_inputs) {
            worker <- av_inter
            members <- worker$members
            av_inter <- worker$peaks
        }

    if (pa@av_inter_params$rmp){
      av_inter <- av_inter[av_inter$pass_flag,]
    }
        if (record_inputs) {
            used_inputs$inter <- .pa_average_used_members(members, av_inter)
            intra <- pa@avSpectra[pa@avSpectra$av_level == "av_intra" &
                pa@avSpectra$grpid == grp_idx_char]
            inputs <- intra$spectrumId[match(samples, as.character(intra$sample))]
            if (anyNA(inputs)) stop("Inter-file input has no derived spectrumId.")
            graph <- list(inputs = inputs, edges = .pa_average_peak_edges(
                worker, av_inter, seq_along(inputs)))
        }

  } else if (av_level=="all") {

    # add additional column used later if filtering applied
    # Average the averaged spectra across everything (ignore intra and inter )
    av_all <- average_spectra(spectra_to_average,
                                cores=1,
                                ppm=pa@av_all_params$ppm,
                                minnum=pa@av_all_params$minnum,
                                sumi=pa@av_all_params$sumi,
                                minfrac=pa@av_all_params$minfrac,
                                snthr=pa@av_all_params$snr,
                                rathr=pa@av_all_params$ra,
                                rathr_pre= pa@av_all_params$ra_pre,
                                snrthr_pre= pa@av_all_params$snr_pre,
                                av_type = pa@av_all_params$av_type,
                                record_inputs = record_inputs, blockSize = blockSize)
        if (record_inputs) {
            worker <- av_all
            members <- worker$members
            av_all <- worker$peaks
        }
    if (pa@av_all_params$rmp){
      av_all <- av_all[av_all$pass_flag,]
    }
        if (record_inputs) {
            used_inputs$all <- .pa_average_used_members(members, av_all)
            graph <- list(links = link_rows, edges = .pa_average_peak_edges(
                worker, av_all, grped_info$index))
        }
  } else {

    stop("Incorrect av_level for averaging fragmentation spectra; use intra, inter or all")

  }

    spectra <- list(av_intra = av_intra, av_inter = av_inter, av_all = av_all)
    if (record_inputs)
        return(list(spectra = spectra, used = used_inputs, graph = graph))
    spectra
}


# Peak rows of a feature's linked spectra, as received by average_spectra().
# info holds the feature-link rows, with their position within the feature in
# index; spectra holds their peak matrices in the same order.
.pa_average_peak_table <- function(info, spectra, filtered) {
  spectra <- plyr::llply(spectra, data.frame)
  df <- data.frame(do.call("rbind", spectra))
  # With filterFragSpectra(rmp = TRUE) every spectrum of a feature can be NULL;
  # an empty table then takes the same path as a feature with no passing peaks
  if (ncol(df) == 0){
    df <- data.frame(mz = numeric(), i = numeric(), pass_flag = numeric())
  }

  colnames(df)[1:2] <- c('mz', 'i')
  df$index <- rep(info$index, vapply(spectra, nrow, integer(1)))
  if ("memberPeak" %in% names(df)) df$memberRank <- df$index

  if (filtered){
    # if prior filtering performed only use those that have passed
    df<-df[df$pass_flag==1,]
  }

  # Grouped feature associations have a source file but no individual peak.
  if (!"sample" %in% names(info)) info$sample <- info$fileid
  columns <- intersect(c('grpid', 'sample', 'cid', 'index', 'inPurity'),
    names(info))
  merge(df, info[, columns, drop = FALSE], by = "index")
}

# Occurrence membership follows the clusters retained by the existing worker.
.pa_average_used_members <- function(members, peaks) {
    if (is.null(peaks) || !nrow(peaks)) return(list())
    members[as.character(peaks$cl)]
}


# Resolve the worker's cluster labels to positions of the actual output peaks
# after rmp, and input occurrence indices to their ranks in the selected group.
.pa_average_peak_edges <- function(result, peaks, ranks = integer()) {
    if (is.null(peaks) || !nrow(peaks))
        return(data.frame(mergedPeak = integer(), rank = integer(),
            memberPeak = integer(), basis = character(), stringsAsFactors = FALSE))
    rows <- result$peak_members
    if (is.null(rows)) stop("Averaging did not record its input peak identities.")
    output <- match(rows$cl, peaks$cl)
    keep <- !is.na(output)
    rank <- match(rows$memberRank[keep], ranks)
    if (anyNA(rank)) stop("Averaged peaks refer to unselected input members.")
    data.frame(mergedPeak = output[keep], rank = rank,
        memberPeak = rows$memberPeak[keep], basis = rows$basis[keep],
        stringsAsFactors = FALSE)
}

average_spectra <- function(spectra, indx = "index", ppm, cores, minnum, sumi,
    minfrac, snthr, snmeth = "median", rathr, rathr_pre = NULL,
    snrthr_pre = NULL, av_type = "median", record_inputs = FALSE,
    tolerance = 0, blockSize = 5000L) {
  if (nrow(spectra)==0){
        if (record_inputs) return(list(peaks = NULL, members = list(),
            effective_cores = 0L))
    return(NULL)
  }


  if (indx=='index'){
    # calculate metrics per scan (if using inter, the index will be sample and the snr and ra will
    # have already have been calculated
    # these will have already been calculated if filterFragSpectra has already been applied
    if ((!'snr' %in% colnames(spectra)) & (!'ra' %in% colnames(spectra))){
      spectra <- plyr::ddply(spectra, indx, set_snr_ra)
    }

  }


  if (!is.null(rathr_pre)){
    spectra <- spectra[spectra$ra>rathr_pre, ]
  }

  if (!is.null(snrthr_pre)){
    spectra <- spectra[spectra$snr>snrthr_pre, ]
  }

  if (nrow(spectra)==0){
        if (record_inputs) return(list(peaks = NULL, members = list(),
            effective_cores = 0L))
    return(NULL)
  }

  # ensure ordered by mz
  spectra <- spectra[order(spectra$mz),]

  mz <- spectra$mz

    # Inter-file inputs refer to peaks in an earlier intra-file average.
    if (record_inputs) {
        input_keys <- if (indx == "sample")
            data.frame(sample = as.character(spectra[[indx]]),
                cluster = as.character(spectra$cl),
                stringsAsFactors = FALSE) else
            spectra[[indx]]
    }

  # Cluster the peaks togther
  clusters <- clustering(mz, clustType = 'hc', cores = cores, ppm = ppm,
    tolerance = tolerance, blockSize = blockSize, record_execution = record_inputs)
  spectra$cl <- if (record_inputs) clusters$clusters else clusters

  averages <- plyr::ddply(spectra, ~ cl,
                          averageCluster, av=av_type, minnum=1,
                          missingV="ignore", totalScans=length(unique(spectra[,indx])), normTIC=FALSE,
                          sumI=sumi)




  averages$frac <- averages$count/averages$total

  averages$snr_pass_flag <- averages$snr > snthr

  averages$minnum_pass_flag <- averages$count >= minnum

  averages$minfrac_pass_flag <- averages$frac >= minfrac

  averages$ra_pass_flag <- averages$ra > rathr

  averages$pass_flag <- (averages$minfrac_pass_flag & averages$snr_pass_flag & averages$ra_pass_flag & averages$minnum_pass_flag)

    if (record_inputs) {
        members <- lapply(
            split(seq_len(nrow(spectra)), spectra$cl), function(i) {
                if (is.data.frame(input_keys))
                    unique(input_keys[i, , drop = FALSE])
                else unique(input_keys[i])
            })
        identity <- c("cl", "memberRank", "memberPeak", "basis")
        peak_members <- if (all(identity %in% names(spectra)))
            spectra[, identity, drop = FALSE]
        return(list(peaks = averages, members = members, peak_members = peak_members,
            effective_cores = clusters$effective_cores))
    }
  return(averages)
}


set_snr_ra <- function(x, snmeth='median'){

  if (snmeth=="median"){
    x$snr <- x$i/median(x$i)
  }else if(snmeth=="mean"){
    x$snr <- x$i/mean(x$i)
  }
  x$ra <- x$i/max(x$i)*100

   return(x)
}

