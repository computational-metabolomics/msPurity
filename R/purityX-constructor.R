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



#' @import plyr
#' @import mzR
#' @import foreach
#' @import parallel
#' @import doSNOW
#' @import stringr
#' @import reshape2
#' @import fastcluster

#' @title Assessing anticipated purity of XCMS features from an LC-MS run
#'
#' @description
#' Constructor for the purityX class.
#'
#' Given an XCMS object get the anticipated precursor purity of the grouped peaks
#'
#' @param xset An \code{XcmsExperiment}, \code{XCMSnExp} or legacy
#'   \code{xcmsSet} containing detected chromatographic peaks. Grouped
#'   prediction also requires feature grouping.
#' @param cores numeric; number of cores to use
#' @param purityType character; Area and average used for the purity predictions. Options are
#'                               "purityFWHMmedian", "purityFWmedian", "purityFWHMmean", "purityFWmean"
#' @param offsets Isolation-window offsets in Da, in historical upper, lower
#'   order. For example, \code{c(0.7, 0.3)} assesses
#'   \code{[target - 0.3, target + 0.7]}. This order differs from purityA.
#' @param fileignore vector; vector of files to ignore for the prediction calculation
#' @param xgroups vector; vector of xcms groups to perform prediction on
#' @param iwNorm boolean; if TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5% (0.05)
#' @param plotP boolean; TRUE if plot of the EIC of feature and associated contamination is the be save to the working directory
#' @param mzRback character; deprecated and ignored. Raw data is read through Spectra, which uses the pwiz reader
#' @param isotopes boolean; TRUE if isotopes are to be removed
#' @param im matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#' @param rtrawColumns Use explicit \code{rt_raw}, \code{rtmin_raw} and
#'   \code{rtmax_raw} peak bounds when supplied. Otherwise compare current
#'   chromatographic-peak and spectrum retention times.
#' @param singleFile numeric; If just a single file for purity is to be calculated (rather than the grouped XCMS peaks). Uses the index of the files in xcmsSet object. If zero this is ignored.
#' @param saveEIC Save chromatograms for grouped predictions to SQLite. Ignored
#'   when \code{singleFile > 0}. SQLite writing is deprecated. Parquet output
#'   stores no chromatograms; extract them with \code{xcms::chromatogram()}.
#' @param sqlitePth character; If saveEIC True, then a path to sqlite database can be used. If NULL then a database will be created in the working directory called eics
#' @param param One complete \code{\link{PurityXParam}} object. Explicit
#'   scientific arguments supplied alongside it must agree. Sample/feature
#'   selection, execution and output options remain on this call.
#'
#' @details Predictions use the selected MS1 observations, including their
#'   surrounding isolation-window signal. Modern inputs retain Spectra subsets
#'   and pending peak processing. Removing interfering peaks before prediction
#'   can change the result. Full-width predictions summarize the chromatographic
#'   peak interval; FWHM predictions use its nearest half-height bounds. Missing
#'   flanks or absent signal leave FWHM predictions unavailable.
#'   Acquisition retention times are retained where available. If the caller
#'   has already applied adjusted times, original acquisition times cannot be
#'   recovered; the legacy EIC column \code{rt_raw} then contains the available
#'   input spectrum times.
#'   Prediction estimates anticipated isolation purity, not fragmentation quality.
#'
#' @return a purityX object containing a dataframe of predicted purity scores
#' @examples
#'
#' msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
#' xset <- readRDS(system.file("extdata", "tests", "xcms", "ms_only_xset_OLD.rds", package="msPurity"))
#' xset@filepaths[1] <- msPths[basename(msPths)=="LCMS_1.mzML"]
#' xset@filepaths[2] <- msPths[basename(msPths)=="LCMS_2.mzML"]
#'
#' ## Keep the example short: predict only the peaks between m/z 100 and 200
#' xset@peaks <- xset@peaks[xset@peaks[, "mz"] >= 100 & xset@peaks[, "mz"] <= 200, ]
#' px <- purityX(xset, singleFile = 1)
#'
#' @export
purityX <- function(xset, purityType = "purityFWHMmedian",
    offsets = c(0.5, 0.5),
    fileignore = NULL, cores = 1, xgroups = NULL,
    iwNorm = FALSE, iwNormFun = NULL, ilim = 0.05, plotP = FALSE,
    mzRback = "pwiz", isotopes = FALSE, im = NULL,
    singleFile = 0, rtrawColumns = FALSE, saveEIC = FALSE,
    sqlitePth = NULL, param = NULL) {
    scalars <- list(purityType = purityType, offsets = offsets,
        iwNorm = iwNorm, iwNormFun = iwNormFun, ilim = ilim,
        isotopes = isotopes, im = im)
    explicit <- intersect(names(match.call()), names(scalars))
    settings <- .px_resolve_parameters(param, scalars, explicit)
    if (!is.numeric(cores) || length(cores) != 1L || !is.finite(cores) ||
            cores < 1 || cores > .Machine$integer.max ||
            cores != as.integer(cores))
        stop("cores must be one positive integer.")
    for (name in c("plotP", "rtrawColumns", "saveEIC")) {
        value <- get(name)
        if (!is.logical(value) || length(value) != 1L || is.na(value))
            stop(name, " must be TRUE or FALSE.")
    }
    .msp_deprecate_mzRback(mzRback)
    if (isTRUE(saveEIC) && singleFile <= 0)
        .msp_deprecate_sqlite("purityX(saveEIC = TRUE)",
            paste("Parquet output stores no chromatograms; extract EICs",
                "when needed with xcms::chromatogram()."))
    if (settings$iwNorm && is.null(settings$iwNormFun)) {
        half_width <- sum(settings$offsets) / 2
        settings$iwNormFun <- iwNormGauss(minOff = -half_width,
            maxOff = half_width)
    }
    prepared <- .px_prepare(xset, singleFile, fileignore, xgroups,
        rtrawColumns)
    chosen <- .px_choose_targets(prepared$targets, singleFile > 0)
    targets <- chosen$targets
    workers <- if (!nrow(targets)) 0L else if (singleFile > 0) 1L else
        min(as.integer(cores), nrow(targets))
    tasks <- split(targets, seq_len(nrow(targets)))
    if (workers > 1L) {
        cl <- parallel::makeCluster(workers)
        on.exit(parallel::stopCluster(cl), add = TRUE)
        predictions <- parallel::parLapply(cl, tasks, .px_predict_one,
            prepared = prepared, settings = settings)
    } else {
        predictions <- lapply(tasks, .px_predict_one,
            prepared = prepared, settings = settings)
    }
    if (plotP && length(predictions)) {
        dir.create("purityXplots", showWarnings = FALSE)
        for (prediction in predictions) .px_plot(prediction)
    }
    if (saveEIC && singleFile <= 0 && length(predictions))
        .px_write_eics(predictions, sqlitePth)

    object <- new("purityX")
    object@predictions <- .px_prediction_table(predictions, prepared$groups,
        settings$purityType, singleFile > 0)
    if (singleFile <= 0) {
        object@purityType <- settings$purityType
        object@cores <- cores
        object@fileignore <- if (is.null(fileignore)) NA else fileignore
    }
    object
}

.px_choose_targets <- function(targets, single) {
    empty <- data.frame(output_id = numeric(), sample = numeric(),
        c_peak_id = numeric(), reason = character())
    if (!nrow(targets)) return(list(targets = targets, excluded = empty))
    if (single) {
        targets$output_id <- targets$id
        return(list(targets = targets, excluded = empty))
    }
    targets$output_id <- targets$grpid
    chosen <- integer()
    for (group in unique(targets$grpid)) {
        rows <- which(targets$grpid == group)
        middle <- median(targets$rt[rows])
        for (sample in unique(targets$sample[rows])) {
            members <- rows[targets$sample[rows] == sample]
            # Preserve the first occurrence when RT distances tie.
            distances <- abs(targets$rt[members] - middle)
            closest <- if (all(is.na(distances))) 1L else which.min(distances)
            chosen <- c(chosen, members[closest])
        }
    }
    excluded <- targets[setdiff(seq_len(nrow(targets)), chosen),
        c("output_id", "sample", "c_peak_id"), drop = FALSE]
    excluded$reason <- rep("alternate_peak", nrow(excluded))
    list(targets = targets[chosen, , drop = FALSE], excluded = excluded)
}

.px_predict_one <- function(target, prepared, settings) {
    metric_names <- c("purityFWHMmedian", "purityFWmedian",
        "pknmFWHMmedian", "pknmFWmedian", "purityFWHMmean", "purityFWmean")
    metrics <- stats::setNames(rep(NA_real_, length(metric_names)),
        metric_names)
    result <- list(target = target, metrics = metrics, trace = data.frame(),
        observations = prepared$tracking[FALSE, , drop = FALSE],
        fwhm = c(NA_real_, NA_real_), status = "invalid_rt_bounds",
        fwhm_status = "unavailable")
    if ("observation_status" %in% names(target) &&
            target$observation_status != "available") {
        result$status <- target$observation_status
        return(result)
    }
    bounds <- unlist(target[c("minscan", "maxscan")], use.names = FALSE)
    if (length(bounds) != 2L || any(!is.finite(bounds)) ||
            bounds[1] > bounds[2])
        return(result)
    tracking <- prepared$tracking
    sample_rows <- which(tracking$sample == target$sample)
    if ("scan_positions" %in% names(target)) {
        use <- sample_rows[match(target$scan_positions[[1L]],
                tracking$scan[sample_rows])]
        if (anyNA(use))
            stop("Prepared purityX observations have missing identities.")
    } else {
        use <- sample_rows[tracking$msLevel[sample_rows] == 1L &
                tracking$scan[sample_rows] >= bounds[1] &
                tracking$scan[sample_rows] <= bounds[2]]
    }
    use <- use[order(tracking$rt_corrected[use], tracking$scan[use])]
    observed <- tracking[use, , drop = FALSE]
    result$observations <- observed
    if (!length(use)) {
        result$status <- "no_ms1_observations"
        return(result)
    }
    arrays <- prepared$scanpeaks[[target$sample]][observed$scan]
    lower <- target$mz - settings$offsets[2]
    upper <- target$mz + settings$offsets[1]
    scores <- vapply(arrays, pcalc, numeric(2L), mzmin = lower,
        mzmax = upper, mztarget = target$mz, iwNorm = settings$iwNorm,
        iwNormFun = settings$iwNormFun, ilim = settings$ilim,
        targetMinMZ = target$mzmin, targetMaxMZ = target$mzmax,
        isotopes = settings$isotopes, im = settings$im)
    intensity <- vapply(arrays, getEic, numeric(1L), target = target)
    window_signal <- vapply(arrays, function(peaks) {
        if (is.null(peaks) || !nrow(peaks)) return(c(0, 0))
        isolated <- peaks[peaks[, 1] > lower & peaks[, 1] < upper,
            , drop = FALSE]
        c(sum(isolated[, 2]), getEic(isolated, target))
    }, numeric(2L))
    trace <- data.frame(scan = observed$scan, purity = scores[1, ],
        pknm = scores[2, ], intensity = intensity,
        sample = observed$sample, msLevel = observed$msLevel,
        acquisitionNum = observed$acquisitionNum,
        retentionTime = observed$retentionTime,
        rt_corrected = observed$rt_corrected,
        isolationIntensity = window_signal[1, ],
        isolatedTargetIntensity = window_signal[2, ])
    positions <- seq_len(nrow(trace))
    fwhm_rows <- calculateFWHM(data.frame(scan = positions,
            intensity = trace$intensity))
    within <- if (all(is.finite(fwhm_rows)))
        positions >= fwhm_rows[1] & positions <= fwhm_rows[2] else
        rep(FALSE, nrow(trace))
    fwhm <- if (all(is.finite(fwhm_rows))) trace$scan[fwhm_rows] else
        c(NA_real_, NA_real_)
    summarise <- function(x, fun) {
        x <- x[!is.na(x)]
        if (length(x)) fun(x) else NA_real_
    }
    metrics[] <- c(summarise(trace$purity[within], median),
        summarise(trace$purity, median),
        summarise(trace$pknm[within], median), summarise(trace$pknm, median),
        summarise(trace$purity[within], mean), summarise(trace$purity, mean))
    result$metrics <- metrics
    result$trace <- trace
    result$fwhm <- fwhm
    result$status <- if (any(!is.na(trace$purity))) "available" else
        "no_usable_purity"
    result$fwhm_status <- if (all(is.finite(fwhm))) "available" else
        "missing_flank_or_signal"
    result
}

.px_prediction_table <- function(predictions, groups, purityType, single) {
    single_fields <- c("purityFWHMmedian", "purityFWmedian",
        "pknmFWHMmedian", "pknmFWmedian")
    if (grepl("mean$", purityType))
        single_fields <- c(single_fields, purityType)
    if (single) {
        out <- as.data.frame(matrix(numeric(), nrow = 0L,
                ncol = length(single_fields) + 1L))
        names(out) <- c("id", single_fields)
        if (!length(predictions)) return(out)
        out <- do.call(rbind, lapply(predictions, function(x) {
            as.data.frame(as.list(c(id = as.numeric(x$target$id),
                        x$metrics[single_fields])))
        }))
        rownames(out) <- NULL
        return(out)
    }
    out <- data.frame(grpid = numeric(), mean = numeric(), median = numeric(),
        sd = numeric(), stde = numeric(), RSD = numeric(), pknm = numeric(),
        i = numeric(), mz = numeric())
    if (!nrow(groups)) return(out)
    group_ids <- unique(groups$grpid)
    prediction_groups <- vapply(predictions, function(x) x$target$grpid,
        numeric(1L))
    members_by_group <- split(predictions, factor(prediction_groups,
            levels = group_ids))
    group_index <- match(groups$grpid, group_ids)
    rows <- lapply(seq_len(nrow(groups)), function(i) {
        group <- groups[i, , drop = FALSE]
        members <- members_by_group[[group_index[i]]]
        values <- vapply(members, function(x) x$metrics[[purityType]],
            numeric(1L))
        values <- values[!is.na(values)]
        stat <- if (length(values)) c(mean(values), median(values), sd(values),
            stde(values), CV(values)) else rep(NA_real_, 5L)
        count_name <- if (startsWith(purityType, "purityFWHM"))
            "pknmFWHMmedian" else "pknmFWmedian"
        counts <- vapply(members, function(x) x$metrics[[count_name]],
            numeric(1L))
        count <- if (any(!is.na(counts))) median(counts, na.rm = TRUE) else
            NA_real_
        as.data.frame(as.list(stats::setNames(c(as.numeric(group$grpid), stat,
                        count, group$i, group$mz), names(out))))
    })
    do.call(rbind, rows)
}


.px_plot <- function(prediction) {
    trace <- prediction$trace
    if (!nrow(trace)) return(invisible(NULL))
    target <- prediction$target
    name <- paste(target$output_id, "sample", target$sample, "mz",
        round(target$mz, 3), "rt", round(target$rt, 3), sep = "_")
    contamination <- trace$isolationIntensity - trace$isolatedTargetIntensity
    upper <- max(c(trace$intensity, contamination), na.rm = TRUE)
    if (!is.finite(upper) || upper <= 0) upper <- 1
    grDevices::png(file.path("purityXplots", paste0(name, ".png")))
    on.exit(grDevices::dev.off())
    graphics::plot(trace$scan, trace$intensity, type = "l", col = "red",
        lwd = 2.5, xlab = "scan number", ylab = "intensity", ylim = c(0, upper))
    graphics::lines(trace$scan, contamination, col = "black", lwd = 2.5)
    graphics::legend("topright", lwd = c(2.5, 2.5),
        legend = c("EIC of feature", "contamination"), col = c("red", "black"))
    if (all(is.finite(prediction$fwhm))) graphics::abline(v = prediction$fwhm)
    invisible(NULL)
}

.px_write_eics <- function(predictions, sqlitePth) {
    path <- if (is.null(sqlitePth)) "eics.sqlite" else sqlitePth
    con <- DBI::dbConnect(RSQLite::SQLite(), path)
    on.exit(DBI::dbDisconnect(con))
    for (prediction in predictions) {
        dfp <- prediction$trace
        if (!nrow(dfp)) next
        dfp$isolationIntensity <- NULL
        dfp$isolatedTargetIntensity <- NULL
        dfp$sample <- as.character(dfp$sample)
        names(dfp)[names(dfp) == "retentionTime"] <- "rt_raw"
        # Retain the existing SQLite table's column order.
        dfp <- dfp[, c("scan", "purity", "pknm", "intensity", "sample",
                "msLevel", "acquisitionNum", "rt_raw", "rt_corrected")]
        dfp$grpid <- as.numeric(prediction$target$grpid)
        dfp$c_peak_id <- as.numeric(prediction$target$c_peak_id)
        dfp$eicidi <- seq_len(nrow(dfp))
        DBI::dbWriteTable(con, "eics", dfp, row.names = FALSE, append = TRUE)
    }
    invisible(NULL)
}

calculateFWHM <- function(df) {
    if (nrow(df) < 3L || any(!is.finite(df$intensity)) ||
            max(df$intensity) <= 0) return(c(NA_real_, NA_real_))
    apex <- which.max(df$intensity)
    left <- which(df$scan < df$scan[apex])
    right <- which(df$scan > df$scan[apex])
    if (!length(left) || !length(right)) return(c(NA_real_, NA_real_))
    half <- df$intensity[apex] / 2
    c(df$scan[left[which.min(abs(df$intensity[left] - half))]],
        df$scan[right[which.min(abs(df$intensity[right] - half))]])
}

getEic <- function(roi_scn, target) {
    if (is.null(roi_scn) || !nrow(roi_scn)) return(0)
    within <- which(roi_scn[, 1] >= target$mzmin & roi_scn[, 1] <= target$mzmax)
    if (!length(within)) return(0)
    chosen <- within[which.min(abs(roi_scn[within, 1] - target$mz))]
    unname(roi_scn[chosen, 2])
}

stde <- function(x) sd(x) / sqrt(length(x))
CV <- function(x) (sd(x) / mean(x)) * 100
