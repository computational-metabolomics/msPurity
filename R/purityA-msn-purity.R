# Return positions in the supplied peak array to record the selected peak.
# MS2 retains its first-choice tie rule; higher levels reject ambiguous choices.
.pa_precursor_candidates <- function(peaks, precursor_mz, precursor_intensity,
    lower, upper, mostIntense = FALSE) {
    candidates <- which(peaks[, 1L] >= lower & peaks[, 1L] <= upper)
    if (!length(candidates) || !is.finite(precursor_mz)) return(integer())
    if (mostIntense) {
        intensity <- peaks[candidates, 2L]
        return(candidates[intensity == max(intensity)])
    }
    if (length(precursor_intensity) == 1L && is.finite(precursor_intensity)) {
        matched <- candidates[round(peaks[candidates, 2L], 1L) ==
                round(precursor_intensity, 1L)]
        if (length(matched)) candidates <- matched
    }
    distance <- abs(peaks[candidates, 1L] - precursor_mz)
    candidates[distance == min(distance)]
}

# Read the supplied parent measurement, including its queued peak processing.
.pa_msn_purity <- function(source, ancestry, pids, mostIntense = FALSE,
    offsets = NA_real_, ilim = 0, isotopes = FALSE, im = NULL,
    configuration = NULL) {
    scans <- ancestry$scans
    rows <- which(scans$ms_level > 2L & scans$pid %in% pids)
    parent <- match(scans$parentLineageScanId[rows], scans$lineageScanId)
    column <- function(name) {
        if (name %in% names(scans)) scans[[name]][rows] else
            rep(NA_real_, length(rows))
    }
    result <- data.frame(pid = scans$pid[rows],
        childLineageScanId = scans$lineageScanId[rows],
        parentPid = scans$pid[parent], parentLineageScanId = scans$lineageScanId[parent],
        precursor_mz = scans$precursor_mz[rows],
        precursor_intensity = column("precursor_intensity"),
        native_isolation_target_mz = column("isolationWindowTargetMz"),
        native_lower_mz = column("isolationWindowLowerMz"),
        native_upper_mz = column("isolationWindowUpperMz"),
        lower_mz = column("isolationWindowLowerMz"),
        upper_mz = column("isolationWindowUpperMz"),
        window_source = rep("recorded_metadata", length(rows)),
        matched_mz = rep(NA_real_, length(rows)),
        parent_peak_index = rep(NA_integer_, length(rows)),
        purity = rep(NA_real_, length(rows)),
        peaks_in_window = rep(NA_integer_, length(rows)),
        status = rep("missing_parent", length(rows)))
    if (is.null(offsets)) offsets <- NA_real_
    windows <- .pa_isolation_windows(result$precursor_mz, result$lower_mz,
        result$upper_mz, result$native_isolation_target_mz,
        if (is.null(configuration)) offsets else NA_real_)
    settings <- rep(list(list(mostIntense = mostIntense, ilim = ilim,
                isotopes = isotopes, im = im)), length(rows))
    if (!is.null(configuration)) {
        assignment <- configuration$assignments
        ids <- assignment$parameter_id[match(result$pid, assignment$pid)]
        settings <- unname(configuration$parameter_sets[ids])
        for (id in unique(ids)) {
            group <- which(ids == id)
            windows[group, ] <- .pa_isolation_windows(
                result$precursor_mz[group], result$lower_mz[group],
                result$upper_mz[group],
                result$native_isolation_target_mz[group],
                configuration$parameter_sets[[id]]$offsets)
        }
    }
    result$lower_mz <- windows$isolationWindowLowerMz
    result$upper_mz <- windows$isolationWindowUpperMz
    result$window_source <- windows$source
    valid_parent <- !is.na(parent) & scans$parentStatus[rows] == "valid"
    valid_parent[valid_parent] <- scans$ms_level[parent[valid_parent]] ==
        scans$ms_level[rows[valid_parent]] - 1L
    result$status[valid_parent] <- windows$status[valid_parent]
    valid <- valid_parent & windows$status == "valid"
    parents <- unique(parent[valid])
    measurements <- source[parents]
    arrays <- if (length(parents)) .msp_peaks(measurements) else list()
    centroided <- Spectra::centroided(measurements)
    for (i in which(valid)) {
        rule <- settings[[i]]
        k <- match(parent[i], parents)
        peaks <- as.matrix(arrays[[k]])
        result$status[i] <- "profile_parent"
        if (identical(centroided[k], FALSE)) next
        result$status[i] <- "empty_parent"
        if (!nrow(peaks)) next
        result$status[i] <- "invalid_parent_peaks"
        if (any(!is.finite(peaks)) || any(peaks[, "mz"] <= 0) ||
                any(peaks[, "intensity"] < 0)) next
        selected <- .pa_precursor_candidates(peaks, result$precursor_mz[i],
            result$precursor_intensity[i], result$lower_mz[i],
            result$upper_mz[i], rule$mostIntense)
        result$status[i] <- "missing_precursor_peak"
        if (!length(selected)) next
        result$status[i] <- "ambiguous_precursor_peak"
        if (length(selected) != 1L) next
        # pcalc identifies its numerator by m/z, so duplicated masses are
        # ambiguous even when recorded intensity distinguishes their rows.
        if (sum(peaks[, "mz"] == peaks[selected, "mz"]) != 1L) next
        result$matched_mz[i] <- peaks[selected, "mz"]
        result$parent_peak_index[i] <- selected
        result$status[i] <- "zero_precursor_intensity"
        if (peaks[selected, "intensity"] == 0) next
        estimate <- pcalc(peaks, mzmin = result$lower_mz[i],
            mzmax = result$upper_mz[i], mztarget = result$matched_mz[i],
            ilim = rule$ilim, isotopes = rule$isotopes, im = rule$im)
        result$status[i] <- "invalid_estimate"
        if (!is.finite(estimate[1L]) || estimate[1L] < 0 || estimate[1L] > 1)
            next
        result$purity[i] <- unname(estimate[1L])
        result$peaks_in_window[i] <- as.integer(estimate[2L])
        result$status[i] <- "assessed"
    }
    variables <- Spectra::spectraVariables(source)
    purity <- if ("inPurity" %in% variables) source$inPurity else
        rep(NA_real_, length(source))
    missing_status <- ifelse(Spectra::msLevel(source) > 2L,
        "not_assessed", "unavailable")
    status <- if ("msPurityPurityStatus" %in% variables)
        source$msPurityPurityStatus else
        ifelse(is.na(purity), missing_status, "assessed")
    purity[rows] <- result$purity
    status[rows] <- result$status
    source$inPurity <- purity
    source$msPurityStagePurity <- purity
    source$msPurityPurityStatus <- status
    parameters <- list(method = "parent_spectrum", mostIntense = mostIntense,
        offsets = offsets, ilim = ilim, isotopes = isotopes, im = im,
        iwNorm = FALSE, interpolation = "none", observations = result)
    if (!is.null(configuration))
        parameters <- parameters[c("method", "iwNorm", "interpolation",
                "observations")]
    list(spectra = source, parameters = parameters)
}
