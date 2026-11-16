# Diagnostic purity uses the original parent measurement, not filtered peaks
# or the MS2 ancestor's purity. It never supplies a filtering decision.
.pa_msn_purity <- function(source, ancestry, pids, ppm, tolerance, offsets) {
    scans <- ancestry$scans
    rows <- which(scans$ms_level > 2L & scans$pid %in% pids)
    parent <- match(scans$parent_id[rows], scans$spectrum_id)
    column <- function(name) {
        if (name %in% names(scans)) scans[[name]][rows] else
            rep(NA_real_, length(rows))
    }
    result <- data.frame(pid = scans$pid[rows],
        parent_pid = scans$pid[parent], parent_id = scans$spectrum_id[parent],
        precursor_mz = scans$precursor_mz[rows],
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
        result$upper_mz, column("isolationWindowTargetMz"), offsets)
    result$lower_mz <- windows$isolationWindowLowerMz
    result$upper_mz <- windows$isolationWindowUpperMz
    result$window_source <- windows$source
    valid_parent <- !is.na(parent) & scans$parent_status[rows] == "valid"
    valid_parent[valid_parent] <- scans$ms_level[parent[valid_parent]] ==
        scans$ms_level[rows[valid_parent]] - 1L
    result$status[valid_parent] <- "invalid_isolation_window"
    valid <- valid_parent & is.finite(result$lower_mz) &
        is.finite(result$upper_mz) & result$lower_mz >= 0 &
        result$lower_mz < result$precursor_mz &
        result$upper_mz > result$precursor_mz
    valid[is.na(valid)] <- FALSE
    parents <- unique(parent[valid])
    original <- Spectra::reset(source[parents])
    arrays <- if (length(parents)) .msp_peaks(original) else list()
    centroided <- Spectra::centroided(original)
    for (i in which(valid)) {
        k <- match(parent[i], parents)
        peaks <- as.matrix(arrays[[k]])
        result$status[i] <- "profile_parent"
        if (identical(centroided[k], FALSE)) next
        result$status[i] <- "empty_parent"
        if (!nrow(peaks)) next
        result$status[i] <- "invalid_parent_peaks"
        if (any(!is.finite(peaks)) || any(peaks[, "mz"] <= 0) ||
                any(peaks[, "intensity"] < 0)) next
        inside <- peaks[, "mz"] > result$lower_mz[i] &
            peaks[, "mz"] < result$upper_mz[i]
        distance <- abs(peaks[, "mz"] - result$precursor_mz[i])
        eligible <- which(inside & distance <=
                abs(result$precursor_mz[i]) * ppm / 1e6 + tolerance)
        result$status[i] <- "missing_precursor_peak"
        if (!length(eligible)) next
        selected <- eligible[distance[eligible] == min(distance[eligible])]
        result$status[i] <- "ambiguous_precursor_peak"
        if (length(selected) != 1L) next
        result$matched_mz[i] <- peaks[selected, "mz"]
        result$parent_peak_index[i] <- selected
        result$status[i] <- "zero_precursor_intensity"
        if (peaks[selected, "intensity"] == 0) next
        estimate <- pcalc(peaks, mzmin = result$lower_mz[i],
            mzmax = result$upper_mz[i], mztarget = result$matched_mz[i])
        result$status[i] <- "invalid_estimate"
        if (!is.finite(estimate[1L]) || estimate[1L] < 0 || estimate[1L] > 1)
            next
        result$purity[i] <- unname(estimate[1L])
        result$peaks_in_window[i] <- as.integer(estimate[2L])
        result$status[i] <- "assessed"
    }
    # Keep MS2 estimates; unassociated higher levels are not assessed.
    purity <- if ("inPurity" %in% Spectra::spectraVariables(source))
        source$inPurity else rep(NA_real_, length(source))
    higher <- Spectra::msLevel(source) > 2L
    purity[higher] <- NA_real_
    status <- ifelse(is.na(purity), "unavailable", "assessed")
    status[higher] <- "not_assessed"
    purity[rows] <- result$purity
    status[rows] <- result$status
    source$msPurityStagePurity <- purity
    source$msPurityPurityStatus <- status
    parameters <- list(method = "parent_spectrum",
        ppm = ppm, tolerance = tolerance, ilim = 0, isotopes = FALSE,
        iwNorm = FALSE, interpolation = "none", observations = result)
    list(spectra = source, parameters = parameters)
}
