.pa_automatic_offsets <- function(offsets) {
    automatic <- (is.numeric(offsets) || is.logical(offsets)) &&
        length(offsets) %in% 1:2 && all(is.na(offsets)) &&
        !any(is.nan(offsets))
    if (!automatic && (!is.numeric(offsets) || length(offsets) != 2L ||
                any(!is.finite(offsets) | offsets < 0) || sum(offsets) <= 0))
        stop("Explicit isolation offsets must be two nonnegative finite ",
            "numbers with positive total width.")
    automatic
}

# Resolve bounds for MS1-based purity and linked-parent diagnostics.
# Native Spectra bounds are absolute masses, not widths around precursorMz.
.pa_isolation_windows <- function(precursor, lower, upper, target,
    offsets = NA_real_) {
    automatic <- .pa_automatic_offsets(offsets)
    if (!automatic) {
        target <- precursor
        lower <- precursor - offsets[1L]
        upper <- precursor + offsets[2L]
    }
    # Absolute bounds suffice when the reported target is unavailable.
    unknown_target <- !is.finite(target)
    offset_reference <- target
    offset_reference[unknown_target] <-
        (lower[unknown_target] + upper[unknown_target]) / 2
    valid <- is.finite(lower) & is.finite(upper) & is.finite(precursor) &
        lower >= 0 & lower < upper & precursor >= lower & precursor <= upper &
        offset_reference >= lower & offset_reference <= upper
    valid[is.na(valid)] <- FALSE
    data.frame(precursorMz = precursor, isolationWindowTargetMz = target,
        isolationWindowLowerMz = lower, isolationWindowUpperMz = upper,
        lowerOffset = offset_reference - lower,
        upperOffset = upper - offset_reference,
        source = rep(if (automatic) "recorded_metadata" else "caller",
            length(precursor)),
        status = ifelse(valid, "valid", "invalid_isolation_window"))
}

.msp_isolation_settings <- function(sp) {
    .pa_isolation_windows(Spectra::precursorMz(sp),
        Spectra::isolationWindowLowerMz(sp),
        Spectra::isolationWindowUpperMz(sp),
        Spectra::isolationWindowTargetMz(sp))
}

# Keep a scalar summary for uniform windows; the table is authoritative when
# conditions differ within a file. Do not retain spectra or arrays in history.
.pa_isolation_history <- function(windows, iwNorm, iwNormFun) {
    offsets <- unique(windows[, c("lowerOffset", "upperOffset")])
    uniform <- nrow(offsets) == 1L && all(is.finite(as.matrix(offsets)))
    generated <- is.null(iwNormFun)
    functions <- if (iwNorm && generated) lapply(seq_len(nrow(windows)),
        function(i) {
            if (windows$status[i] != "valid") return(NULL)
            iwNormGauss(minOff = -windows$lowerOffset[i],
                maxOff = windows$upperOffset[i])
        }) else NULL
    list(offsets = if (uniform) unname(unlist(offsets)) else
            c(NA_real_, NA_real_),
        offsets_source = if (nrow(windows)) windows$source[1L] else "not_used",
        offsets_method = if (any(windows$source == "caller"))
            "explicit_offsets" else "spectra_per_spectrum",
        windows = windows, iwNorm = iwNorm,
        iwNormFun = if (!iwNorm) NULL else if (!generated) iwNormFun else
            if (uniform) functions[[1L]] else NULL,
        iwNormFunctions = if (iwNorm && generated && !uniform)
            functions else NULL,
        iwNormFun_source = if (!iwNorm) "not_used" else if (generated)
            "generated_default" else "caller")
}

.pa_purity_window <- function(row, windows, iwNormFun) {
    window <- windows[match(row$seqNum, windows$seqNum), , drop = FALSE]
    if (is.null(iwNormFun) && window$status == "valid")
        iwNormFun <- iwNormGauss(minOff = -window$lowerOffset,
            maxOff = window$upperOffset)
    list(valid = identical(window$status, "valid"),
        minoff = row$precursorMZ - window$isolationWindowLowerMz,
        maxoff = window$isolationWindowUpperMz - row$precursorMZ,
        iwNormFun = iwNormFun)
}
