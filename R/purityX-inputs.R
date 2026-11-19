# Prepare chromatographic memberships and retained MS1 observations together.
# Peak positions are references into the input's full chromPeaks table.
.px_prepare <- function(xset, singleFile = 0, fileignore = NULL,
    xgroups = NULL, rtrawColumns = FALSE) {
    if (!requireNamespace("xcms", quietly = TRUE))
        stop("purityX() requires the xcms package.")
    modern <- .xcms_is_modern(xset)
    if (!modern && !is(xset, "xcmsSet"))
        stop("xset must contain chromatographic peaks in an xcmsSet, ",
            "XCMSnExp or XcmsExperiment.")
    peaks <- as.data.frame(if (modern) xcms::chromPeaks(xset)
        else xset@peaks)
    required <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", "sample")
    if (!all(required %in% names(peaks)))
        stop("Chromatographic peaks are missing required mass or RT columns.")
    nsamples <- if (is(xset, "XcmsExperiment"))
        nrow(MsExperiment::sampleData(xset)) else length(.xcms_files(xset))
    .px_indices(singleFile, nsamples, "singleFile", zero = TRUE)
    if (length(singleFile) != 1L)
        stop("singleFile must be one sample index, or zero for grouped peaks.")
    if (!is.null(fileignore))
        .px_indices(fileignore, nsamples, "fileignore")
    if (length(rtrawColumns) != 1L || is.na(rtrawColumns) ||
            !is.logical(rtrawColumns))
        stop("rtrawColumns must be TRUE or FALSE.")
    if (rtrawColumns && !all(c("rt_raw", "rtmin_raw", "rtmax_raw") %in%
                names(peaks)))
        stop("rtrawColumns requires rt_raw, rtmin_raw and rtmax_raw columns.")
    .px_indices(peaks$sample, nsamples, "chromatographic peak samples")
    peaks$c_peak_id <- seq_len(nrow(peaks))
    peaks$peak_id <- rownames(peaks)
    filled <- if (modern) {
        .f4f_is_filled(xset, nrow(peaks))
    } else {
        if (!"sn" %in% names(peaks))
            stop("Legacy xcmsSet peaks must contain the sn column.")
        is.na(peaks$sn)
    }
    groups <- data.frame(grpid = integer(), feature_id = character(),
        i = numeric(), mz = numeric())
    intensity_method <- "not_used"
    if (singleFile > 0) {
        targets <- peaks[peaks$sample == singleFile, , drop = FALSE]
        targets$id <- seq_len(nrow(targets))
    } else {
        fd <- if (modern) xcms::featureDefinitions(xset) else xset@groups
        members <- if (modern) as.list(fd$peakidx) else xset@groupidx
        if (length(members) != nrow(fd))
            stop("Feature memberships do not match the feature table.")
        if (is.null(xgroups)) xgroups <- seq_len(nrow(fd))
        .px_indices(xgroups, nrow(fd), "xgroups")
        selected <- seq_len(nrow(fd))[seq_len(nrow(fd)) %in% xgroups]
        peak_idx <- unlist(members[selected], use.names = FALSE)
        if (is.null(peak_idx)) peak_idx <- integer()
        .px_indices(peak_idx, nrow(peaks), "feature peak memberships")
        peak_idx <- as.integer(peak_idx)
        targets <- peaks[peak_idx, , drop = FALSE]
        targets$grpid <- rep(selected, lengths(members[selected]))
        if (length(selected)) {
            values <- if (modern)
                xcms::featureValues(xset, value = "into", method = "medret")
            else xcms::groupval(xset, value = "into", method = "medret")
            intensity_method <- if (modern) "featureValues" else "groupval"
            values <- values[selected, , drop = FALSE]
            if (is(xset, "XcmsExperiment") && ncol(values) != nsamples) {
                # Shared origins can collapse featureValues() columns. Apply
                # its medret rule to each logical sample, in membership order.
                values <- matrix(NA_real_, nrow = length(selected),
                    ncol = nsamples)
                for (i in seq_along(selected)) {
                    feature <- selected[i]
                    idx <- members[[feature]]
                    idx <- idx[order(abs(peaks$rt[idx] - fd$rtmed[feature]))]
                    first <- match(seq_len(nsamples), peaks$sample[idx])
                    values[i, ] <- peaks$into[idx[first]]
                }
                intensity_method <- "sample_link_medret"
            }
            feature_ids <- rownames(fd)
            if (is.null(feature_ids))
                feature_ids <- rep(NA_character_, nrow(fd))
            groups <- data.frame(grpid = selected,
                feature_id = feature_ids[selected],
                i = apply(values, 1L,
                    stats::median, na.rm = TRUE),
                mz = as.numeric(fd[selected, "mzmed"]))
        }
    }
    reason <- rep(NA_character_, nrow(targets))
    reason[filled[targets$c_peak_id]] <- "filled_peak"
    if (modern) {
        level <- xcms::chromPeakData(xset)$ms_level
        if (!is.null(level))
            reason[which(level[targets$c_peak_id] != 1L)] <- "non_ms1_peak"
    }
    if (!singleFile && length(fileignore))
        reason[targets$sample %in% fileignore] <- "ignored_sample"
    exclusions <- targets[!is.na(reason), , drop = FALSE]
    exclusions$reason <- reason[!is.na(reason)]
    targets <- targets[is.na(reason), , drop = FALSE]
    rownames(targets) <- NULL
    sample_ids <- unique(as.integer(targets$sample))
    observations <- .px_observations(xset, sample_ids, nsamples)
    tracking <- observations$tracking
    selection <- lapply(seq_len(nrow(targets)), function(i) {
        t <- targets[i, , drop = FALSE]
        d <- tracking[tracking$sample == t$sample, , drop = FALSE]
        lo <- if (rtrawColumns) t$rtmin_raw else t$rtmin
        hi <- if (rtrawColumns) t$rtmax_raw else t$rtmax
        rt <- if (rtrawColumns) d$retentionTime else d$rt_corrected
        if (!is.finite(lo) || !is.finite(hi) || lo > hi)
            return(list(positions = integer(), status = "invalid_rt_bounds"))
        positions <- d$scan[which(d$msLevel == 1L & is.finite(rt) &
                    rt >= lo & rt <= hi)]
        list(positions = positions, status = if (length(positions))
                "available" else "no_ms1_observations")
    })
    targets$scan_positions <- lapply(selection, `[[`, "positions")
    targets$observation_status <- vapply(selection, `[[`, "", "status")
    targets$minscan <- vapply(targets$scan_positions, function(i) {
        if (length(i)) min(i) else NA_integer_
    }, integer(1))
    targets$maxscan <- vapply(targets$scan_positions, function(i) {
        if (length(i)) max(i) else NA_integer_
    }, integer(1))
    info <- list(class = class(xset), source_files = observations$files,
        sample_names = observations$sample_names,
        feature_intensity_method = intensity_method,
        rt_convention = if (rtrawColumns) "explicit_raw_peak_bounds"
        else if (observations$adjusted) "adjusted_peak_and_spectrum_times"
        else "input_peak_and_spectrum_times",
        peak_bounds = targets[, intersect(c("c_peak_id", "peak_id", "sample",
                    "grpid", "id", required, "rt_raw", "rtmin_raw",
                    "rtmax_raw"),
                names(targets)), drop = FALSE],
        memberships = targets[, intersect(c("grpid", "id", "sample",
                    "c_peak_id", "peak_id"), names(targets)), drop = FALSE])
    list(targets = targets, scanpeaks = observations$scanpeaks,
        tracking = tracking, groups = groups, inputInfo = info,
        exclusions = exclusions)
}

.px_indices <- function(x, n, name, zero = FALSE) {
    if (!is.numeric(x) || anyNA(x) || any(!is.finite(x)) ||
            any(x != trunc(x)) || any(x < if (zero) 0L else 1L) ||
            any(x > n))
        stop(name, " must contain valid integer positions.")
    invisible(x)
}

.px_tracking <- function(sp, sample, source_position = seq_along(sp),
    adjusted = Spectra::rtime(sp)) {
    d <- Spectra::spectraData(sp)
    data.frame(sample = rep.int(as.integer(sample), nrow(d)),
        scan = seq_len(nrow(d)), sourcePosition = as.integer(source_position),
        msLevel = as.integer(d$msLevel),
        retentionTime = as.numeric(d$rtime),
        rt_corrected = as.numeric(adjusted),
        acquisitionNum = if ("acquisitionNum" %in% names(d))
            as.integer(d$acquisitionNum) else rep(NA_integer_, nrow(d)),
        spectrumId = if ("spectrumId" %in% names(d))
            as.character(d$spectrumId) else rep(NA_character_, nrow(d)),
        dataOrigin = as.character(Spectra::dataOrigin(sp)))
}

.px_observations <- function(x, samples, n) {
    scanpeaks <- tracking <- vector("list", n)
    files <- rep(NA_character_, n)
    sample_names <- rep(NA_character_, n)
    adjusted <- FALSE
    if (is(x, "XcmsExperiment")) {
        sp <- MsExperiment::spectra(x)
        links <- MsExperiment::spectraSampleIndex(x, duplicates = "keep")
        positions <- rep(seq_along(links), lengths(links))
        sample_link <- unlist(links, use.names = FALSE)
        adjusted <- xcms::hasAdjustedRtime(x)
        rt <- xcms::rtime(x)
        sample_names <- rownames(MsExperiment::sampleData(x))
        for (i in samples) {
            idx <- positions[sample_link == i]
            s <- sp[idx]
            tracking[[i]] <- .px_tracking(s, i, idx, rt[idx])
            arrays <- vector("list", length(s))
            ms1 <- which(Spectra::msLevel(s) == 1L)
            arrays[ms1] <- .msp_peaks(s[ms1])
            scanpeaks[[i]] <- arrays
            origin <- unique(Spectra::dataOrigin(s))
            if (length(origin) == 1L) files[i] <- origin
        }
    } else if (is(x, "XCMSnExp")) {
        files <- MSnbase::fileNames(x)
        sample_names <- MSnbase::sampleNames(x)
        sample <- MSnbase::fromFile(x)
        raw_rt <- xcms::rtime(x, adjusted = FALSE)
        current_rt <- xcms::rtime(x)
        adjusted <- xcms::hasAdjustedRtime(x)
        metadata <- MSnbase::fData(x)
        # Extract observations without subsetting the chromatographic results.
        raw_input <- as(x, "OnDiskMSnExp")
        for (i in samples) {
            idx <- which(sample == i)
            ms1 <- which(MSnbase::msLevel(x)[idx] == 1L)
            measured <- MSnbase::spectra(raw_input[idx[ms1]],
                BPPARAM = .msp_bpparam())
            arrays <- vector("list", length(idx))
            arrays[ms1] <- lapply(measured, function(s) {
                cbind(mz = MSnbase::mz(s), intensity = MSnbase::intensity(s))
            })
            scanpeaks[[i]] <- arrays
            tracking[[i]] <- data.frame(sample = rep.int(i, length(idx)),
                scan = seq_along(idx), sourcePosition = as.integer(idx),
                msLevel = as.integer(MSnbase::msLevel(x)[idx]),
                retentionTime = as.numeric(raw_rt[idx]),
                rt_corrected = as.numeric(current_rt[idx]),
                acquisitionNum = as.integer(MSnbase::acquisitionNum(x)[idx]),
                spectrumId = if ("spectrumId" %in% names(metadata))
                    as.character(metadata$spectrumId[idx])
                else rep(NA_character_, length(idx)),
                dataOrigin = rep(files[i], length(idx)))
        }
    } else {
        files <- x@filepaths
        sample_names <- xcms::sampnames(x)
        for (i in samples) {
            sp <- .msp_read(files[i])
            d <- .px_tracking(sp, i)
            raw <- x@rt$raw[[i]]
            corrected <- x@rt$corrected[[i]]
            if (length(raw) != length(corrected))
                stop("Raw and corrected RTs differ in length for sample ",
                    i, ".")
            if (length(raw) && !identical(raw, corrected)) {
                d$rt_corrected <- .f4f_interp(d$retentionTime, raw, corrected)
                adjusted <- adjusted || !identical(raw, corrected)
            }
            tracking[[i]] <- d
            arrays <- vector("list", length(sp))
            ms1 <- which(Spectra::msLevel(sp) == 1L)
            arrays[ms1] <- .msp_peaks(sp[ms1])
            scanpeaks[[i]] <- arrays
        }
    }
    present <- Filter(Negate(is.null), tracking)
    tracking <- if (length(present)) do.call(rbind, present) else
        data.frame(sample = integer(), scan = integer(),
            sourcePosition = integer(), msLevel = integer(),
            retentionTime = numeric(), rt_corrected = numeric(),
            acquisitionNum = integer(), spectrumId = character(),
            dataOrigin = character())
    rownames(tracking) <- NULL
    list(scanpeaks = scanpeaks, tracking = tracking, files = files,
        sample_names = sample_names, adjusted = adjusted)
}
