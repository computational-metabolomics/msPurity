# Project live XCMS inputs into the converter's vocabulary. Retain the
# source peak matrix and sparse feature membership. Quantification of one
# feature uses its medret representative for each sample, so a dense
# features-by-samples matrix is never needed.
.mzs_live_source <- function(pa, x, metadata) {
    reserved <- length(pa@avSpectra)
    pa <- .pa_ms2_view(pa)
    modern <- .xcms_is_modern(x)
    if (!modern && !is(x, "xcmsSet"))
        .mzs_abort("semantic", "Expected an XcmsExperiment, XCMSnExp or xcmsSet.")
    paths <- pa@fileList
    xpaths <- .xcms_files(x)
    if (pa@f4f_link_type == "individual" &&
        !identical(unname(basename(paths)), unname(basename(xpaths))) &&
        !identical(names(paths), unname(basename(xpaths))))
        .mzs_abort("semantic", "createDatabase() could not match the files of 'pa' and 'xcmsObj'.")
    uneven <- length(paths) > length(xpaths) && pa@f4f_link_type == "group"
    file_index <- if (uneven) match(basename(xpaths), basename(paths)) else seq_along(xpaths)
    if (anyNA(file_index))
        .mzs_abort("semantic", "An XCMS file is absent from the purityA file list.")
    classes <- rep(NA_character_, length(paths))
    labels <- .xcms_classes(x)
    if (!is.null(labels)) classes[file_index] <- as.character(labels)
    files <- data.frame(fileid = seq_along(paths), filename = basename(paths),
        filepth = unname(paths), nm_save = names(paths) %||% rep(NA_character_, length(paths)),
        class = classes)
    peaks <- if (modern) xcms::chromPeaks(x) else x@peaks
    sample <- as.integer(peaks[, "sample"])
    fileid <- file_index[sample]
    peak_slice <- function(i) {
        cp <- as.data.frame(peaks[i, , drop = FALSE])
        cp$cid <- i
        cp$fileid <- fileid[i]
        cp$i <- cp$sample <- NULL
        names(cp)[names(cp) == "into"] <- "_into"
        names(cp) <- gsub("[.-]", "_", names(cp))
        cp
    }
    g <- if (modern) xcms::featureDefinitions(x) else as.data.frame(x@groups)
    members <- if (modern) as.list(g$peakidx) else x@groupidx
    if (modern) g$peakidx <- NULL
    g <- as.data.frame(g)
    g$grp_name <- if (modern) .xcms_feature_names(x) else xcms::groupnames(x)
    g$grpid <- seq_len(nrow(g))
    names(g)[match(c("mzmed", "rtmed"), names(g))] <- c("mz", "rt")
    names(g) <- gsub("[.-]", "_", make.names(names(g), unique = TRUE))
    n <- lengths(members)
    cid <- as.integer(unlist(members, use.names = FALSE))
    # Match getGroupPeakLink()'s historical per-sample index, including
    # its handling of inputs where peaks from two samples are interleaved.
    idi <- sequence(as.integer(table(sample)))
    links <- data.frame(cid = cid, grpid = rep(g$grpid, n), idi = idi[cid])
    best <- lapply(seq_along(members), function(j) {
        idx <- members[[j]]
        rt <- peaks[idx, "rt"]
        center <- if (modern) g$rt[j] else stats::median(rt)
        idx <- idx[order(abs(rt - center))]
        idx[!duplicated(sample[idx])]
    })
    links$bestpeak <- unlist(Map(`%in%`, members, best), use.names = FALSE)
    value <- function(group, assay) {
        out <- rep(NA_real_, length(group))
        for (which_group in split(seq_along(group), group)) {
            b <- best[[group[which_group[1L]]]]
            at <- match(assay[which_group], fileid[b])
            out[which_group] <- peaks[b[at], "into"]
        }
        out
    }
    classes <- setdiff(names(g), .MZS_GROUP_FIXED)
    for (col in classes) {
        v <- suppressWarnings(as.numeric(g[[col]]))
        if (any(!is.na(g[[col]]) & (is.na(v) | v != round(v) | v < 0 |
                                     v > g$npeaks), na.rm = TRUE))
            .mzs_abort("semantic", "Unrecognised XCMS feature column '", col, "' (R-098).")
    }
    scan <- .db_scan_metadata(purityTable(pa))
    scalar <- .mzs_scalar_index(pa, g, scan, reserved)
    if (!is.null(scalar)) scan <- plyr::rbind.fill(scan, scalar$meta)
    if (is.list(metadata) && !anyNA(metadata))
        for (c in intersect(c("polarity", "instrument", "instrument_type"), names(metadata)))
            scan[[c]] <- metadata[[c]]
    # Keep the SQLite-compatible lexical precision of old source fields.
    scan <- .mzs_sqlite_stored(list(df = scan, name = "s_peak_meta", pk = "pid",
        pk_type = "INTEGER", fks = list(fileid = list(new_name = "fileid"))))
    src <- list(fileinfo = files, c_peaks = data.frame(cid = seq_along(sample), fileid = fileid),
        c_peak_groups = g, c_peak_X_c_peak_group = links, s_peak_meta = scan,
        s_peaks = .pa_empty_scan_peaks(scan$pid))
    by_peak <- pa@f4f_link_type == "individual"
    key <- if (by_peak) "cid" else "grpid"
    link <- unique(pa@grped_df[, c("pid", key), drop = FALSE])
    link[[if (by_peak) "cXp_id" else "gXp_id"]] <- seq_len(nrow(link))
    src[[if (by_peak) "c_peak_X_s_peak_meta" else "c_peak_group_X_s_peak_meta"]] <- link
    src <- .mzs_type_src(src)
    src$scalar <- scalar
    src$peak_slice <- peak_slice
    src$peak_sn <- "sn" %in% colnames(peaks)
    src$has_peakidx <- modern
    src$abundance <- value
    src$sample_columns <- list(samples = seq_along(paths), classes = classes)
    src$lexical <- list(fileid = as.character(files$fileid), grpid = as.character(g$grpid),
        cid = as.character(src$c_peaks$cid), pid = as.character(scan$pid))
    if (!isTRUE(pa@filter_frag_params$allfrag))
        attr(src$s_peaks, "raw_scan_peaks") <- .pa_scan_peak_count(pa, .mzs_scans(src)$pid)
    src
}
