# msPurity results -> a Parquet results dataset.
#
# createDatabase(format = "parquet") stages a temporary SQLite database and
# enriches it from the live purityA object; convertSqliteToParquet() reads an
# existing one. Both map the `.mzs_read_sqlite()` table set:
#
#   fileinfo                      -> runs of the study (or external sources)
#   averaged spectra              -> native runs: av_intra_<run> per source
#                                    run, av_inter, av_all
#   their contributing scans      -> merge_member
#
# MS2 scans are referenced, never copied.

# Averaged-spectrum peak columns: s_peaks column -> (peak column, type).
# `frac` (count / total) and the composite `pass_flag` are not stored; `ra`
# cannot be recomputed from the averaged peak, so it is kept.
.MZS_AV_PEAKS <- list(
    snr = c("sn", "f64"),
    count = c("contributor_count", "i32"),
    total = c("contributor_total", "i32"),
    rsd = c("intensity_rsd", "f64"),
    cl = c("cluster", "i32"),
    ra = c("x_mspurity_ra", "f64"),
    inPurity = c("x_mspurity_in_purity", "f64"),
    snr_pass_flag = c("x_mspurity_snr_pass_flag", "bool"),
    minnum_pass_flag = c("x_mspurity_minnum_pass_flag", "bool"),
    minfrac_pass_flag = c("x_mspurity_minfrac_pass_flag", "bool"),
    ra_pass_flag = c("x_mspurity_ra_pass_flag", "bool"))

# Composite pass rules for scans and averaged spectra.
.MZS_PASS_CRITERIA <- list(
    scan = paste("purity_pass_flag & intensity_pass_flag & ra_pass_flag &",
                 "snr_pass_flag"),
    averaged = paste("minfrac_pass_flag & snr_pass_flag & ra_pass_flag &",
                     "minnum_pass_flag"))

.MZS_PRECURSOR_DEFINITION <- paste(
    "An averaged spectrum's precursor m/z and retention time are those of",
    "its feature (xcms featureDefinitions() mzmed and rtmed, or",
    "peakTable() mz and rt), as createDatabase() records them.")

# Metadata column types of an averaged spectrum.
.MZS_AV_META_TYPES <- c(ms_level = "i32", time = "f64", scan_polarity = "i32",
                        spectrum_representation = "str", id = "str",
                        data_origin = "str", selected_ion_mz = "f64",
                        x_mspurity_precursor_purity = "f64")

# ---------------------------------------------------------------------------
# Source files and runs.
# ---------------------------------------------------------------------------

#' The source files, each matched to a run.
#'
#' Duplicate file names are refused unless `fileMap` covers every file. With
#' a study, each file maps to its study run by path or `fileMap`; otherwise
#' each becomes an external source with a run id minted from its checksum
#' (or path, if unreadable) and ordinal.
#'
#' @return data.frame: fileid, filename, filepth, class, run_id, source,
#'     data_origin, readable, checksum; attributes `sources` (entries) and
#'     `study` (path).
#'
#' @noRd
.mzs_files <- function(fi, study, studyKey, fileMap, dest) {
    fi <- fi[order(fi$fileid), , drop = FALSE]
    base <- basename(fi$filepth)
    dup <- unique(base[duplicated(base)])
    covered <- length(fileMap) && all(fi$filepth %in% names(fileMap))
    if (length(dup) && !covered)
        .mzs_abort("semantic", "Source files share a file name, which ",
                   "msPurity keys its own tables on: ",
                   paste0(fi$filepth[base %in% dup], collapse = ", "),
                   ". Pass 'fileMap' naming the run of every file (R-096).",
                   data = list(paths = fi$filepth[base %in% dup]))
    readable <- !is.na(fi$filepth) & file.exists(fi$filepth)
    checksum <- rep(NA_character_, nrow(fi))
    checksum[readable] <- vapply(fi$filepth[readable], .mzs_sha256,
                                 character(1))
    mapped <- function(i) {
        p <- fi$filepth[i]
        if (!is.null(fileMap) && p %in% names(fileMap))
            unname(fileMap[[p]]) else NA_character_
    }

    if (!is.null(study)) {
        sm <- .mzs_read_manifest(study)
        meta <- .mzs_read_spectra_meta(study, sm, "data_origin")
        origin <- unique(meta[, c("run_id", "data_origin")])
        norm <- function(p) normalizePath(p, mustWork = FALSE)
        run <- vapply(seq_len(nrow(fi)), function(i) {
            r <- mapped(i)
            if (!is.na(r))
                return(r)
            hit <- unique(origin$run_id[!is.na(origin$data_origin) &
                                        norm(origin$data_origin) ==
                                        norm(fi$filepth[i])])
            if (length(hit) != 1L)
                .mzs_abort("semantic", "Source file '", fi$filepth[i],
                           "' matches ", length(hit), " runs of the study ",
                           "'", study, "'. Pass 'fileMap' naming its run.",
                           data = list(path = fi$filepth[i]))
            hit
        }, character(1))
        miss <- setdiff(run, .mzs_runs_frame(sm)$run_id)
        if (length(miss))
            .mzs_abort("semantic", "No run(s) ", paste(miss, collapse = ", "),
                       " in the study '", study, "'.")
        files <- data.frame(
            fileid = fi$fileid, filename = fi$filename,
            filepth = fi$filepth, class = fi$class, run_id = run,
            source = studyKey,
            data_origin = origin$data_origin[match(run, origin$run_id)],
            readable = readable, checksum = checksum,
            stringsAsFactors = FALSE)
        attr(files, "sources") <- list(
            .mzs_dataset_source(studyKey, study, dest, "study", unique(run)))
        attr(files, "study") <- study
        attr(files, "study_manifest") <- sm
    } else {
        run <- vapply(seq_len(nrow(fi)), function(i) {
            r <- mapped(i)
            if (!is.na(r))
                return(r)
            .mzs_mint_run_id(if (readable[i]) checksum[i]
                             else normalizePath(fi$filepth[i],
                                                mustWork = FALSE), i)
        }, character(1))
        if (anyDuplicated(run))
            .mzs_abort("semantic", "Two source files map to one run id: ",
                       paste(unique(run[duplicated(run)]), collapse = ", "),
                       ".")
        files <- data.frame(
            fileid = fi$fileid, filename = fi$filename,
            filepth = fi$filepth, class = fi$class, run_id = run,
            source = run, data_origin = fi$filepth, readable = readable,
            checksum = checksum, stringsAsFactors = FALSE)
        attr(files, "sources") <- lapply(seq_len(nrow(fi)), function(i)
            .mzs_external_source(run[i], location = fi$filepth[i]))
    }
    files
}

#' A deterministic run id from a file's identity and ordinal, never its name.
#'
#' @noRd
.mzs_mint_run_id <- function(identity, ordinal) {
    paste0("run-", substr(.mzs_sha256_text(paste0(identity, "\n", ordinal)),
                          1L, 12L))
}

#' MS2 scan references: source, run, study spectrum_id_ and mzML native id.
#'
#' Study scans are matched by native id, then scan number. Without a study,
#' spectrum_id_ is NA.
#'
#' @param scans data.frame with fileid, seqNum, acquisitionNum.
#'
#' @noRd
.mzs_scan_refs <- function(files, scans) {
    k <- match(scans$fileid, files$fileid)
    out <- data.frame(source = files$source[k], run_id = files$run_id[k],
                      spectrum_id_ = NA_real_, native_id = NA_character_,
                      stringsAsFactors = FALSE)
    for (j in which(files$readable)) {
        i <- which(scans$fileid == files$fileid[j])
        if (!length(i))
            next
        h <- tryCatch(.msp_header(.msp_read(files$filepth[j])),
                      error = function(e) NULL)
        if (is.null(h) || is.null(h$spectrumId))
            next
        out$native_id[i] <- as.character(
            h$spectrumId[match(scans$seqNum[i], h$seqNum)])
    }
    study <- attr(files, "study")
    if (!is.null(study)) {
        sm <- attr(files, "study_manifest")
        cols <- c("id", "scan_number", "acquisition_num_")
        meta <- .mzs_read_spectra_meta(study, sm, cols,
                                       unique(out$run_id))
        sn <- meta$scan_number
        sn[is.na(sn)] <- meta$acquisition_num_[is.na(sn)]
        by_id <- match(paste(out$run_id, out$native_id),
                       paste(meta$run_id, meta$id))
        by_id[is.na(out$native_id)] <- NA_integer_
        by_scan <- match(paste(out$run_id, scans$acquisitionNum),
                         paste(meta$run_id, sn))
        hit <- ifelse(is.na(by_id), by_scan, by_id)
        miss <- !is.na(out$run_id) & is.na(hit)
        if (any(miss))
            .mzs_abort("semantic", sum(miss), " MS2 scan(s) are not in the ",
                       "study run they were recorded against, e.g. run ",
                       out$run_id[miss][1L], " scan ",
                       scans$acquisitionNum[miss][1L], ".")
        out$spectrum_id_ <- as.numeric(meta$spectrum_id_[hit])
        nat <- is.na(out$native_id)
        out$native_id[nat] <- as.character(meta$id[hit[nat]])
    }
    out
}

# ---------------------------------------------------------------------------
# Derived spectra.
# ---------------------------------------------------------------------------

#' Scan polarity as +1, -1 or NA.
#'
#' @noRd
.mzs_scan_polarity <- function(x) {
    x <- tolower(trimws(as.character(x)))
    out <- rep(NA_integer_, length(x))
    out[x %in% c("positive", "pos", "+", "1")] <- 1L
    out[x %in% c("negative", "neg", "-", "-1", "0")] <- -1L
    out
}

#' Averaged spectra in write order, with their run and spectrum id.
#'
#' @noRd
.mzs_derived_index <- function(src, files, uid) {
    m <- src$s_peak_meta
    type <- m$spectrum_type
    av <- m[!is.na(type) & type %in% c("intra", "inter", "all"), ,
            drop = FALSE]
    if (!nrow(av))
        return(av)
    ## Without any "intra" rows, "all" rows with a file id are intra-file.
    av$method <- ifelse(av$spectrum_type == "all" & !is.na(av$fileid) &
                        !any(m$spectrum_type %in% "intra"), "intra",
                        av$spectrum_type)
    k <- match(av$fileid, files$fileid)
    intra <- av$method == "intra"
    if (any(intra & is.na(k)))
        .mzs_abort("semantic", "Averaged spectra name file(s) absent from ",
                   "fileinfo.")
    av$run_id <- ifelse(intra, paste0("av_intra_", files$run_id[k]),
                        paste0("av_", av$method))
    av$data_origin <- ifelse(intra, files$data_origin[k],
                             paste0("mspurity:", uid, ":av_", av$method))
    av$source_run <- ifelse(intra, files$run_id[k], NA_character_)
    av$source_key <- ifelse(intra, files$source[k], NA_character_)
    ord <- order(match(av$method, c("intra", "inter", "all")),
                 match(av$source_run, files$run_id), av$grpid, av$pid)
    av <- av[ord, , drop = FALSE]
    av$spectrum_id_ <- seq_len(nrow(av))
    av$native_id <- paste0("msPurity:av_", av$method, ":", av$grpid)
    rownames(av) <- NULL
    av
}

#' The native runs holding the averaged spectra, for `.mzs_commit_new()`.
#'
#' @noRd
.mzs_derived_runs <- function(src, av) {
    if (!nrow(av))
        return(list(runs = list(), results_runs = list(),
                    peak_columns = list()))
    pk <- src$s_peaks[src$s_peaks$pid %in% av$pid, , drop = FALSE]
    pk <- pk[order(pk$pid, pk$mz, pk$cl), , drop = FALSE]
    unknown <- setdiff(names(pk), c("sid", "fileid", "mz", "i", "type",
                                    "scan", "grpid", "pid", "frac",
                                    "pass_flag", "purity_pass_flag",
                                    "intensity_pass_flag",
                                    names(.MZS_AV_PEAKS)))
    if (length(unknown))
        .mzs_abort("unsupported", "Averaged spectra carry per-peak ",
                   "column(s) this converter does not map: ",
                   paste(unknown, collapse = ", "), " (R-081).")
    present <- names(.MZS_AV_PEAKS)[names(.MZS_AV_PEAKS) %in% names(pk) &
        vapply(names(.MZS_AV_PEAKS), function(c)
            c %in% names(pk) && !all(is.na(pk[[c]])), logical(1))]
    runs <- list()
    results_runs <- list()
    peak_columns <- list()
    for (r in unique(av$run_id)) {
        a <- av[av$run_id == r, , drop = FALSE]
        rows <- split(seq_len(nrow(pk)), factor(pk$pid, levels = a$pid))
        col <- function(c, f) lapply(rows, function(i) {
            v <- f(pk[[c]][i])
            if (all(is.na(v))) NULL else v
        })
        peaks <- list(
            mz = .mzs_peak_list(unname(lapply(rows, function(i)
                as.numeric(pk$mz[i])))),
            intensity = .mzs_peak_list(unname(lapply(rows, function(i)
                as.numeric(pk$i[i])))))
        for (c in present) {
            spec <- .MZS_AV_PEAKS[[c]]
            f <- switch(spec[[2L]], i32 = as.integer, bool = as.logical,
                        as.numeric)
            peaks[[spec[[1L]]]] <- .mzs_peak_list(unname(col(c, f)),
                                                  spec[[2L]])
        }
        meta <- data.frame(
            ms_level = rep(2L, nrow(a)),
            time = a$retention_time / 60,
            scan_polarity = .mzs_scan_polarity(a$polarity),
            spectrum_representation = rep("MS:1000127", nrow(a)),
            id = a$native_id,
            data_origin = a$data_origin,
            selected_ion_mz = a$precursor_mz,
            x_mspurity_precursor_purity = a$inPurity,
            stringsAsFactors = FALSE)
        ## createDatabase()'s `metadata` fills the instrument fields.
        for (c in c("instrument", "instrument_type"))
            if (c %in% names(av) && any(!is.na(av[[c]])))
                meta[[paste0("x_mspurity_", c)]] <- a[[c]]
        runs[[length(runs) + 1L]] <- list(run_id = r, meta = meta,
                                          peaks = peaks)
        m1 <- a$method[1L]
        results_runs[[r]] <- if (m1 == "intra")
            list(method = "intra", scope = "run", source = a$source_key[1L],
                 source_run_id = a$source_run[1L])
        else
            list(method = m1, scope = "dataset", source = NULL,
                 source_run_id = NULL)
        peak_columns[[r]] <- setdiff(names(peaks), c("mz", "intensity"))
    }
    list(runs = runs, results_runs = results_runs,
         peak_columns = peak_columns)
}

#' The merge_member table: the scans each averaged spectrum drew on.
#'
#' All-file averages draw on the feature's passing scans, intra-file averages
#' on those of their file, inter-file averages on the feature's intra-file
#' averages. A scan linked twice is listed twice. Membership is complete when
#' it matches the recorded contributor count.
#'
#' @param members data.frame of contributions: grpid, pid, fileid, one row
#'     per scan-to-feature link.
#'
#' @noRd
.mzs_merge_member <- function(src, av, members, scans, scan_refs) {
    if (!nrow(av))
        return(list(table = NULL, incomplete = 0L))
    avp <- src$s_peaks[src$s_peaks$pid %in% av$pid & !is.na(src$s_peaks$total),
                       c("pid", "total")]
    total <- tapply(avp$total, avp$pid, max)
    out <- lapply(seq_len(nrow(av)), function(i) {
        a <- av[i, ]
        if (a$method == "inter") {
            mem <- av[av$method == "intra" & av$grpid == a$grpid, ,
                      drop = FALSE]
            if (!nrow(mem))
                return(NULL)
            return(data.frame(
                merged_spectrum_id_ = a$spectrum_id_,
                merged_run_id = a$run_id, merged_native_id = a$native_id,
                member_source = "self", member_run_id = mem$run_id,
                member_spectrum_id_ = mem$spectrum_id_,
                member_native_id = mem$native_id,
                stringsAsFactors = FALSE))
        }
        sel <- which(members$grpid == a$grpid &
                     (a$method == "all" | members$fileid %in% a$fileid))
        ## Passing scans if they match the contributor count, else all.
        want <- total[as.character(a$pid)]
        pass <- sel[members$passes[sel] %in% TRUE]
        if (length(pass) && (is.na(want) || length(pass) == want ||
                             length(sel) != want))
            sel <- pass
        if (!length(sel))
            return(NULL)
        k <- match(members$pid[sel], scans$pid)
        o <- order(match(scans$fileid[k], sort(unique(scans$fileid))),
                   scans$acquisitionNum[k], members$pid[sel])
        k <- k[o]
        data.frame(
            merged_spectrum_id_ = a$spectrum_id_,
            merged_run_id = a$run_id, merged_native_id = a$native_id,
            member_source = scan_refs$source[k],
            member_run_id = scan_refs$run_id[k],
            member_spectrum_id_ = scan_refs$spectrum_id_[k],
            member_native_id = scan_refs$native_id[k],
            x_mspurity_member_scan_number = scans$acquisitionNum[k],
            stringsAsFactors = FALSE)
    })
    mm <- do.call(.mzs_rbind_fill, out)
    if (is.null(mm) || !nrow(mm))
        return(list(table = NULL, incomplete = 0L))
    mm$merged_source <- "self"
    mm$member_rank <- as.integer(stats::ave(mm$merged_spectrum_id_,
                                            mm$merged_spectrum_id_,
                                            FUN = seq_along))
    want <- suppressWarnings(as.integer(total[as.character(
        av$pid[match(mm$merged_spectrum_id_, av$spectrum_id_)])]))
    got <- as.integer(stats::ave(mm$member_rank, mm$merged_spectrum_id_,
                                 FUN = length))
    mm$members_complete <- !is.na(want) & got == want
    list(table = mm,
         incomplete = sum(!tapply(mm$members_complete,
                                  mm$merged_spectrum_id_, all)))
}

#' Row-bind data.frames whose columns differ, filling with NA.
#'
#' @noRd
.mzs_rbind_fill <- function(...) {
    dfs <- Filter(function(d) !is.null(d) && nrow(d), list(...))
    if (!length(dfs))
        return(NULL)
    cols <- unique(unlist(lapply(dfs, names)))
    do.call(rbind, lapply(dfs, function(d) {
        for (c in setdiff(cols, names(d)))
            d[[c]] <- rep(NA, nrow(d))
        d[, cols, drop = FALSE]
    }))
}

#' MS2 scans as recorded: s_peak_meta rows that are not averages.
#'
#' @noRd
.mzs_scans <- function(src) {
    m <- src$s_peak_meta
    m[is.na(m$spectrum_type) | m$spectrum_type == "scan", , drop = FALSE]
}

#' MS2 scan-to-feature links: grpid, pid and cid (NA for feature-width links).
#'
#' @noRd
.mzs_scan_links <- function(src) {
    if (!is.null(src$c_peak_X_s_peak_meta) &&
        nrow(src$c_peak_X_s_peak_meta)) {
        link <- merge(src$c_peak_X_s_peak_meta[, c("pid", "cid")],
                      src$c_peak_X_c_peak_group[, c("cid", "grpid")],
                      by = "cid")
        return(link[order(link$pid, link$cid, link$grpid), , drop = FALSE])
    }
    g <- src$c_peak_group_X_s_peak_meta
    if (is.null(g))
        return(data.frame(pid = integer(), cid = integer(),
                          grpid = integer()))
    data.frame(pid = g$pid, cid = NA_integer_, grpid = g$grpid)
}

#' The s_peak_meta pid of each scan peak in s_peaks.
#'
#' s_peaks `pid` does not identify the scan, so match on file id and `scan`
#' (the acquisition number).
#'
#' @noRd
.mzs_scan_peak_pid <- function(src, sp) {
    m <- .mzs_scans(src)
    k <- match(paste(sp$fileid, sp$scan), paste(m$fileid, m$acquisitionNum))
    m$pid[k]
}

#' Candidate contributions per feature from the database: every
#' scan-to-feature link, with `passes` TRUE if the scan has a passing peak
#' (NA when scan peaks carry no flags).
#'
#' @noRd
.mzs_members_from_db <- function(src, scans) {
    link <- .mzs_scan_links(src)
    sp <- src$s_peaks[src$s_peaks$type %in% "scan", , drop = FALSE]
    passes <- rep(NA, nrow(link))
    if ("pass_flag" %in% names(sp) && any(!is.na(sp$pass_flag))) {
        pid <- .mzs_scan_peak_pid(src, sp)
        passes <- link$pid %in% unique(pid[sp$pass_flag %in% TRUE])
    }
    data.frame(grpid = link$grpid, pid = link$pid,
               fileid = scans$fileid[match(link$pid, scans$pid)],
               passes = passes)
}

#' Contributing scans per feature from the live purityA object: grped_df rows
#' with a peak passing filterFragSpectra() (any peak if unfiltered).
#'
#' @noRd
.mzs_members_from_pa <- function(pa) {
    g <- pa@grped_df
    if (!nrow(g))
        return(data.frame(grpid = integer(), pid = integer(),
                          fileid = integer()))
    filtered <- length(pa@filter_frag_params) > 0L
    keep <- logical(nrow(g))
    ms2 <- groupedSpectra(pa)
    for (grp in unique(as.character(g$grpid))) {
        rows <- which(as.character(g$grpid) == grp)
        spectra <- ms2[[grp]]
        for (j in seq_along(rows)) {
            s <- spectra[[j]]
            keep[rows[j]] <- !is.null(s) && NROW(s) > 0L &&
                (!filtered || any(s[, ncol(s)] == 1))
        }
    }
    g <- g[keep, , drop = FALSE]
    fileid <- if ("sample" %in% names(g)) g$sample else g$fileid
    data.frame(grpid = as.integer(as.character(g$grpid)),
               pid = as.integer(g$pid),
               fileid = as.integer(as.character(fileid)),
               passes = TRUE)
}

# ---------------------------------------------------------------------------
# Raw MS2 scans: msPurity's per-scan and per-peak annotation.
# ---------------------------------------------------------------------------

.MZS_SCAN_COLUMNS <- c(
    precursorMZ = "precursor_mz", precursorRT = "precursor_retention_time",
    precursorIntensity = "precursor_intensity",
    retentionTime = "retention_time", acquisitionNum = "scan_number",
    precursorScanNum = "precursor_scan_number",
    precursorNearest = "precursor_nearest", aMz = "a_mz",
    aPurity = "a_purity", apkNm = "a_peak_count", iMz = "i_mz",
    iPurity = "i_purity", ipkNm = "i_peak_count", inPurity = "in_purity",
    inPkNm = "in_peak_count", purity_pass_flag = "purity_pass_flag")

# MassBank-style s_peak_meta columns. Only polarity and instrument fields are
# ever filled.
.MZS_LIBRARY_META <- c("name", "collision_energy", "ms_level", "accession",
                       "resolution", "polarity", "fragmentation_type",
                       "precursor_type", "instrument_type", "instrument",
                       "copyright", "column", "mass_accuracy", "mass_error",
                       "origin", "splash", "retention_index", "inchikey_id",
                       "sourceid")

#' x_mspurity_scan (assessed MS2 scans with precursor purity) and, if
#' filterFragSpectra() flagged peaks, x_mspurity_scan_peak.
#'
#' @noRd
.mzs_scan_tables <- function(src, scans, scan_refs) {
    n <- nrow(scans)
    scan <- data.frame(
        scan_annotation_id_ = seq_len(n),
        ms2_source = scan_refs$source, ms2_run_id = scan_refs$run_id,
        ms2_spectrum_id_ = scan_refs$spectrum_id_,
        ms2_native_id = scan_refs$native_id,
        x_mspurity_seq_num = scans$seqNum,
        stringsAsFactors = FALSE)
    for (c in intersect(names(.MZS_SCAN_COLUMNS), names(scans)))
        scan[[.MZS_SCAN_COLUMNS[[c]]]] <- scans[[c]]
    sp <- src$s_peaks[src$s_peaks$type %in% "scan", , drop = FALSE]
    ann <- intersect(c("snr", "ra", "purity_pass_flag",
                       "intensity_pass_flag", "ra_pass_flag",
                       "snr_pass_flag"), names(sp))
    ann <- ann[vapply(ann, function(c) any(!is.na(sp[[c]])), logical(1))]
    peak <- NULL
    if (length(ann) && nrow(sp)) {
        pid <- .mzs_scan_peak_pid(src, sp)
        if (anyNA(pid))
            .mzs_abort("semantic", sum(is.na(pid)), " scan peak(s) in ",
                       "s_peaks name no scan of s_peak_meta (by file and ",
                       "scan number).")
        peak <- data.frame(
            scan_annotation_id_ = match(pid, scans$pid),
            peak_rank = as.integer(stats::ave(seq_len(nrow(sp)), pid,
                                              FUN = seq_along)),
            mz = sp$mz, intensity = sp$i,
            snr = if ("snr" %in% ann) sp$snr else NA_real_,
            relative_intensity = if ("ra" %in% ann) sp$ra else NA_real_,
            stringsAsFactors = FALSE)
        for (f in c("purity_pass_flag", "intensity_pass_flag",
                    "ra_pass_flag", "snr_pass_flag"))
            peak[[f]] <- if (f %in% ann) sp[[f]] else NA
    }
    list(scan = scan, peak = peak, raw_peaks = nrow(sp),
         annotated = !is.null(peak))
}

# ---------------------------------------------------------------------------
# The conversion.
# ---------------------------------------------------------------------------

#' tool_provenance steps. Parameters are known only on the live route.
#'
#' @noRd
.mzs_mspurity_steps <- function(src, live, version) {
    has_av <- any(src$s_peak_meta$spectrum_type %in% c("intra", "inter",
                                                        "all"))
    has_filter <- any(c("pass_flag", "snr_pass_flag") %in%
                      names(src$s_peaks))
    step <- function(name, params, ok = TRUE)
        list(tool_name = "msPurity", tool_version = version,
             step_name = name, parameters = params,
             completeness = if (is.null(params)) "absent"
                            else if (ok) "complete" else "partial")
    if (is.null(live)) {
        steps <- list(step("purityA", NULL), step("frag4feature", NULL))
        if (has_filter)
            steps <- c(steps, list(step("filterFragSpectra", NULL)))
        if (has_av)
            steps <- c(steps, list(step("averageFragSpectra", NULL)))
        return(c(steps, list(step("createDatabase", NULL))))
    }
    pa <- live$pa
    prm <- if (methods::.hasSlot(pa, "params")) pa@params else list()
    list(
        step("purityA", prm$purityA %||%
                 list(mzRback = pa@mzRback, cores = pa@cores),
             ok = !is.null(prm$purityA)),
        step("frag4feature", prm$frag4feature %||%
                 list(f4f_link_type = pa@f4f_link_type),
             ok = !is.null(prm$frag4feature)),
        if (length(pa@filter_frag_params))
            step("filterFragSpectra", pa@filter_frag_params),
        if (length(pa@av_intra_params))
            step("averageIntraFragSpectra", pa@av_intra_params),
        if (length(pa@av_inter_params))
            step("averageInterFragSpectra", pa@av_inter_params),
        if (length(pa@av_all_params))
            step("averageAllFragSpectra", pa@av_all_params),
        step("createDatabase", list(
            metadata = if (is.list(live$metadata)) live$metadata
                       else .mzs_object())))
}

#' Coverage manifest: the disposition of each msPurity construct.
#'
#' @noRd
.mzs_mspurity_coverage <- function(src, live) {
    has <- function(t) !is.null(src[[t]]) && nrow(src[[t]]) > 0L
    link <- if (has("c_peak_X_s_peak_meta")) "c_peak_X_s_peak_meta"
            else "c_peak_group_X_s_peak_meta"
    cov <- c(fileinfo = "mapped",
             c_peaks = "mapped",
             c_peak_groups = "mapped_partial",
             c_peak_X_c_peak_group = "mapped",
             s_peak_meta = "mapped_partial",
             s_peaks = "mapped_partial",
             source = "mapped",
             metab_compound = if (has("metab_compound")) "mapped"
                              else "not_present_in_source",
             sm_matches = if (has("sm_matches")) "mapped"
                          else "not_present_in_source",
             l_s_peak_meta = if (has("l_s_peak_meta")) "mapped_partial"
                             else "not_present_in_source",
             xcms_match = if (has("xcms_match")) "dropped"
                          else "not_present_in_source",
             parameters = if (is.null(live)) "not_present_in_source"
                          else "mapped")
    cov[[link]] <- "mapped"
    cov
}

#' Loss ledger rows for the msPurity tables.
#'
#' @noRd
.mzs_mspurity_losses <- function(src, live, scan_tabs, incomplete) {
    pk <- src$s_peaks
    av <- pk[!pk$type %in% "scan", , drop = FALSE]
    out <- list(
        if ("frac" %in% names(av) && nrow(av))
            list(construct = "s_peaks.frac", disposition =
                     "dropped_by_configuration", affected_rows = nrow(av),
                 reason = paste("contributor_count / contributor_total,",
                                "recomputable exactly (R-060).")),
        if ("pass_flag" %in% names(pk))
            list(construct = "s_peaks.pass_flag",
                 disposition = "dropped_by_configuration",
                 affected_rows = sum(!is.na(pk$pass_flag)),
                 reason = paste("A composite whose rule differs between",
                                "scans and averaged spectra; its component",
                                "flags are stored and both rules recorded",
                                "in provenance (pass_criteria).")),
        if (!scan_tabs$annotated && scan_tabs$raw_peaks)
            list(construct = "s_peaks (type = scan): mz, i",
                 disposition = "no_such_concept_in_target",
                 affected_rows = scan_tabs$raw_peaks,
                 reason = paste("Raw scan peaks re-read from the mzML files;",
                                "raw scans are referenced, not copied",
                                "(R-009).")),
        if ("sid" %in% names(pk))
            list(construct = "s_peaks.sid",
                 disposition = "dropped_by_configuration",
                 affected_rows = nrow(pk),
                 reason = "A positional row number with no identity."),
        if (is.null(live))
            list(construct = "processing parameters",
                 disposition = "not_present_in_source",
                 target_ref = "tool_provenance",
                 reason = paste("createDatabase() does not write msPurity's",
                                "purity, linking, filtering or averaging",
                                "parameters (R-090).")),
        if (incomplete > 0L)
            list(construct = "averaged spectrum membership",
                 disposition = "not_present_in_source",
                 affected_rows = incomplete, target_ref = "merge_member",
                 reason = paste("msPurity records how many scans each",
                                "average drew on, not which; membership is",
                                "reconstructed from the scan-to-feature",
                                "links and marked incomplete where the",
                                "two disagree.")),
        list(construct = "source.name", disposition =
                 "no_such_concept_in_target", affected_rows = 1,
             reason = "The database's own name; the conversion row records its location."))
    m <- src$s_peak_meta
    for (c in intersect(.MZS_LIBRARY_META, names(m))) {
        n <- sum(!is.na(m[[c]]))
        mapped <- c %in% c("polarity", "instrument", "instrument_type")
        if (mapped && n)
            next
        out[[length(out) + 1L]] <- list(
            construct = paste0("s_peak_meta.", c),
            disposition = if (n) "no_such_concept_in_target"
                          else "not_present_in_source",
            affected_rows = n,
            reason = if (n) "Library-record metadata with no place in mzStack-4."
                     else "Always empty in msPurity's output.")
    }
    out
}

#' Loss ledger rows for the feature tables.
#'
#' @noRd
.mzs_feature_losses <- function(src, ftabs, live) {
    grp <- src$c_peak_groups
    out <- list(
        if ("peakidx" %in% names(grp))
            list(construct = "c_peak_groups.peakidx",
                 disposition = "dropped_by_configuration",
                 affected_rows = nrow(grp),
                 target_ref = "chromatographic_peak_feature",
                 reason = paste("The feature's peaks, recorded in full in",
                                "chromatographic_peak_feature.")),
        if (length(ftabs$sample_columns$classes))
            list(construct = paste0("c_peak_groups.",
                                    ftabs$sample_columns$classes,
                                    collapse = ", "),
                 disposition = "dropped_by_configuration",
                 affected_rows = nrow(grp),
                 reason = paste("Peak counts per sample class, recomputable",
                                "from chromatographic_peak_feature and",
                                "assay.")),
        list(construct = "c_peak_X_c_peak_group.cXg_id",
             disposition = "dropped_by_configuration",
             affected_rows = nrow(src$c_peak_X_c_peak_group),
             reason = "A positional row number of the link table."),
        if ("sn" %in% names(src$c_peaks))
            list(construct = "c_peaks.sn", disposition = "preserved_opaque",
                 affected_rows = nrow(src$c_peaks),
                 target_ref = "chromatographic_peak.x_mspurity_sn",
                 reason = paste("signal_to_noise is float32 in mzStack-4;",
                                "the exact value is carried in",
                                "x_mspurity_sn.")),
        list(construct = "s_peak_meta.inPurity (linked scans)",
             disposition = "preserved_opaque",
             affected_rows = nrow(.mzs_scan_links(src)),
             target_ref = "x_mspurity_scan.in_purity",
             reason = paste("spectrum_feature.precursor_purity is float32 in",
                            "mzStack-4; the exact value is in",
                            "x_mspurity_scan.in_purity.")),
        if (!is.null(live))
            list(construct = "grped_df.precurMtchPPM",
                 disposition = "preserved_opaque",
                 target_ref = paste0("spectrum_feature.",
                                     "x_mspurity_precursor_mz_error_exact"),
                 reason = paste("precursor_mz_error is float32 in mzStack-4;",
                                "the exact value is carried beside it.")))
    out
}

#' Reasons for wholly-null columns.
#'
#' @noRd
.MZS_NULL_REASONS <- list(
    "*.ms2_usi" = list(disposition = "not_present_in_source",
                       reason = "No source declares a USI collection."),
    "*.ms2_native_id" = list(
        disposition = "not_present_in_source",
        reason = paste("The source files could not be read at conversion,",
                       "and no study dataset records the native ids.")),
    "*.member_native_id" = list(
        disposition = "not_present_in_source",
        reason = paste("The source files could not be read at conversion,",
                       "and no study dataset records the native ids.")),
    "tool_provenance.parameters" = list(
        disposition = "not_present_in_source",
        reason = "createDatabase() does not write processing parameters."),
    "conversion.source_uri" = list(
        disposition = "no_such_concept_in_source",
        reason = "The live-object route converts no file."),
    "conversion.source_checksum" = list(
        disposition = "no_such_concept_in_source",
        reason = "The live-object route converts no file."),
    "conversion.source_checksum_algorithm" = list(
        disposition = "no_such_concept_in_source",
        reason = "The live-object route converts no file."),
    "abundance.chromatographic_peak_id_" = list(
        disposition = "dropped_by_configuration",
        reason = paste("Omitted by default (mzStack-4 section 5.4): recoverable",
                       "through chromatographic_peak_feature.")),
    "sample.sample_class" = list(
        disposition = "not_present_in_source",
        reason = "No sample class was set in xcms."),
    "x_mspurity_scan.purity_pass_flag" = list(
        disposition = "not_present_in_source",
        reason = paste("filterFragSpectra() did not assess every scan",
                       "(allfrag = FALSE).")),
    "spectrum_feature.precursor_mz_error" = list(
        disposition = "not_present_in_source",
        reason = "createDatabase() does not record the precursor match error."),
    "*.member_spectrum_id_" = list(
        disposition = "not_present_in_source",
        reason = paste("The scans are in external sources (no study",
                       "dataset), so no spectrum id is assigned (R-026).")),
    "*.ms2_spectrum_id_" = list(
        disposition = "not_present_in_source",
        reason = paste("The scans are in external sources (no study",
                       "dataset), so no spectrum id is assigned (R-026).")))

#' Convert a typed msPurity table set into a Parquet results dataset.
#'
#' @param live `NULL` for the SQLite route, or `list(pa, xcmsObj, metadata)`.
#'
#' @noRd
.mzs_mspurity_convert <- function(src, path, study, studyKey, fileMap,
                                  overwrite, route, live, fn, started,
                                  inputs, parameters, source_uri,
                                  checksum = NA_character_) {
    uid <- .mzs_uid()
    version <- if (is.null(live))
        sub("^.* ", "", src$source$parsing_software[1L] %||% NA_character_)
    else .mzs_pkg_version("msPurity")
    files <- .mzs_files(src$fileinfo, study, studyKey, fileMap, path)
    scans <- .mzs_scans(src)
    scan_refs <- .mzs_scan_refs(files, scans)
    av <- .mzs_derived_index(src, files, uid)
    derived <- .mzs_derived_runs(src, av)
    members <- if (!is.null(live)) .mzs_members_from_pa(live$pa)
               else .mzs_members_from_db(src, scans)
    mm <- .mzs_merge_member(src, av, members, scans, scan_refs)
    st <- .mzs_scan_tables(src, scans, scan_refs)
    ftabs <- .mzs_feature_tables(src, files)
    sf <- .mzs_spectrum_feature(
        src, scans, scan_refs, av, ftabs,
        match_details = if (!is.null(live)) .mzs_match_details(live$pa))

    tables <- c(ftabs$tables, list(spectrum_feature = sf,
                                   x_mspurity_scan = st$scan))
    evd <- if (is.null(live)) .mzs_sqlite_evidence(src, files, av, scans,
                                                     scan_refs, ftabs)
    tables <- c(tables, evd$tables)
    if (!is.null(st$peak))
        tables$x_mspurity_scan_peak <- st$peak
    if (!is.null(mm$table))
        tables$merge_member <- mm$table

    tables$conversion <- .mzs_conversion_row(
        route = paste0("mspurity:", route),
        source_format = if (is.null(live)) "msPurity createDatabase SQLite"
                        else "msPurity purityA",
        source_format_version = version,
        coverage = .mzs_mspurity_coverage(src, live),
        source_uri = source_uri, checksum = checksum)
    tables$tool_provenance <- .mzs_tool_provenance(
        .mzs_mspurity_steps(src, live, version))
    tables$loss_ledger <- .mzs_loss_ledger(c(
        .mzs_mspurity_losses(src, live, st, mm$incomplete),
        .mzs_feature_losses(src, ftabs, live),
        evd$losses,
        .mzs_null_column_losses(tables, .MZS_NULL_REASONS)))
    tables$source_identifier <- .mzs_bind_source_ids(
        .mzs_source_ids("spectra", av$spectrum_id_, "s_peak_meta", "pid",
                        src$lexical$pid[match(av$pid, src$s_peak_meta$pid)],
                        "positional"),
        .mzs_source_ids("x_mspurity_scan", st$scan$scan_annotation_id_,
                        "s_peak_meta", "pid",
                        src$lexical$pid[match(scans$pid,
                                              src$s_peak_meta$pid)],
                        "positional"),
        .mzs_source_ids("feature", ftabs$tables$feature$feature_id_,
                        "c_peak_groups", "grpid",
                        src$lexical$grpid[match(ftabs$groups$grpid,
                                                src$c_peak_groups$grpid)],
                        "positional"),
        .mzs_source_ids("feature", ftabs$tables$feature$feature_id_,
                        "c_peak_groups", "grp_name", ftabs$groups$grp_name,
                        "name"),
        .mzs_source_ids("chromatographic_peak",
                        ftabs$tables$chromatographic_peak$chromatographic_peak_id_,
                        "c_peaks", "cid",
                        src$lexical$cid[match(ftabs$peaks$cid,
                                              src$c_peaks$cid)],
                        "positional"),
        .mzs_source_ids("assay", seq_len(nrow(files)), "fileinfo", "fileid",
                        src$lexical$fileid[match(files$fileid,
                                                 src$fileinfo$fileid)],
                        "positional"),
        .mzs_source_ids("assay", seq_len(nrow(files)), "fileinfo", "filepth",
                        files$filepth, "path"),
        .mzs_source_ids("assay", seq_len(nrow(files)), "fileinfo", "nm_save",
                        src$fileinfo$nm_save[match(files$fileid,
                                                   src$fileinfo$fileid)],
                        "name"),
        .mzs_source_ids("sample", seq_len(nrow(files)), "fileinfo",
                        "filename", files$filename, "name"),
        .mzs_source_ids("spectrum_feature",
                        sf$association_id_[seq_along(attr(sf, "link_id"))],
                        attr(sf, "link_table"), attr(sf, "link_column"),
                        attr(sf, "link_id"), "positional"),
        evd$source_ids)

    plan <- list(
        uid = uid, role = "results",
        runs = derived$runs,
        run_types = c(.MZS_AV_META_TYPES, x_mspurity_instrument = "str",
                      x_mspurity_instrument_type = "str"),
        results_runs = derived$results_runs,
        peak_columns = derived$peak_columns,
        tables = tables,
        fill_policy = "dense",
        sources = c(attr(files, "sources"), evd$sources),
        selection = if (!is.null(evd)) list(evidence = .mzs_selection(NA)),
        activity = list(
            action = "convert", started = started, fn = fn,
            inputs = c(lapply(attr(files, "sources"), function(x)
                list(source = x$key)), inputs),
            parameters = c(parameters, list(
                route = paste0("mspurity:", route),
                fileMap = if (length(fileMap)) as.list(fileMap),
                pass_criteria = .MZS_PASS_CRITERIA,
                precursor_definition = .MZS_PRECURSOR_DEFINITION))))
    .mzs_commit_new(path, plan, overwrite)
    invisible(path)
}

#' createDatabase(format = "parquet").
#'
#' @noRd
.createDatabase_parquet <- function(pa, xcmsObj, xsa, outDir, grpPeaklist,
                                    dbName, metadata, xset, study, studyKey,
                                    fileMap, overwrite) {
    .mzs_require()
    started <- .mzs_now()
    if (!is.null(xsa))
        .mzs_abort("unsupported", "createDatabase(format = \"parquet\") does ",
                   "not map CAMERA annotation ('xsa').",
                   data = list(constructs = "xsa"))
    if (is.data.frame(grpPeaklist))
        .mzs_abort("unsupported", "createDatabase(format = \"parquet\") ",
                   "takes features from 'xcmsObj'; a custom 'grpPeaklist' ",
                   "is not mapped.", data = list(constructs = "grpPeaklist"))
    if (is.na(dbName))
        dbName <- paste0("lcmsms_data-", format(Sys.time(),
                                                "%Y-%m-%d-%H%M%S"),
                         ".parquet")
    path <- file.path(outDir, dbName)
    if (!is.null(study))
        study <- normalizePath(study, mustWork = TRUE)
    .mzs_check_destination(path, overwrite, study)

    tmp <- tempfile("mzs-createDatabase-")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    db <- .createDatabase_sqlite(pa = pa, xcmsObj = xcmsObj, xsa = NULL,
                                 outDir = tmp, grpPeaklist = NA,
                                 dbName = "staging.sqlite",
                                 metadata = metadata, xset = xset)
    if (is.null(db))
        .mzs_abort("semantic", "createDatabase() could not match the files ",
                   "of 'pa' and 'xcmsObj'.")
    src <- .mzs_read_sqlite(db)
    message("Writing the Parquet dataset '", path, "'")
    .mzs_mspurity_convert(
        src, path, study, studyKey, fileMap, overwrite, route = "live-object",
        live = list(pa = pa, xcmsObj = xcmsObj, metadata = metadata),
        fn = "createDatabase", started = started, inputs = list(),
        parameters = list(), source_uri = NA_character_)
    invisible(normalizePath(path))
}

#' Precursor match details per MS2 scan from frag4feature().
#'
#' @noRd
.mzs_match_details <- function(pa) {
    g <- pa@grped_df
    need <- c("pid", "grpid", "precurMtchMZ", "precurMtchPPM", "precurMtchRT")
    if (!nrow(g) || !all(need %in% names(g)))
        return(NULL)
    data.frame(pid = as.integer(g$pid),
               cid = if ("cid" %in% names(g)) as.integer(g$cid)
                     else NA_integer_,
               grpid = as.integer(as.character(g$grpid)),
               precurMtchMZ = as.numeric(g$precurMtchMZ),
               precurMtchPPM = as.numeric(g$precurMtchPPM),
               precurMtchRT = as.numeric(g$precurMtchRT))
}

#' Spectral-matching results in a database: evidence, scores, compounds and
#' coverage.
#'
#' Coverage lists matched queries only, since attempted queries are not
#' recorded. Library spectra are referenced by accession, one external source
#' per library.
#'
#' @noRd
.mzs_sqlite_evidence <- function(src, files, av, scans, scan_refs, ftabs) {
    sm <- src$sm_matches
    if (is.null(sm) || !nrow(sm))
        return(NULL)
    grp_ids <- ftabs$groups$grpid
    ka <- match(sm$qpid, av$pid)
    ks <- match(sm$qpid, scans$pid)
    if (any(is.na(ka) & is.na(ks)))
        .mzs_abort("semantic", "sm_matches names query spectra absent from ",
                   "s_peak_meta.")
    self <- !is.na(ka)
    link <- .mzs_scan_links(src)
    scan_feature <- tapply(link$grpid, link$pid, function(g)
        if (length(unique(g)) == 1L) g[1L] else NA)
    lib <- .mzs_na_empty(sm$library_source_name)
    lib[is.na(lib)] <- "unknown"
    key <- paste0("library_", gsub("[^A-Za-z0-9._-]+", "_", lib))
    matches <- data.frame(
        qkey = paste0("q", sm$qpid), order = seq_len(nrow(sm)),
        query_source = ifelse(self, "self", scan_refs$source[ks]),
        query_run_id = ifelse(self, av$run_id[ka], scan_refs$run_id[ks]),
        query_spectrum_id_ = ifelse(self, av$spectrum_id_[ka],
                                    scan_refs$spectrum_id_[ks]),
        query_native_id = ifelse(self, av$native_id[ka],
                                 scan_refs$native_id[ks]),
        feature_id_ = match(ifelse(self, av$grpid[ka],
                                   scan_feature[as.character(sm$qpid)]),
                            grp_ids),
        reference_source = key, reference_run_id = NA_character_,
        reference_spectrum_id_ = NA_real_,
        reference_native_id = .mzs_na_empty(sm$library_accession),
        dpc = sm$dpc, rdpc = sm$rdpc, cdpc = sm$cdpc, mcount = sm$mcount,
        allcount = sm$allcount, mpercent = sm$mpercent,
        library_rt = sm$library_rt, query_rt = sm$query_rt,
        library_precursor_mz = sm$library_precursor_mz,
        query_precursor_mz = sm$query_precursor_mz,
        library_precursor_ion_purity = sm$library_precursor_ion_purity,
        query_precursor_ion_purity = sm$query_precursor_ion_purity,
        library_accession = sm$library_accession,
        library_precursor_type = sm$library_precursor_type,
        library_entry_name = sm$library_entry_name, inchikey = sm$inchikey,
        library_source_name = sm$library_source_name,
        library_compound_name = sm$library_compound_name,
        mid = sm$mid, stringsAsFactors = FALSE)
    mc <- src$metab_compound %||% data.frame(inchikey_id = character())
    cmp <- .mzs_compound_rows(data.frame(
        inchikey = mc$inchikey_id,
        chemical_name = .mzs_na_empty(mc$name),
        chemical_formula = .mzs_na_empty(mc$molecular_formula),
        smiles = .mzs_na_empty(mc$smiles),
        theoretical_neutral_mass = suppressWarnings(as.numeric(mc$exact_mass)),
        average_molecular_weight = suppressWarnings(as.numeric(
            mc$molecular_weight)),
        compound_class = .mzs_na_empty(mc$compound_class),
        pubchem = .mzs_na_empty(mc$pubchem_id),
        chemspider = .mzs_na_empty(mc$chemspider_id),
        stringsAsFactors = FALSE))
    matches$compound_id_ <- cmp$id(matches$inchikey)
    matches$compound_match_level <- .mzs_match_level(matches$inchikey)
    attempted <- unique(matches[, c("qkey", "query_source", "query_run_id",
                                    "query_spectrum_id_",
                                    "query_native_id")])
    ranked <- .mzs_rank_matches(matches, attempted)
    rows <- .mzs_evidence_rows(ranked$matches, ranked$coverage, 1L)
    ev_mid <- ranked$matches$mid[order(ranked$matches$qkey,
                                       ranked$matches$rank,
                                       ranked$matches$order)]
    syn <- NULL
    if (!is.null(cmp$compound) && !is.null(mc$other_names)) {
        on <- .mzs_na_empty(mc$other_names[match(cmp$raw,
                                                 .mzs_na_empty(mc$inchikey_id))])
        if (any(!is.na(on)))
            syn <- data.frame(compound_id_ = cmp$compound$compound_id_[!is.na(on)],
                              synonym = on[!is.na(on)],
                              stringsAsFactors = FALSE)
    }
    tables <- list(evidence = rows$evidence,
                   evidence_score = rows$evidence_score,
                   coverage = rows$coverage,
                   software = data.frame(
                       software_id_ = 1L, name = "msPurity",
                       version = sub("^.* ", "", src$source$parsing_software[1L] %||%
                                         NA_character_),
                       parameters = NA_character_, stringsAsFactors = FALSE))
    if (!is.null(cmp$compound))
        tables$compound <- cmp$compound
    if (!is.null(cmp$xref))
        tables$compound_xref <- cmp$xref
    if (!is.null(syn))
        tables$compound_synonym <- syn
    list(tables = tables,
         sources = lapply(unique(key), function(k)
             .mzs_external_source(k, role = "library",
                                  version = lib[match(k, key)])),
         source_ids = .mzs_bind_source_ids(
             .mzs_source_ids("evidence", rows$evidence$evidence_id_,
                             "sm_matches", "mid",
                             src$lexical$mid[match(ev_mid, sm$mid)],
                             "positional"),
             if (!is.null(cmp$compound))
                 .mzs_source_ids("compound", cmp$compound$compound_id_,
                                 "metab_compound", "inchikey_id", cmp$raw,
                                 "name")),
         losses = list(
             list(construct = "spectral matching scope",
                  disposition = "not_present_in_source",
                  target_ref = "coverage",
                  reason = paste("The database records the matches made, not",
                                 "which spectra matching was attempted on:",
                                 "coverage lists the matched queries only.")),
             if (!is.null(src$l_s_peak_meta) && nrow(src$l_s_peak_meta))
                 list(construct = "l_s_peak_meta",
                      disposition = "no_such_concept_in_target",
                      affected_rows = nrow(src$l_s_peak_meta),
                      target_ref = "evidence.reference_native_id",
                      reason = paste("A copy of the matched library entries;",
                                     "they are referenced by accession, and",
                                     "their metadata belongs to the library.")),
             if (!is.null(src$xcms_match) && nrow(src$xcms_match))
                 list(construct = "xcms_match",
                      disposition = "dropped_by_configuration",
                      affected_rows = nrow(src$xcms_match),
                      reason = paste("A join of sm_matches with c_peak_groups,",
                                     "recomputable from evidence and",
                                     "feature.")),
             if (nrow(mc) && any(!is.na(.mzs_na_empty(mc$created_at))))
                 list(construct = "metab_compound.created_at, updated_at",
                      disposition = "no_such_concept_in_target",
                      affected_rows = nrow(mc),
                      reason = "Row timestamps of the database."),
             list(construct = "sm_matches.mid",
                  disposition = "preserved_opaque",
                  affected_rows = nrow(sm), target_ref = "source_identifier",
                  reason = "Recorded as a positional source identifier (R-017).")))
}
