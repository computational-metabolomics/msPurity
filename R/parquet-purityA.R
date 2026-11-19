# msPurity results -> a Parquet results dataset.
#
# createDatabase(format = "parquet") builds msPurity's table set in memory,
# as the SQLite writer would store it, and adds what only the live purityA
# object holds (parquet-msn.R); convertSqliteToParquet() reads the table set
# from an existing database. Both map it the same way:
#
#   fileinfo                      -> runs of the study (or external sources)
#   averaged spectra              -> native runs: av_intra_<run> per source
#                                    run, av_inter, av_all
#   their contributing scans      -> merge_member
#
# Scans are referenced, never copied.

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

.MZS_PRECURSOR_DEFINITION_PARAM <- paste(
    "An average computed with AverageFragSpectraParam has the median",
    "precursor m/z and retention time of the scans that contributed to it,",
    "of their last stage for MS3/MS4 averages, whose earlier stages are in",
    "x_mspurity_precursor_path; its spectrum_feature link gives the MS2",
    "precursor of that path.")

.MZS_PRECURSOR_DEFINITION <- paste(
    "An averaged spectrum's precursor m/z and retention time are those of",
    "its feature (xcms featureDefinitions() mzmed and rtmed, or",
    "peakTable() mz and rt), as createDatabase() records them.")

# Metadata column types of an averaged spectrum. `rtime` is the retention
# time in seconds, which is stored as `time` in minutes.
.MZS_AV_META_TYPES <- c(ms_level = "i32", rtime = "f64", scan_polarity = "i32",
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
    ## Checksums identify files that are not in a study.
    checksum <- rep(NA_character_, nrow(fi))
    if (is.null(study))
        checksum[readable] <- vapply(fi$filepth[readable], .mzs_sha256,
                                     character(1))
    mapped <- function(i) {
        p <- fi$filepth[i]
        if (!is.null(fileMap) && p %in% names(fileMap))
            unname(fileMap[[p]]) else NA_character_
    }

    if (!is.null(study)) {
        sm <- MsBackendParquet::readManifest(study)
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

#' Scan references: source, run, study spectrum_id_, native id, and how the
#' scan was matched to the study (`link_basis`).
#'
#' Without a study, scans are external references by native id, which the
#' live object records or the mzML files give. With one, a scan is its study
#' spectrum: the input's own key when the input is that study (`dataset`),
#' else matched by native id, then by scan number at the same position and
#' MS level, all read from the study's metadata.
#'
#' @param scans data.frame with pid, fileid, seqNum, acquisitionNum and,
#'     optionally, ms_level.
#'
#' @param ann live route only: `.mzs_scan_annotations()` of the scans.
#'
#' @param fast whether the live object's spectra are the study's own.
#'
#' @noRd
.mzs_scan_refs <- function(files, scans, ann = NULL, fast = FALSE) {
    k <- match(scans$fileid, files$fileid)
    out <- data.frame(source = files$source[k], run_id = files$run_id[k],
                      spectrum_id_ = NA_real_, native_id = NA_character_,
                      link_basis = NA_character_, stringsAsFactors = FALSE)
    if (!is.null(ann))
        out$native_id <- ann$native_id
    study <- attr(files, "study")
    if (is.null(study)) {
        need <- is.na(out$native_id)
        for (j in which(files$readable)) {
            i <- which(scans$fileid == files$fileid[j] & need)
            if (!length(i))
                next
            h <- tryCatch(.msp_header(.msp_read(files$filepth[j])),
                          error = function(e) NULL)
            if (is.null(h) || is.null(h$spectrumId))
                next
            out$native_id[i] <- as.character(
                h$spectrumId[match(scans$seqNum[i], h$seqNum)])
        }
        return(out)
    }
    level <- if (is.null(scans$ms_level)) rep(2L, nrow(scans)) else
        scans$ms_level
    if (fast && !is.null(ann) && !anyNA(ann$dataset_spectrum_id_)) {
        if (any(ann$dataset_run_id != out$run_id, na.rm = TRUE))
            .mzs_abort("semantic", "The spectra of 'pa' come from runs of the ",
                       "study other than those its files map to.")
        out$spectrum_id_ <- ann$dataset_spectrum_id_
        out$link_basis <- "dataset"
        return(out)
    }
    sm <- attr(files, "study_manifest")
    meta <- .mzs_read_spectra_vars(study, sm, c("spectrumId",
        "acquisitionNum", "scanIndex", "msLevel"), unique(out$run_id))
    key <- function(r, v) ifelse(is.na(v), NA_character_,
                                 paste(r, v, sep = "\r"))
    unique_key <- function(x) {
        x[duplicated(x) | duplicated(x, fromLast = TRUE)] <- NA_character_
        x
    }
    by_id <- match(key(out$run_id, out$native_id),
                   key(meta$run_id, meta$spectrumId), incomparables = NA)
    by_scan <- match(key(out$run_id, scans$acquisitionNum),
                     unique_key(key(meta$run_id, meta$acquisitionNum)),
                     incomparables = NA)
    ## A scan number names the spectrum only at the same file position.
    moved <- !is.na(by_scan) & !is.na(scans$seqNum) &
        !is.na(meta$scanIndex[by_scan]) &
        meta$scanIndex[by_scan] != scans$seqNum
    by_scan[moved] <- NA_integer_
    hit <- ifelse(is.na(by_id), by_scan, by_id)
    out$link_basis <- ifelse(!is.na(by_id), "native_id",
                             ifelse(!is.na(by_scan), "scan_number",
                                    NA_character_))
    miss <- !is.na(out$run_id) & is.na(hit)
    if (any(miss))
        .mzs_abort("semantic", sum(miss), " scan(s) are not in the ",
                   "study run they were recorded against, e.g. run ",
                   out$run_id[miss][1L], " scan ",
                   scans$acquisitionNum[miss][1L], ".",
                   data = list(run_id = out$run_id[miss][1L]))
    wrong <- !is.na(hit) & !is.na(meta$msLevel[hit]) &
        meta$msLevel[hit] != level
    if (any(wrong))
        .mzs_abort("semantic", sum(wrong), " scan(s) match a study ",
                   "spectrum of another MS level, e.g. run ",
                   out$run_id[wrong][1L], " scan ",
                   scans$acquisitionNum[wrong][1L], ".",
                   data = list(run_id = out$run_id[wrong][1L]))
    out$spectrum_id_ <- as.numeric(meta$spectrum_id_[hit])
    nat <- is.na(out$native_id)
    out$native_id[nat] <- as.character(meta$spectrumId[hit[nat]])
    out
}

#' Parent and root references of each scan.
#'
#' The parent of an MS3 or MS4 scan, and the MS2 root of its path, are the
#' scans the live object's lineage resolved, so they carry those scans'
#' references. The parent of an MS2 scan is its MS1 survey scan, which
#' msPurity does not retain: it is the study's MS1 spectrum whose
#' acquisitionNum is the scan's precursor scan number, and unknown without
#' a study.
#'
#' @return data.frame of the `parent_` and `root_` reference columns, with
#'     `parent_row` and `root_row`, the positions of the parent and root
#'     among `scans`.
#'
#' @noRd
.mzs_parent_refs <- function(files, scans, refs, ann = NULL) {
    n <- nrow(scans)
    na_chr <- rep(NA_character_, n)
    out <- data.frame(parent_source = na_chr, parent_run_id = na_chr,
                      parent_spectrum_id_ = rep(NA_real_, n),
                      parent_native_id = na_chr, root_source = na_chr,
                      root_run_id = na_chr,
                      root_spectrum_id_ = rep(NA_real_, n),
                      root_native_id = na_chr,
                      parent_row = rep(NA_integer_, n),
                      root_row = rep(NA_integer_, n),
                      stringsAsFactors = FALSE)
    level <- if (is.null(scans$ms_level)) rep(2L, n) else scans$ms_level
    ref_cols <- c("source", "run_id", "spectrum_id_", "native_id")
    if (!is.null(ann)) {
        pr <- match(ann$parent_pid, scans$pid)
        rr <- match(ann$root_pid, scans$pid)
        msn <- which(level > 2L & !is.na(pr))
        for (c in ref_cols)
            out[[paste0("parent_", c)]][msn] <- refs[[c]][pr[msn]]
        out$parent_row[msn] <- pr[msn]
        rooted <- which(level > 2L & !is.na(rr))
        for (c in ref_cols)
            out[[paste0("root_", c)]][rooted] <- refs[[c]][rr[rooted]]
        out$root_row[rooted] <- rr[rooted]
    }
    study <- attr(files, "study")
    ms2 <- which(level == 2L & !is.na(scans$precursorScanNum))
    if (!is.null(study) && length(ms2)) {
        meta <- .mzs_read_spectra_vars(study, attr(files, "study_manifest"),
            c("spectrumId", "acquisitionNum", "msLevel"),
            unique(refs$run_id[ms2]))
        meta <- meta[meta$msLevel %in% 1L & !is.na(meta$acquisitionNum), ,
                     drop = FALSE]
        mk <- paste(meta$run_id, meta$acquisitionNum, sep = "\r")
        mk[duplicated(mk) | duplicated(mk, fromLast = TRUE)] <- NA_character_
        hit <- match(paste(refs$run_id[ms2], scans$precursorScanNum[ms2],
                           sep = "\r"), mk, incomparables = NA)
        found <- !is.na(hit)
        i <- ms2[found]
        out$parent_source[i] <- refs$source[i]
        out$parent_run_id[i] <- refs$run_id[i]
        out$parent_spectrum_id_[i] <- meta$spectrum_id_[hit[found]]
        out$parent_native_id[i] <- as.character(meta$spectrumId[hit[found]])
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
    averaged <- c("intra", "inter", "all")
    if (any(!is.na(type) & !type %in% c("scan", averaged)))
        .mzs_abort("semantic", "Unsupported stored spectrum type.")
    pk <- src$s_peaks[src$s_peaks$type %in% averaged, , drop = FALSE]
    expected <- type[match(pk$pid, m$pid)]
    if (anyNA(expected) || any(pk$type != expected))
        .mzs_abort("semantic", "Averaged peaks have missing or contradictory ",
            "spectrum metadata.")
    av <- m[!is.na(type) & type %in% averaged, ,
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

#' Peak rows with unambiguous identity and original stored average type.
#' @noRd
.mzs_average_peaks <- function(src, av) {
    pk <- src$s_peaks[src$s_peaks$pid %in% av$pid, , drop = FALSE]
    if (!nrow(pk)) return(pk)
    if (anyDuplicated(av$pid) || !"type" %in% names(pk))
        .mzs_abort("semantic", "Average peak membership is ambiguous.")
    if (anyNA(pk$type) || any(!nzchar(trimws(as.character(pk$type)))))
        .mzs_abort("semantic",
            "Average peak membership has missing type information.")
    pk <- pk[pk$type != "scan", , drop = FALSE]
    expected <- av$spectrum_type[match(pk$pid, av$pid)]
    if (anyNA(expected) || any(!pk$type %in% c("all", "intra", "inter")) ||
            any(pk$type != expected))
        .mzs_abort("semantic",
            "Averaged peaks disagree with their stored spectrum type.")
    pk
}

#' The peaks of the averaged spectra, as `.mzs_derived_runs()` writes them:
#' those of s_peaks and those of averages computed with
#' AverageFragSpectraParam, by spectrum, m/z and cluster.
#'
#' @noRd
.mzs_av_peaks <- function(src, av, live = NULL) {
    legacy <- if (is.null(av$configured)) av else
        av[!av$configured %in% TRUE, , drop = FALSE]
    pk <- .mzs_rbind_fill(.mzs_average_peaks(src, legacy), live$peaks)
    if (is.null(pk))
        return(.mzs_average_peaks(src, legacy))
    pk[order(pk$pid, pk$mz, pk$cl), , drop = FALSE]
}

#' The native runs holding the averaged spectra, for `.mzs_commit_new()`.
#'
#' @param pk `.mzs_av_peaks()`.
#'
#' @noRd
.mzs_derived_runs <- function(av, pk) {
    if (!nrow(av))
        return(list(runs = list(), results_runs = list(),
                    peak_columns = list()))
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
            ms_level = if (is.null(a$ms_level)) rep(2L, nrow(a))
                       else as.integer(a$ms_level),
            rtime = a$retention_time,
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
        ## How averages computed with AverageFragSpectraParam were made.
        if (any(av$configured %in% TRUE))
            for (c in names(.MZS_AV_CONFIGURED_TYPES))
                meta[[c]] <- a[[sub("^x_mspurity_", "", c)]]
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
#' Where the averaging provenance is recorded (`gk`), the members are the
#' inputs it records, selected and contributing. Otherwise they are
#' reconstructed: all-file averages draw on the feature's passing scans,
#' intra-file averages on those of their file, inter-file averages on the
#' feature's intra-file averages; a scan linked twice is listed twice, and
#' membership is complete when it matches the recorded contributor count.
#'
#' @param members data.frame of contributions: grpid, pid, fileid, one row
#'     per scan-to-feature link.
#'
#' @param gk `.mzs_av_graph()`, or `NULL`.
#'
#' @noRd
.mzs_merge_member <- function(src, av, members, scans, scan_refs,
                              gk = NULL) {
    if (!nrow(av))
        return(list(table = NULL, incomplete = 0L))
    recorded <- if (is.null(gk)) NULL
                else .mzs_graph_merge_member(av, gk, scans, scan_refs)
    keep <- !(av$configured %in% TRUE) &
        !(if (is.null(gk)) FALSE else gk$use)
    av_all <- av
    av <- av[keep, , drop = FALSE]
    if (!nrow(av)) {
        if (is.null(recorded))
            return(list(table = NULL, incomplete = 0L))
        return(list(table = recorded, incomplete = 0L))
    }
    avp <- .mzs_average_peaks(src, av)
    avp <- avp[!is.na(avp$total), c("pid", "total"), drop = FALSE]
    total <- tapply(avp$total, avp$pid, max)
    out <- lapply(seq_len(nrow(av)), function(i) {
        a <- av[i, ]
        if (a$method == "inter") {
            mem <- av_all[av_all$method == "intra" & av_all$grpid == a$grpid &
                          !(av_all$configured %in% TRUE), , drop = FALSE]
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
        return(list(table = recorded, incomplete = 0L))
    mm$merged_source <- "self"
    mm$member_rank <- as.integer(stats::ave(mm$merged_spectrum_id_,
                                            mm$merged_spectrum_id_,
                                            FUN = seq_along))
    want <- suppressWarnings(as.integer(total[as.character(
        av$pid[match(mm$merged_spectrum_id_, av$spectrum_id_)])]))
    got <- as.integer(stats::ave(mm$member_rank, mm$merged_spectrum_id_,
                                 FUN = length))
    mm$members_complete <- !is.na(want) & got == want
    list(table = .mzs_rbind_fill(mm, recorded),
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

#' Resolve scan peaks by the file-position convention of msPurity writers.
#' @noRd
.mzs_scan_peak_pid <- function(src, sp) {
    m <- .mzs_scans(src)
    key <- paste(m$fileid, m$seqNum)
    if (anyDuplicated(key))
        .mzs_abort("semantic",
            "Ambiguous file-position metadata for scan peaks.")
    if (nrow(sp) && anyNA(m$seqNum))
        .mzs_abort("semantic",
            "Scan peak ownership requires file-position metadata.")
    k <- match(paste(sp$fileid, sp$scan), key)
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
    pa <- .pa_ms2_view(pa)
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
    fileid <- if ("fileid" %in% names(g)) g$fileid else g$sample
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

#' x_mspurity_scan (assessed scans with their purity and lineage) and, if
#' filterFragSpectra() flagged peaks, x_mspurity_scan_peak.
#'
#' Peaks flagged in s_peaks come first. On the live route the flags of the
#' other fragmentation spectra, those above MS2 and the feature-linked MS2
#' scans of an `allfrag = FALSE` filter, come from fragSpectra.
#'
#' @param parents `.mzs_parent_refs()` of the scans.
#'
#' @param ann,pa live route only.
#'
#' @noRd
.mzs_scan_tables <- function(src, scans, scan_refs, parents = NULL,
                             ann = NULL, pa = NULL) {
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
    scan$ms_level <- if (is.null(scans$ms_level)) rep(2L, n) else
        as.integer(scans$ms_level)
    scan$link_basis <- scan_refs$link_basis
    if (!is.null(parents)) {
        for (c in grep("^(parent|root)_", names(parents), value = TRUE))
            if (!c %in% c("parent_row", "root_row"))
                scan[[c]] <- parents[[c]]
        scan$parent_scan_annotation_id_ <- parents$parent_row
        scan$root_scan_annotation_id_ <- parents$root_row
    }
    if (!is.null(ann)) {
        for (c in c("parent_status", "path_status", "parent_reference_type",
                    "stage_purity", "purity_status", "activation_type",
                    "activation_type_basis", "scan_filter_analyzer"))
            scan[[c]] <- ann[[c]]
    }
    sp <- src$s_peaks[src$s_peaks$type %in% "scan", , drop = FALSE]
    ann_cols <- intersect(c("snr", "ra", "purity_pass_flag",
                            "intensity_pass_flag", "ra_pass_flag",
                            "snr_pass_flag"), names(sp))
    ann_cols <- ann_cols[vapply(ann_cols, function(c) any(!is.na(sp[[c]])),
                                logical(1))]
    peak <- NULL
    if (length(ann_cols) && nrow(sp)) {
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
            snr = if ("snr" %in% ann_cols) sp$snr else NA_real_,
            relative_intensity = if ("ra" %in% ann_cols) sp$ra else NA_real_,
            stringsAsFactors = FALSE)
        for (f in c("purity_pass_flag", "intensity_pass_flag",
                    "ra_pass_flag", "snr_pass_flag"))
            peak[[f]] <- if (f %in% ann_cols) sp[[f]] else NA
    }
    annotated <- !is.null(peak)
    if (!is.null(pa)) {
        done <- if (annotated) scans$pid[unique(peak$scan_annotation_id_)]
                else integer()
        lp <- .mzs_live_scan_peaks(pa, setdiff(scans$pid, done))
        if (!is.null(lp) && nrow(lp)) {
            extra <- data.frame(
                scan_annotation_id_ = match(lp$pid, scans$pid),
                peak_rank = as.integer(stats::ave(seq_len(nrow(lp)), lp$pid,
                                                  FUN = seq_along)),
                mz = lp$mz, intensity = lp$i, snr = lp$snr,
                relative_intensity = lp$ra,
                purity_pass_flag = lp$purity_pass_flag,
                intensity_pass_flag = lp$intensity_pass_flag,
                ra_pass_flag = lp$ra_pass_flag,
                snr_pass_flag = lp$snr_pass_flag,
                stringsAsFactors = FALSE)
            extra <- extra[!is.na(extra$scan_annotation_id_), , drop = FALSE]
            peak <- if (is.null(peak)) extra else rbind(peak, extra)
        }
    }
    ## Raw scan peaks the direct route counted rather than read.
    uncopied <- as.integer(attr(src$s_peaks, "raw_scan_peaks") %||% 0L)
    list(scan = scan, peak = peak, raw_peaks = nrow(sp) + uncopied,
         annotated = annotated)
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
                 affected_rows = nrow(pk) +
                     as.integer(attr(pk, "raw_scan_peaks") %||% 0L),
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
                       "dataset), so no spectrum id is assigned (R-026).")),
    "x_mspurity_scan.link_basis" = list(
        disposition = "not_present_in_source",
        reason = paste("No study dataset: the scans are external references",
                       "by native id.")),
    "x_mspurity_scan.parent_source" = list(
        disposition = "not_present_in_source",
        reason = paste("An MS2 scan's MS1 survey scan is referenced only in a",
                       "study dataset, and no MS3/MS4 scan was assessed.")),
    "x_mspurity_scan.parent_run_id" = list(
        disposition = "not_present_in_source",
        reason = paste("An MS2 scan's MS1 survey scan is referenced only in a",
                       "study dataset, and no MS3/MS4 scan was assessed.")),
    "x_mspurity_scan.parent_spectrum_id_" = list(
        disposition = "not_present_in_source",
        reason = paste("Parent scans are given study spectrum ids only in a",
                       "study dataset (R-026).")),
    "x_mspurity_scan.parent_native_id" = list(
        disposition = "not_present_in_source",
        reason = paste("An MS2 scan's MS1 survey scan is referenced only in a",
                       "study dataset, and no MS3/MS4 scan was assessed.")),
    "x_mspurity_scan.root_source" = list(
        disposition = "not_present_in_source",
        reason = "No MS3/MS4 scan was assessed; an MS2 scan is its own root."),
    "x_mspurity_scan.root_run_id" = list(
        disposition = "not_present_in_source",
        reason = "No MS3/MS4 scan was assessed; an MS2 scan is its own root."),
    "x_mspurity_scan.root_spectrum_id_" = list(
        disposition = "not_present_in_source",
        reason = paste("No MS3/MS4 scan was assessed in a study dataset; an",
                       "MS2 scan is its own root.")),
    "x_mspurity_scan.root_native_id" = list(
        disposition = "not_present_in_source",
        reason = "No MS3/MS4 scan was assessed; an MS2 scan is its own root."),
    "x_mspurity_scan.parent_scan_annotation_id_" = list(
        disposition = "not_present_in_source",
        reason = "Only MS3/MS4 scans have a parent among the assessed scans."),
    "x_mspurity_scan.root_scan_annotation_id_" = list(
        disposition = "not_present_in_source",
        reason = "Only MS3/MS4 scans have a root among the assessed scans."),
    "x_mspurity_scan.parent_status" = list(
        disposition = "not_present_in_source",
        reason = paste("The source holds no acquisition lineage: the SQLite",
                       "schema has none, and older purityA objects did not",
                       "record it.")),
    "x_mspurity_scan.path_status" = list(
        disposition = "not_present_in_source",
        reason = paste("The source holds no acquisition lineage: the SQLite",
                       "schema has none, and older purityA objects did not",
                       "record it.")),
    "x_mspurity_scan.parent_reference_type" = list(
        disposition = "not_present_in_source",
        reason = paste("The source holds no acquisition lineage: the SQLite",
                       "schema has none, and older purityA objects did not",
                       "record it.")),
    "x_mspurity_scan.stage_purity" = list(
        disposition = "not_present_in_source",
        reason = paste("No MS3/MS4 stage purity was assessed; MS2 precursor",
                       "purity is in in_purity.")),
    "x_mspurity_scan.purity_status" = list(
        disposition = "not_present_in_source",
        reason = paste("The source records no purity assessment status: the",
                       "SQLite schema has none, and older purityA objects",
                       "did not record it.")),
    "x_mspurity_scan.activation_type" = list(
        disposition = "not_present_in_source",
        reason = paste("Neither the source spectra nor their filter strings",
                       "name an activation method.")),
    "x_mspurity_scan.activation_type_basis" = list(
        disposition = "not_present_in_source",
        reason = paste("Neither the source spectra nor their filter strings",
                       "name an activation method.")),
    "x_mspurity_merge_peak.member_peak_rank" = list(
        disposition = "not_present_in_source",
        reason = paste("The averaging provenance does not record which",
                       "peak of these inputs a peak came from.")),
    "x_mspurity_scan.scan_filter_analyzer" = list(
        disposition = "not_present_in_source",
        reason = "The filter strings name no analyzer."))

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
    raw_peaks <- src$s_peaks[src$s_peaks$type %in% "scan", , drop = FALSE]
    if (is.null(live) && nrow(raw_peaks)) {
        pid <- .mzs_scan_peak_pid(src, raw_peaks)
        if (any(!is.na(pid) & (is.na(raw_peaks$pid) | pid != raw_peaks$pid)))
            warning(warningCondition(paste0(
                "Stored scan peak identities disagree with their ",
                "file positions; ",
                "identities are remapped, but historical flags and ",
                "removed peaks ",
                "cannot be repaired automatically."),
            class = "msPurity_parquet_identity_warning"))
    }
    ## The MS2 scans of the table set, then, from the live object, those
    ## above MS2, which SQLite has no place for.
    pa <- live$pa
    scans$ms_level <- rep(2L, nrow(scans))
    msn <- if (!is.null(pa)) .mzs_msn_scans(pa)
    if (!is.null(msn)) {
        msn$ms_level <- as.integer(Spectra::msLevel(pa@spectra)[
            match(msn$pid, pa@spectra$pid)])
        scans <- .mzs_rbind_fill(scans, msn)
    }
    ann <- if (!is.null(pa)) .mzs_scan_annotations(pa, scans$pid)
    ## The input's own dataset keys are the study's when purityA() read the
    ## study itself: the same dataset (by uid, or by path for one without),
    ## with runs unchanged since (checked when the study was resolved).
    read <- if (!is.null(pa)) .pa_study(pa)
    same <- function(a, b) !is.null(a) && !is.null(b) && !is.na(a) &&
        !is.na(b) && identical(a, b)
    fast <- !is.null(read) && !is.null(study) && is.null(fileMap) &&
        (same(read$uid, attr(files, "study_manifest")$uid) ||
         (is.null(attr(files, "study_manifest")$uid) &&
          identical(read$path, normalizePath(study, mustWork = FALSE))))
    scan_refs <- .mzs_scan_refs(files, scans, ann, fast)
    parents <- .mzs_parent_refs(files, scans, scan_refs, ann)
    av <- .mzs_derived_index(src, files, uid)
    live_av <- if (!is.null(pa))
        .mzs_live_averages(pa, files, uid, live$metadata)
    av <- .mzs_av_combine(av, live_av, files)
    av_peaks <- .mzs_av_peaks(src, av, live_av)
    derived <- .mzs_derived_runs(av, av_peaks)
    gk <- if (!is.null(pa)) .mzs_av_graph(pa, av, av_peaks)
    members <- if (!is.null(live)) .mzs_members_from_pa(live$pa)
               else .mzs_members_from_db(src, scans)
    mm <- .mzs_merge_member(src, av, members, scans, scan_refs, gk)
    merge_peak <- if (!is.null(gk)) .mzs_merge_peak(av, gk)
    st <- .mzs_scan_tables(src, scans, scan_refs, parents, ann, pa)
    ftabs <- .mzs_feature_tables(src, files)
    sf <- .mzs_spectrum_feature(
        src, scans, scan_refs, av, ftabs,
        match_details = if (!is.null(live)) .mzs_match_details(live$pa))
    msn_feature <- if (!is.null(pa))
        .mzs_msn_feature(pa, scans, st$scan$scan_annotation_id_, ftabs)

    tables <- c(ftabs$tables, list(spectrum_feature = sf,
                                   x_mspurity_scan = st$scan))
    if (!is.null(msn_feature))
        tables$x_mspurity_msn_feature <- msn_feature
    evd <- if (is.null(live)) .mzs_sqlite_evidence(src, files, av, scans,
                                                     scan_refs, ftabs)
    tables <- c(tables, evd$tables)
    if (!is.null(st$peak))
        tables$x_mspurity_scan_peak <- st$peak
    if (!is.null(mm$table))
        tables$merge_member <- mm$table
    if (!is.null(merge_peak))
        tables$x_mspurity_merge_peak <- merge_peak

    tables$conversion <- .mzs_conversion_row(
        route = paste0("mspurity:", route),
        source_format = if (is.null(live)) "msPurity createDatabase SQLite"
                        else "msPurity purityA",
        source_format_version = version,
        coverage = .mzs_mspurity_coverage(src, live),
        source_uri = source_uri, checksum = checksum)
    tables$tool_provenance <- .mzs_tool_provenance(
        .mzs_mspurity_steps(src, live, version))
    ## The processing records behind each step, when the object has them.
    history <- if (!is.null(pa)) .mzs_history_json(pa) else list()
    if (length(history) && !is.null(tables$tool_provenance))
        tables$tool_provenance$x_mspurity_history <- vapply(
            tables$tool_provenance$step_name, function(s)
                history[[s]] %||% NA_character_, "", USE.NAMES = FALSE)
    tables$loss_ledger <- .mzs_loss_ledger(c(
        .mzs_mspurity_losses(src, live, st, mm$incomplete),
        .mzs_msn_losses(files, scans, pa),
        .mzs_av_losses(av, gk, merge_peak, pa),
        .mzs_feature_losses(src, ftabs, live),
        evd$losses,
        .mzs_null_column_losses(tables, .MZS_NULL_REASONS)))
    tables$source_identifier <- .mzs_bind_source_ids(
        .mzs_source_ids("spectra", av$spectrum_id_[av$pid > 0],
                        "s_peak_meta", "pid",
                        src$lexical$pid[match(av$pid[av$pid > 0],
                                              src$s_peak_meta$pid)],
                        "positional"),
        if (any(av$configured %in% TRUE))
            .mzs_source_ids("spectra",
                            av$spectrum_id_[av$configured %in% TRUE],
                            "avSpectra", "spectrumId",
                            av$spectrum_key[av$configured %in% TRUE],
                            "name"),
        .mzs_source_ids("x_mspurity_scan",
                        st$scan$scan_annotation_id_[scans$ms_level == 2L],
                        "s_peak_meta", "pid",
                        src$lexical$pid[match(scans$pid[scans$ms_level == 2L],
                                              src$s_peak_meta$pid)],
                        "positional"),
        .mzs_source_ids("x_mspurity_scan",
                        st$scan$scan_annotation_id_[scans$ms_level > 2L],
                        "spectra", "pid",
                        scans$pid[scans$ms_level > 2L], "positional"),
        if (!is.null(msn_feature))
            .mzs_source_ids("x_mspurity_msn_feature",
                            msn_feature$msn_feature_id_, "grped_df", "pid",
                            scans$pid[match(msn_feature$scan_annotation_id_,
                                            st$scan$scan_annotation_id_)],
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
                      x_mspurity_instrument_type = "str",
                      .MZS_AV_CONFIGURED_TYPES),
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
                precursor_definition = .MZS_PRECURSOR_DEFINITION,
                precursor_definition_param = if (any(av$configured %in% TRUE))
                    .MZS_PRECURSOR_DEFINITION_PARAM))))
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
    study <- .pa_study_resolve(pa, study)
    .mzs_check_destination(path, overwrite, study)

    src <- .mzs_src(.pa_ms2_view(pa), xcmsObj, xset, metadata)
    message("Writing the Parquet dataset '", path, "'")
    .mzs_mspurity_convert(
        src, path, study, studyKey, fileMap, overwrite, route = "live-object",
        live = list(pa = pa, xcmsObj = xcmsObj, metadata = metadata),
        fn = "createDatabase", started = started, inputs = list(),
        parameters = list(), source_uri = NA_character_)
    invisible(normalizePath(path))
}

#' The msPurity table set of `pa` and `xcmsObj`, built in memory as
#' createDatabase()'s SQLite writer would store it.
#'
#' The option `msPurity.parquet.route = "bridge"` writes and reads back a
#' temporary SQLite database instead, so the two can be compared.
#'
#' @noRd
.mzs_src <- function(pa, xcmsObj, xset, metadata) {
    if (!is.na(xset))
        xcmsObj <- xset
    route <- getOption("msPurity.parquet.route", "direct")
    if (identical(route, "bridge")) {
        tmp <- tempfile("mzs-createDatabase-")
        dir.create(tmp)
        on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
        db <- .createDatabase_sqlite(pa = pa, xcmsObj = xcmsObj, xsa = NULL,
                                     outDir = tmp, grpPeaklist = NA,
                                     dbName = "staging.sqlite",
                                     metadata = metadata)
        prep <- if (!is.null(db)) db
    } else {
        prep <- .db_prepare(pa, .db_peaklist(xcmsObj), xcmsObj, xsa = NULL,
                            metadata = metadata, scanPeaks = "count")
    }
    if (is.null(prep))
        .mzs_abort("semantic", "createDatabase() could not match the files ",
                   "of 'pa' and 'xcmsObj'.")
    if (identical(route, "bridge")) .mzs_read_sqlite(prep)
    else .mzs_src_from_tables(prep)
}

#' Precursor match details per MS2 scan from frag4feature().
#'
#' @noRd
.mzs_match_details <- function(pa) {
    g <- .pa_ms2_view(pa)@grped_df
    need <- c("pid", "grpid", "precurMtchMZ", "precurMtchPPM", "precurMtchRT")
    if (!nrow(g) || !all(need %in% names(g)))
        return(NULL)
    out <- data.frame(pid = as.integer(g$pid),
               cid = if ("cid" %in% names(g)) as.integer(g$cid)
                     else NA_integer_,
               grpid = as.integer(as.character(g$grpid)),
               precurMtchMZ = as.numeric(g$precurMtchMZ),
               precurMtchPPM = as.numeric(g$precurMtchPPM),
               precurMtchRT = as.numeric(g$precurMtchRT))
    for (c in intersect(names(.MZS_LINK_DETAILS), names(g)))
        out[[c]] <- g[[c]]
    out
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
