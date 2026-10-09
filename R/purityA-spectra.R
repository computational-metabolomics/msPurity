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

# Spectra slots of purityA.
#
# The spectra, fragSpectra and avSpectra slots hold the results; the legacy
# slots (puritydf, grped_ms2, all_frag_scans and av_spectra) are filled from
# the same computation while options(msPurity.legacySlots) is TRUE, the
# default. Every method reads its inputs from the Spectra slots, and the
# accessors rebuild the legacy shapes from them.

# Peak variables of fragSpectra once filterFragSpectra() has run.
.PA_FRAG_PEAK_VARS <- c("mz", "intensity", "snr", "ra", "intensity_pass_flag",
                        "ra_pass_flag", "snr_pass_flag", "pass_flag", "sid")

# Columns of a filtered matrix in grped_ms2.
.PA_FRAG_LEGACY_COLS <- c("mz", "i", "snr", "ra", "purity_pass_flag",
                          "intensity_pass_flag", "ra_pass_flag",
                          "snr_pass_flag", "pass_flag")

# Whether the legacy slots are filled.
.pa_legacy <- function() isTRUE(getOption("msPurity.legacySlots", TRUE))

# An empty legacy data frame slot, as the class prototype has it.
.pa_empty_frame <- function() methods::new("data.frame")

# Row names of a data frame as one string, exactly as stored: integer ("i:")
# or character ("c:"), or NA for R's compact form of 1..n. Subsetting can
# store 1..n as an explicit integer vector, which identical() tells apart
# from the compact form.
.pa_encode_rownames <- function(df) {
    rn <- .row_names_info(df, type = 0L)
    if (is.integer(rn) && length(rn) == 2L && is.na(rn[1]))
        return(NA_character_)
    paste0(if (is.integer(rn)) "i:" else "c:", paste(rn, collapse = ","))
}

.pa_decode_rownames <- function(df, code) {
    if (is.na(code))
        return(df)
    rn <- strsplit(substring(code, 3), ",", fixed = TRUE)[[1]]
    attr(df, "row.names") <- if (startsWith(code, "i:")) as.integer(rn) else rn
    df
}

# ---------------------------------------------------------------------------
# spectra: MS/MS scans with the purity results as spectra variables.
# ---------------------------------------------------------------------------

# The spectra slot from the raw spectra of all files and puritydf. Each
# puritydf row is matched to its scan by file position and scan index.
.pa_spectra_new <- function(raw, fileid, puritydf) {
    if (!nrow(puritydf))
        return(Spectra::Spectra(Spectra::MsBackendMemory()))
    idx <- match(paste(as.character(puritydf$fileid), puritydf$seqNum),
                 paste(fileid, Spectra::scanIndex(raw)))
    sp <- raw[idx]
    sp <- .pa_set_purity_vars(sp, puritydf)
    sp
}

# Add or replace puritydf columns as spectra variables, remembering the
# column order of puritydf.
.pa_set_purity_vars <- function(sp, df) {
    for (cn in colnames(df))
        sp[[cn]] <- df[[cn]]
    sp@metadata$puritydfColumns <- colnames(df)
    sp
}

# The spectra slot of an object saved before the Spectra slots existed. The
# scans hold no peaks; they are read from fileList when needed.
.pa_spectra_from_table <- function(puritydf, fileList) {
    n <- nrow(puritydf)
    if (!n)
        return(Spectra::Spectra(Spectra::MsBackendMemory()))
    col <- function(nm, f = identity) {
        if (nm %in% colnames(puritydf)) f(puritydf[[nm]]) else rep(NA, n)
    }
    fl <- unname(as.character(fileList))
    fid <- suppressWarnings(as.integer(as.character(col("fileid"))))
    svars <- data.frame(msLevel = rep(2L, n),
                        rtime = col("retentionTime", as.numeric),
                        precursorMz = col("precursorMZ", as.numeric),
                        scanIndex = col("seqNum", as.integer),
                        precScanNum = col("precursorScanNum", as.integer),
                        dataOrigin = if (length(fl)) fl[fid] else NA_character_,
                        stringsAsFactors = FALSE)
    empty <- rep(list(matrix(numeric(0), ncol = 2,
                             dimnames = list(NULL, c("mz", "intensity")))), n)
    sp <- .msp_memory(svars, empty)
    .pa_set_purity_vars(sp, puritydf)
}

# puritydf rebuilt from the spectra slot.
.pa_purity_table <- function(sp) {
    cols <- sp@metadata$puritydfColumns
    if (!length(sp) || is.null(cols))
        return(.pa_empty_frame())
    as.data.frame(Spectra::spectraData(sp, columns = cols), optional = TRUE)
}

# Whether the spectra slot can read its own peaks.
.pa_has_peaks <- function(sp) {
    if (!length(sp) || !is(sp@backend, "MsBackendMzR"))
        return(FALSE)
    all(file.exists(unique(Spectra::dataStorage(sp))))
}

# The raw scans of the given pids, in memory. Peaks come from the spectra
# slot when its files are readable, otherwise from fileList, matched by file
# position and scan index.
.pa_raw_scans <- function(pa, pids) {
    sp <- pa@spectra
    rows <- match(pids, sp$pid)
    if (anyNA(rows))
        stop("Scans ", paste(head(pids[is.na(rows)]), collapse = ", "),
             " are not in the purityA object.")
    sp <- sp[rows]
    if (!.pa_has_peaks(sp)) {
        files <- unname(as.character(pa@fileList))
        raw <- .msp_read(files)
        o <- Spectra::dataOrigin(raw)
        fid <- match(o, unique(o))
        key <- paste(as.character(sp$fileid), sp$seqNum)
        idx <- match(key, paste(fid, Spectra::scanIndex(raw)))
        if (anyNA(idx))
            stop("Scans could not be found in the files of fileList.")
        pk <- .msp_peaks(raw[idx])
        sv <- Spectra::spectraData(sp)
        sv$mz <- NULL
        sv$intensity <- NULL
        out <- .msp_memory(as.data.frame(sv, optional = TRUE), pk)
        out@metadata <- sp@metadata
        return(out)
    }
    Spectra::setBackend(sp, Spectra::MsBackendMemory())
}

# ---------------------------------------------------------------------------
# fragSpectra: the scans used downstream, each stored once.
# ---------------------------------------------------------------------------

# fragSpectra from grped_ms2 and grped_df. Each scan takes its peaks from the
# first link that refers to it.
.pa_frag_from_grouped <- function(pa, ms2, grped_df = pa@grped_df) {
    if (!nrow(grped_df))
        return(Spectra::Spectra(Spectra::MsBackendMemory()))
    grp <- as.character(grped_df$grpid)
    pos <- stats::ave(seq_along(grp), grp, FUN = seq_along)
    mats <- lapply(seq_along(grp), function(i) {
        g <- ms2[[grp[i]]]
        if (is.null(g) || pos[i] > length(g)) NULL else g[[pos[i]]]
    })
    first <- !duplicated(grped_df$pid)
    pids <- grped_df$pid[first]
    mats <- mats[first]
    tab <- .pa_purity_table(pa@spectra)
    sv <- data.frame(msLevel = rep(2L, length(pids)))
    if (nrow(tab)) {
        k <- match(pids, tab$pid)
        sv$rtime <- as.numeric(tab$retentionTime[k])
        sv$precursorMz <- as.numeric(tab$precursorMZ[k])
        sv$scanIndex <- as.integer(tab$seqNum[k])
        sv$fileid <- tab$fileid[k]
        sv$seqNum <- tab$seqNum[k]
    }
    sv$pid <- pids
    filtered <- length(pa@filter_frag_params) > 0
    if (!filtered) {
        mats <- lapply(mats, function(m) {
            if (is.null(m)) return(NULL)
            m <- m[, 1:2, drop = FALSE]
            colnames(m) <- c("mz", "intensity")
            m
        })
        return(.msp_memory(sv, mats))
    }
    # Whether the raw scan had a single peak: setFlagMatrix() names the row
    # of such a spectrum (see .pa_frag_matrix()).
    sv$single_peak <- vapply(mats, function(m)
        !is.null(m) && !is.null(rownames(m)), logical(1))
    sv$purity_pass_flag <- vapply(seq_along(mats), function(i) {
        m <- mats[[i]]
        if (!is.null(m) && nrow(m)) return(as.logical(m[1, "purity_pass_flag"]))
        as.logical(grped_df$purity_pass_flag[first][i] %||% NA)
    }, logical(1))
    mats <- lapply(mats, function(m) {
        if (is.null(m)) return(NULL)
        m <- m[, setdiff(.PA_FRAG_LEGACY_COLS, "purity_pass_flag"), drop = FALSE]
        colnames(m)[2] <- "intensity"
        cbind(m, sid = rep(NA_real_, nrow(m)))
    })
    .msp_memory(sv, mats, .PA_FRAG_PEAK_VARS)
}

# fragSpectra for filterFragSpectra(allfrag = TRUE): every MS/MS scan, read
# from the raw data and flagged with setFlagMatrix() using its own purity
# flag. sid numbers the peaks of all scans in pid order, as all_frag_scans
# does.
.pa_frag_all <- function(pa, filter_frag_params, puritydf) {
    raw <- .pa_raw_scans(pa, puritydf$pid)
    pk <- .msp_peaks(raw)
    offset <- c(0, cumsum(vapply(pk, nrow, integer(1))))
    flag <- puritydf$purity_pass_flag
    unflagged <- filter_frag_params
    unflagged$rmp <- FALSE
    mats <- lapply(seq_along(pk), function(i) {
        m <- pk[[i]]
        if (!nrow(m)) return(NULL)
        x <- setFlagMatrix(cbind(m, flag[i]), unflagged)
        x <- cbind(x[, setdiff(.PA_FRAG_LEGACY_COLS, "purity_pass_flag"),
                     drop = FALSE],
                   sid = offset[i] + seq_len(nrow(m)))
        colnames(x)[2] <- "intensity"
        if (isTRUE(filter_frag_params$rmp))
            x <- x[x[, "pass_flag"] == 1, , drop = FALSE]
        x
    })
    sv <- data.frame(msLevel = rep(2L, length(pk)),
                     rtime = as.numeric(puritydf$retentionTime),
                     precursorMz = as.numeric(puritydf$precursorMZ),
                     scanIndex = as.integer(puritydf$seqNum))
    sv$fileid <- puritydf$fileid
    sv$seqNum <- puritydf$seqNum
    sv$pid <- puritydf$pid
    sv$purity_pass_flag <- as.logical(flag)
    sv$single_peak <- vapply(pk, nrow, integer(1)) == 1L
    .msp_memory(sv, mats, .PA_FRAG_PEAK_VARS)
}

# all_frag_scans as filterFragSpectra(allfrag = TRUE) has always written it,
# for the frozen SQLite output. Its pid column counts scans in file order,
# skipping scans without peaks, so after such a scan it names the wrong
# scan; the scan column is right. Kept as it was on purpose.
.pa_allfrag_frozen <- function(pa) {
    if (nrow(pa@all_frag_scans))
        return(pa@all_frag_scans)
    prm <- pa@filter_frag_params
    if (!isTRUE(prm$allfrag))
        return(.pa_empty_frame())
    scanpeaksFrag <- getScanPeaks(pa)
    puritydf <- purityTable(pa)
    puritydf$purity_pass_flag <- puritydf$inPurity > prm$plim
    scanpeaksFrag <- merge(puritydf[, c("pid", "inPurity", "purity_pass_flag")],
                           scanpeaksFrag, by = "pid")
    scanpeaksFrag <- scanpeaksFrag[, c("pid", "sid", "fileid", "scan", "mz",
                                       "i", "type", "purity_pass_flag"),
                                   drop = FALSE]
    plyr::ddply(scanpeaksFrag, ~pid, setFlagMatrix, filter_frag_params = prm)
}

# fragSpectra from the all_frag_scans of a saved object. The pid of each row
# is corrected from its file and scan, and the purity and overall flags are
# recomputed for that scan. Peaks that rmp = TRUE removed on the wrong flag
# cannot be restored.
.pa_frag_from_allfrag <- function(pa, afs) {
    tab <- .pa_purity_table(pa@spectra)
    key <- paste(as.character(afs$fileid), as.character(afs$scan))
    afs$pid <- tab$pid[match(key, paste(as.character(tab$fileid), tab$seqNum))]
    afs$purity_pass_flag <- tab$purity_pass_flag[match(afs$pid, tab$pid)]
    flags <- c("purity_pass_flag", "intensity_pass_flag", "ra_pass_flag",
               "snr_pass_flag")
    afs$pass_flag <- rowSums(afs[, flags]) == 4
    pids <- tab$pid
    by_pid <- split(seq_len(nrow(afs)), factor(afs$pid, levels = pids))
    mats <- lapply(by_pid, function(r) {
        if (!length(r)) return(NULL)
        x <- afs[r, , drop = FALSE]
        cbind(mz = x$mz, intensity = x$i, snr = x$snr, ra = x$ra,
              intensity_pass_flag = as.numeric(x$intensity_pass_flag),
              ra_pass_flag = as.numeric(x$ra_pass_flag),
              snr_pass_flag = as.numeric(x$snr_pass_flag),
              pass_flag = as.numeric(x$pass_flag),
              sid = as.numeric(x$sid))
    })
    sv <- data.frame(msLevel = rep(2L, length(pids)),
                     rtime = as.numeric(tab$retentionTime),
                     precursorMz = as.numeric(tab$precursorMZ),
                     scanIndex = as.integer(tab$seqNum))
    sv$fileid <- tab$fileid
    sv$seqNum <- tab$seqNum
    sv$pid <- pids
    sv$purity_pass_flag <- as.logical(tab$purity_pass_flag)
    .msp_memory(sv, unname(mats), .PA_FRAG_PEAK_VARS)
}

# One legacy grped_ms2 matrix from a fragSpectra peak matrix.
.pa_frag_matrix <- function(m, purity_flag, filtered, rmp, single = FALSE) {
    if (!filtered || !"snr" %in% colnames(m)) {
        m <- m[, c("mz", "intensity"), drop = FALSE]
        dimnames(m) <- list(NULL, c("mz", "intensity"))
        return(m)
    }
    if (rmp && !nrow(m))
        return(NULL)
    out <- cbind(m[, "mz"], m[, "intensity"], m[, "snr"], m[, "ra"],
                 rep(as.numeric(purity_flag), nrow(m)),
                 m[, "intensity_pass_flag"], m[, "ra_pass_flag"],
                 m[, "snr_pass_flag"], m[, "pass_flag"])
    # setFlagMatrix() names the row of a spectrum that had a single peak
    # after the intensity column it divides.
    one <- nrow(out) == 1 && (is.na(single) || isTRUE(single))
    dimnames(out) <- list(if (one) "intensity", .PA_FRAG_LEGACY_COLS)
    out
}

# grped_ms2 rebuilt from fragSpectra and grped_df.
.pa_grouped_legacy <- function(pa) {
    g <- pa@grped_df
    frag <- pa@fragSpectra
    if (!nrow(g) || !length(frag))
        return(list())
    filtered <- length(pa@filter_frag_params) > 0
    rmp <- filtered && isTRUE(pa@filter_frag_params$rmp)
    pk <- .msp_peak_list(frag)
    svars <- Spectra::spectraVariables(frag)
    pflag <- if ("purity_pass_flag" %in% svars)
                 frag$purity_pass_flag else rep(NA, length(frag))
    single <- if ("single_peak" %in% svars)
                  frag$single_peak else rep(NA, length(frag))
    k <- match(g$pid, frag$pid)
    mats <- lapply(seq_along(k), function(i)
        .pa_frag_matrix(pk[[k[i]]], pflag[k[i]], filtered, rmp, single[k[i]]))
    g$.row <- seq_len(nrow(g))
    out <- plyr::dlply(g, ~grpid, function(x) mats[x$.row])
    # filterFragSpectra() rebuilds the list with lapply(), which drops the
    # plyr attributes.
    if (filtered) lapply(out, identity) else out
}

# all_frag_scans rebuilt from fragSpectra.
.pa_allfrag_legacy <- function(pa) {
    if (!isTRUE(pa@filter_frag_params$allfrag) || !length(pa@fragSpectra))
        return(.pa_empty_frame())
    frag <- pa@fragSpectra
    frag <- frag[order(frag$pid)]
    pk <- .msp_peak_list(frag)
    n <- vapply(pk, nrow, integer(1))
    if (!sum(n))
        return(.pa_empty_frame())
    tab <- .pa_purity_table(pa@spectra)
    scan_levels <- unique(as.character(tab$seqNum))
    all <- do.call(rbind, pk)
    rownames(all) <- NULL
    rep_sv <- function(v) rep(v, n)
    flag <- function(nm) as.logical(unname(all[, nm]))
    col <- function(nm) unname(all[, nm])
    data.frame(sid = as.integer(col("sid")),
               pid = rep_sv(frag$pid),
               fileid = rep_sv(as.integer(as.character(frag$fileid))),
               scan = factor(rep_sv(as.character(frag$seqNum)),
                             levels = scan_levels),
               mz = col("mz"), i = col("intensity"),
               snr = col("snr"), ra = col("ra"),
               type = rep("scan", sum(n)),
               purity_pass_flag = rep_sv(frag$purity_pass_flag),
               intensity_pass_flag = flag("intensity_pass_flag"),
               ra_pass_flag = flag("ra_pass_flag"),
               snr_pass_flag = flag("snr_pass_flag"),
               pass_flag = flag("pass_flag"),
               stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# avSpectra: one spectrum per feature, averaging level and sample.
# ---------------------------------------------------------------------------

# avSpectra from the nested av_spectra list. A NULL entry for a sample is
# kept as an empty spectrum flagged is_null; a computed intra level with no
# samples is kept as one such spectrum with no sample.
.pa_av_to_spectra <- function(av, grped_df) {
    rows <- list()
    peaks <- list()
    cols <- NULL
    add <- function(g, lvl, smp, df) {
        is_null <- is.null(df)
        if (!is_null && is.null(cols))
            cols <<- colnames(df)
        # Column types can differ between spectra (cl is integer or double).
        types <- if (is_null) NA_character_
                 else paste(vapply(df, function(x) class(x)[1], ""),
                            collapse = ",")
        # Filtering with rmp = TRUE leaves gaps in the row names.
        rn <- if (is_null) NA_character_ else .pa_encode_rownames(df)
        rows[[length(rows) + 1L]] <<- data.frame(
            grpid = g, av_level = lvl, sample = smp, is_null = is_null,
            col_types = types, row_names = rn, stringsAsFactors = FALSE)
        peaks[length(peaks) + 1L] <<- list(df)
    }
    for (g in names(av)) {
        a <- av[[g]]
        if (!is.null(a$av_intra)) {
            if (!length(a$av_intra))
                add(g, "av_intra", NA_character_, NULL)
            for (s in names(a$av_intra))
                add(g, "av_intra", s, a$av_intra[[s]])
        }
        for (lvl in c("av_inter", "av_all"))
            if (!is.null(a[[lvl]]))
                add(g, lvl, NA_character_, a[[lvl]])
    }
    if (!length(rows)) {
        sp <- Spectra::Spectra(Spectra::MsBackendMemory())
        sp@metadata$groups <- names(av)
        return(sp)
    }
    sv <- do.call(rbind, rows)
    sv$msLevel <- 2L
    grp <- as.character(grped_df$grpid)
    med <- function(x, g) stats::median(x[grp == g], na.rm = TRUE)
    sv$precursorMz <- vapply(sv$grpid, med, numeric(1),
                             x = as.numeric(grped_df$precurMtchMZ))
    sv$rtime <- vapply(sv$grpid, med, numeric(1),
                       x = as.numeric(grped_df$retentionTime))
    pcols <- c("mz", "intensity", setdiff(cols, c("mz", "i")))
    mats <- lapply(peaks, function(df) {
        if (is.null(df) || !nrow(df)) return(NULL)
        df$intensity <- df$i
        m <- vapply(pcols, function(cn) as.numeric(df[[cn]]),
                    numeric(nrow(df)))
        matrix(m, ncol = length(pcols), dimnames = list(NULL, pcols))
    })
    sp <- .msp_memory(sv, mats, pcols)
    sp@metadata$groups <- names(av)
    sp@metadata$columns <- cols
    sp
}

# One averaged spectrum as its legacy data frame.
.pa_av_frame <- function(m, cols, types, row_names = NA) {
    types <- strsplit(types, ",", fixed = TRUE)[[1]]
    out <- lapply(seq_along(cols), function(k) {
        cn <- cols[k]
        # unname(): a column of a one-row matrix is named after the column.
        v <- unname(m[, if (cn == "i") "intensity" else cn])
        switch(types[k], integer = as.integer(v), logical = as.logical(v), v)
    })
    names(out) <- cols
    out <- as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
    .pa_decode_rownames(out, row_names)
}

# av_spectra rebuilt from avSpectra. The lists are built with the same plyr
# calls the averaging code uses, so they carry the same attributes; the
# SQLite writer, and possibly other code, relies on them.
.pa_av_legacy <- function(sp, grped_df = NULL) {
    groups <- sp@metadata$groups
    if (is.null(groups))
        return(list())
    cols <- sp@metadata$columns
    pk <- .msp_peak_list(sp)
    sv <- if (length(sp))
        as.data.frame(Spectra::spectraData(sp, columns = c(
            "grpid", "av_level", "sample", "is_null", "col_types",
            "row_names")), optional = TRUE)
    else data.frame(grpid = character(), av_level = character())
    frame <- function(i) {
        if (isTRUE(sv$is_null[i])) NULL
        else .pa_av_frame(pk[[i]], cols, sv$col_types[i], sv$row_names[i])
    }
    # The sample values in the type grped_df holds them in.
    sample_type <- if (!is.null(grped_df) && "sample" %in% colnames(grped_df))
                       class(grped_df$sample)[1] else "character"
    intra <- function(rows) {
        if (!length(rows))
            return(NULL)
        smp <- sv$sample[rows]
        keep <- !is.na(smp)
        rows <- rows[keep]
        smp <- methods::as(smp[keep], sample_type)
        plyr::dlply(data.frame(sample = smp, .i = rows), ~sample,
                    function(x) frame(x$.i))
    }
    level <- function(rows) if (length(rows)) frame(rows[1]) else NULL
    out <- plyr::alply(groups, 1, function(g) {
        r <- which(sv$grpid == g)
        lvl <- sv$av_level[r]
        list(av_intra = intra(r[lvl == "av_intra"]),
             av_inter = level(r[lvl == "av_inter"]),
             av_all = level(r[lvl == "av_all"]))
    })
    names(out) <- groups
    out
}

# ---------------------------------------------------------------------------
# Keeping both representations in step.
# ---------------------------------------------------------------------------

# Set the legacy slots from the Spectra slots, or empty them.
.pa_sync_legacy <- function(pa, puritydf = NULL, grped_ms2 = NULL,
                            all_frag_scans = NULL, av_spectra = NULL) {
    if (!.pa_legacy()) {
        pa@puritydf <- .pa_empty_frame()
        pa@grped_ms2 <- list()
        pa@all_frag_scans <- .pa_empty_frame()
        pa@av_spectra <- list()
        return(pa)
    }
    if (!is.null(puritydf)) pa@puritydf <- puritydf
    if (!is.null(grped_ms2)) pa@grped_ms2 <- grped_ms2
    if (!is.null(all_frag_scans)) pa@all_frag_scans <- all_frag_scans
    if (!is.null(av_spectra)) pa@av_spectra <- av_spectra
    pa
}

# Bring an object up to date before a method uses it.
.pa_update <- function(pa) {
    if (!.pa_is_current(pa))
        pa <- updateObject(pa)
    pa
}

.pa_is_current <- function(pa) {
    all(vapply(c("spectra", "fragSpectra", "avSpectra", "params"),
               function(s) methods::.hasSlot(pa, s), logical(1)))
}

#' Update a purityA object to the current class definition
#'
#' Objects saved before msPurity stored spectra in \code{Spectra} objects
#' lack the \code{spectra}, \code{fragSpectra} and \code{avSpectra} slots.
#' \code{updateObject()} builds them from the legacy slots, without access to
#' the original files: scan information comes from \code{puritydf}, the peaks
#' of grouped scans from \code{grped_ms2} (or \code{all_frag_scans}), and the
#' averaged spectra from \code{av_spectra}. The legacy slots are kept. Every
#' purityA method updates its input this way, so calling it is only needed to
#' convert saved objects explicitly.
#'
#' @param object purityA object.
#' @param ... unused.
#' @param verbose logical; unused.
#'
#' @return An updated purityA object.
#' @aliases updateObject,purityA-method
#' @seealso \code{\link{purityTable}}
#' @export
setMethod("updateObject", "purityA", function(object, ..., verbose = FALSE) {
    if (.pa_is_current(object))
        return(object)
    old <- function(s, default) {
        if (methods::.hasSlot(object, s)) methods::slot(object, s) else default
    }
    pa <- new("purityA")
    for (s in setdiff(methods::slotNames("purityA"),
                      c("spectra", "fragSpectra", "avSpectra")))
        if (methods::.hasSlot(object, s))
            methods::slot(pa, s) <- methods::slot(object, s)
    pa@spectra <- .pa_spectra_from_table(old("puritydf", data.frame()),
                                         old("fileList", character()))
    afs <- old("all_frag_scans", data.frame())
    if (isTRUE(pa@filter_frag_params$allfrag) && nrow(afs)) {
        pa@fragSpectra <- .pa_frag_from_allfrag(pa, afs)
    } else {
        pa@fragSpectra <- .pa_frag_from_grouped(pa, old("grped_ms2", list()))
    }
    av <- old("av_spectra", list())
    if (length(av))
        pa@avSpectra <- .pa_av_to_spectra(av, pa@grped_df)
    methods::validObject(pa)
    pa
})

setValidity("purityA", function(object) {
    if (!.pa_is_current(object) || !length(object@spectra))
        return(TRUE)
    msg <- character()
    pid <- object@spectra$pid
    if (anyDuplicated(pid))
        msg <- c(msg, "pid values of the spectra slot are not unique")
    g <- object@grped_df
    if (nrow(g) && "pid" %in% colnames(g) && !all(g$pid %in% pid))
        msg <- c(msg, "grped_df refers to pid values not in the spectra slot")
    if (length(msg)) msg else TRUE
})

# ---------------------------------------------------------------------------
# Accessors.
# ---------------------------------------------------------------------------

#' Results of a purityA object
#'
#' Accessors for the results stored in a purityA object. Each returns either
#' a \code{Spectra} object (\code{legacy = FALSE}) or the shape the
#' corresponding legacy slot has always had (\code{legacy = TRUE}, the
#' default). Use these instead of reading slots directly: the legacy slots
#' are only filled while \code{options(msPurity.legacySlots = TRUE)}, the
#' default for this release.
#'
#' \itemize{
#'   \item \code{purityTable()}: the MS/MS scans with their precursor ion
#'   purity, as the \code{puritydf} data frame or as \code{Spectra} with the
#'   purity results as spectra variables.
#'   \item \code{groupedSpectra()}: the MS/MS scans linked to XCMS features by
#'   \code{frag4feature()}, as the \code{grped_ms2} list (one list of peak
#'   matrices per feature) or as \code{Spectra} holding each scan once. The
#'   links themselves are in \code{pa@grped_df}, keyed by \code{pid}.
#'   \item \code{allFragSpectra()}: every MS/MS scan with the flags of
#'   \code{filterFragSpectra(allfrag = TRUE)}, as the \code{all_frag_scans}
#'   data frame or as \code{Spectra} with the flags as peak variables.
#'   \item \code{averagedSpectra()}: the averaged spectra, as the nested
#'   \code{av_spectra} list or as \code{Spectra} with one spectrum per
#'   feature, averaging level and sample.
#' }
#'
#' @param pa purityA object.
#' @param legacy logical; TRUE for the legacy shape, FALSE for \code{Spectra}.
#'
#' @return A data frame, list or \code{Spectra} object.
#' @examples
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "9_averageAllFragSpectra_with_filter_pa.rds",
#'                           package = "msPurity"))
#' head(purityTable(pa))
#' purityTable(pa, legacy = FALSE)
#' averagedSpectra(pa, legacy = FALSE)
#' @name purityA-accessors
#' @aliases purityTable groupedSpectra allFragSpectra averagedSpectra
#'   purityTable,purityA-method groupedSpectra,purityA-method
#'   allFragSpectra,purityA-method averagedSpectra,purityA-method
#' @export purityTable groupedSpectra allFragSpectra averagedSpectra
NULL

setMethod("purityTable", "purityA", function(pa, legacy = TRUE) {
    pa <- .pa_update(pa)
    if (legacy) .pa_purity_table(pa@spectra) else pa@spectra
})

setMethod("groupedSpectra", "purityA", function(pa, legacy = TRUE) {
    pa <- .pa_update(pa)
    if (legacy)
        return(.pa_grouped_legacy(pa))
    frag <- pa@fragSpectra
    frag[frag$pid %in% pa@grped_df$pid]
})

setMethod("allFragSpectra", "purityA", function(pa, legacy = TRUE) {
    pa <- .pa_update(pa)
    if (legacy)
        return(.pa_allfrag_legacy(pa))
    if (isTRUE(pa@filter_frag_params$allfrag)) pa@fragSpectra
    else Spectra::Spectra(Spectra::MsBackendMemory())
})

setMethod("averagedSpectra", "purityA", function(pa, legacy = TRUE) {
    pa <- .pa_update(pa)
    if (legacy) .pa_av_legacy(pa@avSpectra, pa@grped_df) else pa@avSpectra
})
