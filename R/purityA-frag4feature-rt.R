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

# Retention time mapping for frag4feature(). MS/MS times are raw, read from
# each mzML file; XCMS peak and feature windows are on the adjusted scale once
# retention time correction has been applied. Every MS/MS time is mapped to
# the adjusted scale of its file before it is compared with a window.

# Columns named inside dplyr::join_by().
utils::globalVariables(c("between", "rt", "mz", ".lo", ".hi", ".mzlo",
                         ".mzhi", "rtPred", "mz1"))

#' File name of a path, splitting on both "/" and "\".
#'
#' Objects saved on Windows keep backslash paths, which basename() does not
#' split on other platforms.
#'
#' @noRd
.f4f_basename <- function(x) sub("^.*[/\\\\]", "", as.character(x))

#' Map times between two scales of one file.
#'
#' Linear interpolation between the knots `from` -> `to`. A time equal to a
#' knot maps exactly to its partner. Outside the knots, the shift at the
#' nearest end is applied, as xcms does when it adjusts spectra beyond the
#' fitted range. Swapping `from` and `to` gives the inverse map.
#'
#' @noRd
.f4f_interp <- function(x, from, to) {
    out <- rep(NA_real_, length(x))
    ok <- !is.na(x)
    keep <- !is.na(from) & !is.na(to)
    from <- from[keep]
    to <- to[keep]
    if (!any(ok))
        return(out)
    if (!length(from))
        stop("A retention time map has no knots.")
    o <- order(from)
    from <- from[o]
    to <- to[o]
    xs <- x[ok]
    if (length(unique(from)) == 1L) {
        out[ok] <- xs + (mean(to) - from[1L])
        return(out)
    }
    y <- stats::approx(from, to, xout = xs, ties = list("ordered", mean),
                       rule = 2)$y
    n <- length(from)
    lo <- xs < from[1L]
    hi <- xs > from[n]
    y[lo] <- xs[lo] + (to[1L] - from[1L])
    y[hi] <- xs[hi] + (to[n] - from[n])
    out[ok] <- y
    out
}

#' Raw and adjusted scan times of each XCMS sample.
#'
#' Returns a list with one `list(raw, adjusted)` per sample, or NULL when no
#' retention time correction was applied. When the object holds the MS/MS
#' spectra their times are among the knots, so mapping an MS/MS time returns
#' the adjusted time xcms assigned to it.
#'
#' @noRd
.f4f_rt_pairs <- function(xcmsObj) {
    if (.xcms_is_modern(xcmsObj)) {
        if (!xcms::hasAdjustedRtime(xcmsObj))
            return(NULL)
        raw <- .xcms_rtime_by_sample(xcmsObj, adjusted = FALSE)
        adj <- .xcms_rtime_by_sample(xcmsObj, adjusted = TRUE)
    } else {
        corrected <- any(vapply(xcmsObj@.processHistory, function(h)
            "Retention time correction" %in% h@type, logical(1)))
        if (!corrected)
            return(NULL)
        raw <- xcmsObj@rt$raw
        adj <- xcmsObj@rt$corrected
    }
    mapply(function(r, a) list(raw = unname(r), adjusted = unname(a)),
           raw, adj, SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

#' Position of each file of `pa` among the samples of `xcmsObj`.
#'
#' Files are matched by name, so the MS/MS files may be any subset of the
#' XCMS samples, in any order. NA marks a file the object does not hold.
#'
#' @noRd
.f4f_file_index <- function(pa_files, xcms_files) {
    xb <- .f4f_basename(xcms_files)
    dup <- unique(xb[duplicated(xb)])
    if (length(dup))
        stop("The XCMS object has several files named ",
             paste(dup, collapse = ", "),
             ", so they cannot be matched to the purityA files by name.")
    match(.f4f_basename(pa_files), xb)
}

#' Build the retention time map of each file of `pa`.
#'
#' The status of each file is one of
#' * "aligned": its adjusted times are in the XCMS object;
#' * "unaligned": no retention time correction was applied, times are used
#'   unchanged;
#' * "external": the mapping comes from `ms2RtMap`;
#' * "estimated": the mapping was estimated from the data
#'   (see .f4f_estimate_maps());
#' * "missing": the features are aligned but nothing maps this file, so its
#'   raw times are used;
#' * "not_converted": convert2RawRT = FALSE, raw times are used.
#'
#' @return list(fun = list of functions, status = character), by pa fileid.
#' @noRd
.f4f_rt_maps <- function(pa_files, sample_idx, pairs, convert2RawRT,
                         ms2RtMap = NULL) {
    n <- length(pa_files)
    identity_fun <- function(x) x
    fun <- rep(list(identity_fun), n)
    status <- rep(if (!convert2RawRT) "not_converted"
                  else if (is.null(pairs)) "unaligned" else "missing", n)
    if (!convert2RawRT || is.null(pairs))
        return(list(fun = fun, status = status))
    for (i in which(!is.na(sample_idx))) {
        p <- pairs[[sample_idx[i]]]
        fun[[i]] <- .f4f_map_fun(p$raw, p$adjusted)
        status[i] <- "aligned"
    }
    if (!is.null(ms2RtMap)) {
        ext <- .f4f_external_maps(ms2RtMap, pa_files)
        for (i in which(status == "missing" & !vapply(ext, is.null, logical(1)))) {
            fun[[i]] <- ext[[i]]
            status[i] <- "external"
        }
    }
    list(fun = fun, status = status)
}

.f4f_map_fun <- function(from, to) {
    force(from)
    force(to)
    function(x) .f4f_interp(x, from, to)
}

#' Mapping functions for the pa files from a user-supplied `ms2RtMap`.
#'
#' `ms2RtMap` is an XCMS object holding the files with adjusted times, a list
#' of functions from raw to adjusted time, or a numeric vector of constant
#' offsets (adjusted = raw + offset). Lists and vectors are named by file.
#'
#' @return A list by pa file, NULL where the map has no entry.
#' @noRd
.f4f_external_maps <- function(ms2RtMap, pa_files) {
    pb <- .f4f_basename(pa_files)
    out <- vector("list", length(pb))
    if (.xcms_is_modern(ms2RtMap) || is(ms2RtMap, "xcmsSet") ||
            is(ms2RtMap, "xsAnnotate")) {
        if (is(ms2RtMap, "xsAnnotate"))
            ms2RtMap <- ms2RtMap@xcmsSet
        pairs <- .f4f_rt_pairs(ms2RtMap)
        if (is.null(pairs))
            stop("ms2RtMap has no adjusted retention times.")
        idx <- .f4f_file_index(pb, .xcms_files(ms2RtMap))
        for (i in which(!is.na(idx)))
            out[[i]] <- .f4f_map_fun(pairs[[idx[i]]]$raw,
                                     pairs[[idx[i]]]$adjusted)
        return(out)
    }
    nms <- .f4f_basename(names(ms2RtMap))
    if (is.null(names(ms2RtMap)) || any(!nzchar(nms)))
        stop("ms2RtMap must be an XCMS object, or a list or vector named ",
             "by file.")
    if (is.numeric(ms2RtMap)) {
        for (i in which(pb %in% nms)) {
            off <- unname(ms2RtMap[match(pb[i], nms)])
            out[[i]] <- local({
                o <- off
                function(x) x + o
            })
        }
        return(out)
    }
    if (is.list(ms2RtMap) && all(vapply(ms2RtMap, is.function, logical(1)))) {
        for (i in which(pb %in% nms))
            out[[i]] <- ms2RtMap[[match(pb[i], nms)]]
        return(out)
    }
    stop("ms2RtMap must be an XCMS object, a list of functions or a numeric ",
         "vector of offsets.")
}

#' Estimate a raw-to-adjusted map for files with no mapping.
#'
#' Landmarks are MS/MS scans that match exactly one feature by m/z and by
#' time, with purity at or above `minPurity`. A first pass looks `window`
#' seconds either side of each feature rtmed and takes the median offset
#' rtmed - rt, after trimming landmarks more than 3 MADs (at least 2 s) from
#' it. Two more passes look 15 and then 8 seconds around the current estimate
#' and refit. With 10 or more landmarks, a robust loess of the offset against
#' time replaces the median offset when it lowers the spread of the residuals
#' by a quarter or more, so a drift is followed and a constant offset stays
#' constant. Outside the landmarks the offset at the nearest end is used.
#' With fewer than 3 landmarks in the first pass, no map is returned.
#'
#' @return list(fun, method, landmarks) for one file.
#' @noRd
.f4f_estimate_map <- function(scans, features, ppm, window, minPurity) {
    s <- scans[!is.na(scans$mz1) & !is.na(scans$rtRaw), , drop = FALSE]
    s <- s[is.na(s$inPurity) | s$inPurity >= minPurity, , drop = FALSE]
    fun <- function(x) x
    method <- "none"
    n <- 0L
    windows <- unique(c(window, 15, 8))
    for (w in windows) {
        m <- .f4f_landmarks(s, features, fun, ppm, w)
        if (nrow(m)) {
            res <- m$rtmed - m$rtPred
            k <- abs(res - stats::median(res)) <=
                max(3 * stats::mad(res), 2)
            m <- m[k, , drop = FALSE]
        }
        if (nrow(m) < 3L)
            break
        n <- nrow(m)
        off <- m$rtmed - m$rtRaw
        rng <- range(m$rtRaw)
        med <- stats::median(off)
        fun <- local({
            o <- med
            function(x) x + o
        })
        method <- "median"
        if (w != window && n >= 10L) {
            d <- data.frame(off = off, rtRaw = m$rtRaw)
            fit <- tryCatch(stats::loess(off ~ rtRaw, data = d, span = 0.75,
                                         degree = 1, family = "symmetric"),
                            error = function(e) NULL)
            spread <- stats::mad(off - med)
            if (!is.null(fit) && spread > 0 &&
                    stats::mad(stats::residuals(fit)) <= 0.75 * spread) {
                fun <- local({
                    f <- fit
                    r <- rng
                    function(x) x + unname(stats::predict(
                        f, data.frame(rtRaw = pmin(pmax(x, r[1L]), r[2L]))))
                })
                method <- "loess"
            }
        }
    }
    list(fun = fun, method = method, landmarks = n)
}

#' Scans matching exactly one feature by m/z and by mapped time.
#'
#' @noRd
.f4f_landmarks <- function(s, features, fun, ppm, window) {
    if (!nrow(s) || !nrow(features))
        return(data.frame(pid = integer(), rtRaw = numeric(),
                          rtPred = numeric(), rtmed = numeric()))
    s$rtPred <- fun(s$rtRaw)
    f <- features[, c("grpid", "mzmed", "rtmed")]
    f$.lo <- f$rtmed - window
    f$.hi <- f$rtmed + window
    f$.mzlo <- f$mzmed * (1 - ppm * 1e-6)
    f$.mzhi <- f$mzmed * (1 + ppm * 1e-6)
    m <- dplyr::inner_join(
        s[, c("pid", "mz1", "rtRaw", "rtPred")], f,
        by = dplyr::join_by(between(rtPred, .lo, .hi),
                            between(mz1, .mzlo, .mzhi)),
        relationship = "many-to-many")
    m <- m[check_ppm(m$mz1, m$mzmed) < ppm, , drop = FALSE]
    as.data.frame(m[!m$pid %in% m$pid[duplicated(m$pid)], , drop = FALSE])
}
