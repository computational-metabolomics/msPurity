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

# Spectra slots of purityD.
#
# The averaged peak list of each file is one spectrum in the avSpectra slot,
# for each stage ("orig", as averaged, and "processed", after filterp,
# subtract and dimsPredictPurity), with every other column as a peak
# variable. The files are the sample data of the experiment slot. Methods
# rebuild the avPeaks list from avSpectra, work on it as before, and store
# the result back; avPeaks itself is kept while
# options(msPurity.legacySlots = TRUE), the default.

# avSpectra from an avPeaks list.
.pd_peaks_to_spectra <- function(avPeaks) {
    rows <- list()
    mats <- list()
    for (stage in c("orig", "processed")) {
        pl <- avPeaks[[stage]]
        for (k in seq_along(pl)) {
            df <- pl[[k]]
            empty <- !is.data.frame(df) || !ncol(df)
            cols <- if (empty) character() else colnames(df)
            rows[[length(rows) + 1L]] <- data.frame(
                stage = stage, file = k,
                name = if (is.null(names(pl))) NA_character_ else names(pl)[k],
                col_names = paste(cols, collapse = ","),
                col_types = if (empty) "" else paste(vapply(
                    df, function(x) class(x)[1], ""), collapse = ","),
                row_names = if (empty) NA_character_
                            else .pa_encode_rownames(df),
                stringsAsFactors = FALSE)
            mats[length(mats) + 1L] <- list(if (empty || !nrow(df)) NULL else {
                pcols <- c("mz", "intensity", setdiff(cols, c("mz", "i")))
                df$intensity <- df$i
                m <- vapply(pcols, function(cn) as.numeric(df[[cn]]),
                            numeric(nrow(df)))
                matrix(m, ncol = length(pcols), dimnames = list(NULL, pcols))
            })
        }
    }
    if (!length(rows))
        return(Spectra::Spectra(Spectra::MsBackendMemory()))
    # Peak variables differ between stages; every spectrum gets all of them,
    # NA where its stage has none.
    pvars <- unique(unlist(lapply(mats, colnames)))
    pvars <- c("mz", "intensity", setdiff(pvars, c("mz", "intensity")))
    mats <- lapply(mats, function(m) {
        if (is.null(m)) return(NULL)
        miss <- setdiff(pvars, colnames(m))
        if (length(miss))
            m <- cbind(m, matrix(NA_real_, nrow(m), length(miss),
                                 dimnames = list(NULL, miss)))
        m[, pvars, drop = FALSE]
    })
    sv <- do.call(rbind, rows)
    sv$msLevel <- 1L
    sp <- .msp_memory(sv, mats, pvars)
    sp@metadata$names <- lapply(avPeaks[c("orig", "processed")], names)
    sp
}

# One avPeaks data frame from its spectrum.
.pd_frame <- function(m, cols, types, row_names) {
    if (!nzchar(cols))
        return(data.frame())
    cols <- strsplit(cols, ",", fixed = TRUE)[[1]]
    types <- strsplit(types, ",", fixed = TRUE)[[1]]
    out <- lapply(seq_along(cols), function(k) {
        # unname(): a column of a one-row matrix is named after the column.
        v <- unname(m[, if (cols[k] == "i") "intensity" else cols[k]])
        switch(types[k], integer = as.integer(v), logical = as.logical(v), v)
    })
    names(out) <- cols
    out <- as.data.frame(out, stringsAsFactors = FALSE, optional = TRUE)
    .pa_decode_rownames(out, row_names)
}

# The avPeaks list rebuilt from avSpectra.
.pd_spectra_to_peaks <- function(sp) {
    nms <- sp@metadata$names
    if (!length(sp) || is.null(nms))
        return(list())
    pk <- .msp_peak_list(sp)
    sv <- as.data.frame(Spectra::spectraData(sp, columns = c(
        "stage", "file", "col_names", "col_types", "row_names")),
        optional = TRUE)
    out <- list()
    for (stage in c("orig", "processed")) {
        i <- which(sv$stage == stage)
        if (!length(i))
            next
        i <- i[order(sv$file[i])]
        out[[stage]] <- lapply(i, function(j) .pd_frame(
            pk[[j]], sv$col_names[j], sv$col_types[j], sv$row_names[j]))
        names(out[[stage]]) <- nms[[stage]]
    }
    out
}

# The MsExperiment of a fileList data frame.
.pd_experiment <- function(fileList) {
    MsExperiment::MsExperiment(
        sampleData = S4Vectors::DataFrame(fileList, check.names = FALSE))
}

# The fileList data frame of an MsExperiment: its sample data, with file
# paths, names, sample types and classes filled in where missing.
.pd_filelist <- function(x) {
    sd <- as.data.frame(MsExperiment::sampleData(x), optional = TRUE)
    # sampleIdx is taken from the row numbers of fileList.
    rownames(sd) <- NULL
    n <- length(x)
    if (is.null(sd$filepth)) {
        sps <- MsExperiment::spectra(x)
        idx <- MsExperiment::spectraSampleIndex(x)
        sd$filepth <- vapply(seq_len(n), function(i)
            Spectra::dataOrigin(sps)[match(i, idx)], "")
    }
    if (is.null(sd$name))
        sd$name <- tools::file_path_sans_ext(basename(sd$filepth))
    if (is.null(sd$sampleType))
        sd$sampleType <- rep("sample", n)
    if (is.null(sd$class))
        sd$class <- rep(NA_character_, n)
    sd
}

# Rebuild avPeaks for a method to work on.
.pd_begin <- function(Object) {
    Object <- .pd_update(Object)
    Object@avPeaks <- .pd_spectra_to_peaks(Object@avSpectra)
    Object
}

# Store avPeaks as avSpectra, keeping avPeaks only for the legacy slots.
.pd_end <- function(Object) {
    Object@avSpectra <- .pd_peaks_to_spectra(Object@avPeaks)
    if (!.pa_legacy())
        Object@avPeaks <- list()
    Object
}

.pd_is_current <- function(Object) {
    all(vapply(c("avSpectra", "experiment"),
               function(s) methods::.hasSlot(Object, s), logical(1)))
}

.pd_update <- function(Object) {
    if (!.pd_is_current(Object))
        Object <- updateObject(Object)
    Object
}

#' Update a purityD object to the current class definition
#'
#' Objects saved before msPurity stored the averaged DIMS peak lists in a
#' \code{Spectra} object lack the \code{avSpectra} and \code{experiment}
#' slots. \code{updateObject()} builds them from \code{avPeaks} and
#' \code{fileList}. Every purityD method updates its input this way.
#'
#' @param object purityD object.
#' @param ... unused.
#' @param verbose logical; unused.
#'
#' @return An updated purityD object.
#' @aliases updateObject,purityD-method
#' @seealso \code{\link{averagedPeaks}}
#' @export
setMethod("updateObject", "purityD", function(object, ..., verbose = FALSE) {
    if (.pd_is_current(object))
        return(object)
    pd <- methods::new("purityD")
    for (s in setdiff(methods::slotNames("purityD"),
                      c("avSpectra", "experiment")))
        if (methods::.hasSlot(object, s))
            methods::slot(pd, s) <- methods::slot(object, s)
    pd@experiment <- .pd_experiment(pd@fileList)
    pd@avSpectra <- .pd_peaks_to_spectra(pd@avPeaks)
    pd
})

#' Averaged peaks of a purityD object
#'
#' The averaged peak list of each file, as the \code{avPeaks} list (a list
#' with elements \code{orig} and \code{processed}, each with one data frame
#' per file; \code{legacy = TRUE}, the default) or as a \code{Spectra}
#' object with one spectrum per file and stage, with the other columns as
#' peak variables (\code{legacy = FALSE}). The legacy \code{avPeaks} slot is
#' only filled while \code{options(msPurity.legacySlots = TRUE)}, the
#' default for this release.
#'
#' @param pd purityD object.
#' @param legacy logical; TRUE for the avPeaks list, FALSE for \code{Spectra}.
#'
#' @return A list or a \code{Spectra} object.
#' @aliases averagedPeaks averagedPeaks,purityD-method
#' @export averagedPeaks
setMethod("averagedPeaks", "purityD", function(pd, legacy = TRUE) {
    pd <- .pd_update(pd)
    if (legacy) .pd_spectra_to_peaks(pd@avSpectra) else pd@avSpectra
})
