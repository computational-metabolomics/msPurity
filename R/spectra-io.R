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

# Raw data access through Spectra.
#
# All raw spectra are read through these helpers. They return the same
# header data frame and peak matrices that mzR::header() and mzR::peaks()
# returned, so the purity and averaging code is unchanged.

# mzR header column -> Spectra spectra variable, where the names differ.
.MSP_HEADER_RENAME <- c(seqNum = "scanIndex",
                        retentionTime = "rtime",
                        precursorScanNum = "precScanNum",
                        precursorMZ = "precursorMz",
                        isolationWindowTargetMZ = "isolationWindowTargetMz")

# Column order of mzR::header().
.MSP_HEADER_COLS <- c("seqNum", "acquisitionNum", "msLevel", "polarity",
                      "peaksCount", "totIonCurrent", "retentionTime",
                      "basePeakMZ", "basePeakIntensity", "collisionEnergy",
                      "electronBeamEnergy", "ionisationEnergy", "lowMZ",
                      "highMZ", "precursorScanNum", "precursorMZ",
                      "precursorCharge", "precursorIntensity", "mergedScan",
                      "mergedResultScanNum", "mergedResultStartScanNum",
                      "mergedResultEndScanNum", "injectionTime",
                      "filterString", "spectrumId", "centroided",
                      "ionMobilityDriftTime", "isolationWindowTargetMZ",
                      "isolationWindowLowerOffset",
                      "isolationWindowUpperOffset", "scanWindowLowerLimit",
                      "scanWindowUpperLimit")

# Read one or more mzML/mzXML files as a Spectra object. Spectra are in file
# order, then scan order, as mzR reads them.
.msp_read <- function(files) {
    files <- as.character(unname(files))
    Spectra::Spectra(files, source = Spectra::MsBackendMzR())
}

# Accept a file path or a Spectra object.
.msp_as_spectra <- function(x) {
    if (is(x, "Spectra")) x else .msp_read(x)
}

# The spectra of one file of a Spectra object, by position of the file in
# the order the files were read.
.msp_split_files <- function(sp) {
    o <- Spectra::dataOrigin(sp)
    unname(split(sp, factor(o, levels = unique(o))))
}

# The mzR::header() data frame of a Spectra object. For a single file this
# equals mzR::header() exactly, with the isolation window offsets derived
# from the window bounds Spectra stores.
.msp_header <- function(sp) {
    d <- Spectra::spectraData(sp)
    have <- colnames(d)
    out <- lapply(.MSP_HEADER_COLS, function(cn) {
        sv <- if (cn %in% names(.MSP_HEADER_RENAME)) .MSP_HEADER_RENAME[[cn]]
              else cn
        if (cn == "isolationWindowLowerOffset")
            return(d$isolationWindowTargetMz - d$isolationWindowLowerMz)
        if (cn == "isolationWindowUpperOffset")
            return(d$isolationWindowUpperMz - d$isolationWindowTargetMz)
        if (sv %in% have) d[[sv]] else rep(NA, nrow(d))
    })
    names(out) <- .MSP_HEADER_COLS
    as.data.frame(out, stringsAsFactors = FALSE)
}

# The mzR::peaks() list of peak matrices (columns mz and intensity).
.msp_peaks <- function(sp) {
    as.list(Spectra::peaksData(sp, columns = c("mz", "intensity")))
}

# An in-memory Spectra from a data frame of spectra variables and a list of
# peak matrices that all have the same columns, mz and intensity first.
# Spectra with only mz and intensity use MsBackendMemory. With more peak
# variables MsBackendDataFrame is used, which keeps each peak variable as one
# compressed list column; MsBackendMemory would keep them as a data frame per
# spectrum, about three times the memory.
.msp_memory <- function(svars, peaks, cols = c("mz", "intensity")) {
    if (!nrow(svars))
        return(Spectra::Spectra(Spectra::MsBackendMemory()))
    d <- S4Vectors::DataFrame(svars, check.names = FALSE)
    if (is.null(d$msLevel))
        d$msLevel <- rep(2L, nrow(d))
    peaks <- lapply(peaks, function(m) {
        if (is.null(m) || !nrow(m))
            return(matrix(numeric(0), ncol = length(cols),
                          dimnames = list(NULL, cols)))
        m <- m[, cols, drop = FALSE]
        storage.mode(m) <- "double"
        dimnames(m) <- list(NULL, cols)
        m
    })
    if (identical(cols, c("mz", "intensity"))) {
        be <- Spectra::backendInitialize(Spectra::MsBackendMemory(), d)
        be <- Spectra::`peaksData<-`(be, value = peaks)
        return(Spectra::Spectra(be))
    }
    for (cn in cols)
        d[[cn]] <- IRanges::NumericList(unname(lapply(peaks, function(m)
            unname(m[, cn]))), compress = TRUE)
    Spectra::Spectra(d, source = Spectra::MsBackendDataFrame(),
                     peaksVariables = cols)
}

# Peak matrices of a Spectra object. For MsBackendDataFrame without queued
# processing the matrices are built from its list columns, as its
# peaksData() does; peaksData() itself would build a data frame per
# spectrum, which is a hundred times slower.
.msp_peak_list <- function(sp, cols = Spectra::peaksVariables(sp)) {
    if (!length(sp))
        return(list())
    if (is(sp@backend, "MsBackendDataFrame") && !length(sp@processingQueue)) {
        lst <- lapply(cols, function(cn)
            as.list(sp@backend@spectraData[[cn]]))
        return(lapply(seq_along(sp), function(i) {
            m <- do.call(cbind, lapply(lst, function(l) as.numeric(l[[i]])))
            if (is.null(dim(m)))
                m <- matrix(m, ncol = length(cols))
            dimnames(m) <- list(NULL, cols)
            m
        }))
    }
    lapply(as.list(Spectra::peaksData(sp, columns = cols)), function(m) {
        if (is.data.frame(m)) {
            m <- as.matrix(m)
            dimnames(m) <- list(NULL, colnames(m))
        }
        m
    })
}

# Isolation window offsets c(lower, upper) for the purity calculation. For a
# readable raw file they come from get_isolation_offsets(), as before, which
# also handles Agilent files that record none. Otherwise they are the first
# offsets recorded in the spectra, the same first occurrence the file scan
# finds.
.msp_isolation_offsets <- function(sp) {
    f <- unique(Spectra::dataStorage(sp))[1]
    if (is(sp@backend, "MsBackendMzR") && !is.na(f) && file.exists(f))
        return(get_isolation_offsets(f))
    h <- .msp_header(sp)
    lo <- h$isolationWindowLowerOffset[!is.na(h$isolationWindowLowerOffset)]
    up <- h$isolationWindowUpperOffset[!is.na(h$isolationWindowUpperOffset)]
    if (!length(lo) || !length(up))
        stop("No isolation window offsets recorded; set offsets explicitly.")
    c(lo[1], up[1])
}

# Peak lists exported from Thermo MSFileReader (msfr) as CSV: one row per
# peak with columns scanid, mz, i and, optionally, background, noise and snr.
# Read as one spectrum per scan, the numeric columns other than scanid as
# peak variables.
.msp_read_msfr <- function(filePth) {
    csv <- utils::read.csv(filePth)
    if (!all(c("scanid", "mz", "i") %in% names(csv)))
        stop("'", filePth, "' is not an msfr peak list: it needs the ",
             "columns scanid, mz and i.", call. = FALSE)
    scans <- unique(csv$scanid)
    rows <- split(seq_len(nrow(csv)), factor(csv$scanid, levels = scans))
    pcols <- c("mz", "intensity",
               setdiff(names(csv), c("scanid", "mz", "i")))
    csv$intensity <- csv$i
    peaks <- lapply(rows, function(r) {
        m <- vapply(pcols, function(cn) as.numeric(csv[[cn]][r]),
                    numeric(length(r)))
        matrix(m, ncol = length(pcols), dimnames = list(NULL, pcols))
    })
    sv <- data.frame(msLevel = rep(1L, length(scans)), scanid = scans,
                     dataOrigin = rep(normalizePath(filePth), length(scans)))
    sp <- .msp_memory(sv, unname(peaks), pcols)
    sp@metadata$columns <- names(csv)[names(csv) != "intensity"]
    sp
}

# The msfr peak list of a Spectra object from .msp_read_msfr(), as read.csv()
# returns it.
.msp_msfr_frame <- function(sp) {
    cols <- sp@metadata$columns
    pk <- .msp_peak_list(sp)
    n <- vapply(pk, nrow, integer(1))
    all <- do.call(rbind, pk)
    out <- lapply(cols, function(cn) {
        if (cn == "scanid") rep(sp$scanid, n)
        else unname(all[, if (cn == "i") "intensity" else cn])
    })
    names(out) <- cols
    as.data.frame(out, optional = TRUE)
}

# mzRback is no longer used; MsBackendMzR always reads with the pwiz reader.
.msp_deprecate_mzRback <- function(mzRback) {
    if (!is.null(mzRback) && !is.na(mzRback) && !identical(mzRback, "pwiz"))
        message("The mzRback argument is deprecated and ignored: raw data ",
                "is read through Spectra, which uses the pwiz reader.")
    invisible(NULL)
}
