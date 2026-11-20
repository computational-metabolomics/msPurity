# A Spectra object may be a reordered subset, or have queued peak processing.
# Its row numbers identify entries in the supplied object only.
.mzs_spectra_library_manifest <- function(version, checksum) {
    if (!is.null(version) && (!is.character(version) || length(version) != 1L ||
                             is.na(version) || !nzchar(trimws(version))))
        .mzs_abort("format", "'libraryVersion' must be a non-empty release name.")
    if (!is.null(checksum)) {
        if (!is.character(checksum) || length(checksum) != 1L || is.na(checksum) ||
            !grepl("^(sha256:)?[[:xdigit:]]{64}$", checksum))
            .mzs_abort("format", "'libraryChecksum' must be a SHA-256 digest.")
        checksum <- paste0("sha256:", tolower(sub("^sha256:", "", checksum)))
    }
    list(role = "library", resolution = "external",
         library = list(version = version %||% "Spectra", digest = checksum))
}

# Keep only precursor masses and row indices for the whole library. Binary
# searches find a conservative window; the scorer checks exact ppm overlap.
.mzs_direct_library <- function(x, a, key) {
    n <- length(x)
    if (isTRUE(a$usePrecursors)) {
        mz <- Spectra::precursorMz(x)
        index <- order(mz, na.last = NA)
        sorted <- mz[index]
        rm(mz)
    }
    candidates <- function(precursor) {
        if (!isTRUE(a$usePrecursors)) return(seq_len(n))
        if (!is.finite(precursor)) return(integer())
        window <- .mzs_precursor_bounds(precursor, a)
        first <- findInterval(window[1L], sorted, left.open = TRUE, checkSorted = FALSE) + 1L
        last <- findInterval(window[2L], sorted, checkSorted = FALSE)
        if (first > last) return(integer())
        sort.int(index[seq.int(first, last)])
    }
    fields <- .MZS_SPECTRA_LIBRARY_MAP[paste0("x_mspurity_", c(
        "name", "compound_name", "precursor_type", "inchikey", "formula",
        "smiles", "inchi", "exact_mass", "instrument", "instrument_type", "splash"))]
    vars <- intersect(Spectra::spectraVariables(x), c(
        "polarity", "precursorMz", "rtime", "dataOrigin", "accession", "spectrumId",
        unname(fields), names(fields),
        "precursorPurity", "x_mspurity_precursor_purity", "x_mspurity_source_name",
        "x_mspurity_retention_time"))
    metadata <- function(ids) {
        d <- Spectra::spectraData(x[ids], columns = vars)
        # Coalesce per entry. A backend can advertise a variable that is
        # missing on some or all spectra.
        value <- function(..., numeric = FALSE) {
            result <- rep(if (numeric) NA_real_ else NA_character_, length(ids))
            for (v in intersect(c(...), vars)) {
                val <- d[[v]]
                if (is.list(val) || is(val, "List"))
                    val <- vapply(as.list(val), function(e) {
                        if (!length(e) || all(is.na(e))) NA_character_
                        else paste(as.character(e), collapse = "; ")
                    }, "")
                missing <- is.na(result)
                val <- val[missing]
                if (numeric) {
                    if (is.factor(val)) val <- as.character(val)
                    result[missing] <- suppressWarnings(as.numeric(val))
                } else result[missing] <- .mzs_na_empty(val)
            }
            result
        }
        field <- function(name, numeric = FALSE) {
            f <- paste0("x_mspurity_", name)
            value(f, fields[[f]], numeric = numeric)
        }
        accession <- value("accession", "spectrumId", "x_mspurity_splash", "splash")
        pol <- as.integer(d$polarity)
        meta <- data.frame(
            reference_source = key, reference_run_id = NA_character_,
            reference_spectrum_id_ = as.numeric(ids), reference_native_id = accession,
            precursor_mz = as.numeric(d$precursorMz),
            retention_time = value("x_mspurity_retention_time", "rtime", numeric = TRUE),
            inPurity = value("x_mspurity_precursor_purity", "precursorPurity", numeric = TRUE),
            accession = accession, name = field("name"), precursor_type = field("precursor_type"),
            inchikey = field("inchikey"), compound_name = field("compound_name"),
            formula = field("formula"), smiles = field("smiles"), inchi = field("inchi"),
            exact_mass = field("exact_mass", numeric = TRUE),
            source_name = value("x_mspurity_source_name", "dataOrigin"),
            polarity = unname(c("0" = "negative", "1" = "positive")[as.character(pol)]),
            instrument_type = field("instrument_type"), instrument = field("instrument"))
        .mzs_filter_library(meta, a, anonymous = TRUE)
    }
    peaks <- function(ids)
        stats::setNames(.mzs_matching_peaks(x[ids]), as.character(ids))
    list(candidates = candidates, metadata = metadata, peaks = peaks)
}

# Hash only anonymous candidate entries, using a platform-independent binary
# representation. Include precursor mass, polarity and ordered fragment peaks.
.mzs_anonymous_accession <- function(mz, polarity, peaks) {
    pol <- match(polarity, c("negative", "positive"), nomatch = 0L)
    bytes <- writeBin(as.double(c(mz, pol, peaks$mz, peaks$i)), raw(),
                      size = 8L, endian = "big")
    paste0("sha256:", tools::sha256sum(bytes = bytes))
}
