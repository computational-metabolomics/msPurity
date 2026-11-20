# One stage's conditions can hold sets of contributor values for an average.
.sm_condition_fields <- setdiff(names(.PA_FILTER_AXES), "msLevel")
.sm_summary_fields <- c(msPurityPrecursorPath = "x_mspurity_precursor_path",
    msPurityAncestorConditions = "x_mspurity_ancestor_conditions",
    msPurityAcquisitionValues = "x_mspurity_acquisition_values",
    msPurityAcquisitionStatus = "x_mspurity_acquisition_status")

.sm_summary <- function(row, field) {
    live <- row[[field]]
    if (is.list(live)) return(live)
    json <- row[[.sm_summary_fields[[field]]]]
    if (length(json) && !is.na(json) && nzchar(json))
        jsonlite::fromJSON(json)
}

.sm_metadata_row <- function(x, i) {
    fields <- c(unlist(.PA_FILTER_AXES), "precursorMz", "activationTypeAccession",
        "massAnalyzerAccession", "collisionEnergyUnitAccession", "filterString",
        names(.sm_summary_fields), .sm_summary_fields, "x_mspurity_activation_type")
    fields <- intersect(fields, Spectra::spectraVariables(x))
    d <- Spectra::spectraData(x[i], columns = fields)
    lapply(as.list(d), `[[`, 1L)
}

.sm_normalize <- function(x, field) {
    if (field %in% c("collisionEnergy", "polarity")) return(as.numeric(x))
    x <- tolower(trimws(as.character(x)))
    x[x %in% c("", "unknown", "mixed")] <- NA_character_
    if (field == "activationType") return(.pa_activation_type(x))
    if (field == "collisionEnergyUnit") {
        x[x %in% c("electronvolt", "electronvolts", "electron volt", "uo:0000266")] <- "ev"
        x[x %in% c("normalized collision energy", "normalised collision energy",
            "ms:1000138")] <- "nce"
        x[x %in% c("percent", "uo:0000187")] <- "%"
    }
    if (field == "massAnalyzer")
        x[x %in% "ms:1000484"] <- "orbitrap"
    x
}

.sm_stage <- function(row, values = NULL, status = NULL) {
    # The package's acquisition reader resolves known aliases, CV activation
    # methods and scan filter strings. Preserve list-valued consensus metadata
    # separately from the scalar metadata passed to that reader.
    fields <- intersect(names(row), c(unlist(.PA_FILTER_AXES),
        "activationTypeAccession", "filterString"))
    scalar <- as.data.frame(row[fields], stringsAsFactors = FALSE)
    canonical <- .pa_filter_metadata(scalar)
    stage <- list(precursorMz = row$precursorMz)
    cv_fields <- c("activationType", "massAnalyzer", "collisionEnergyUnit")
    for (field in .sm_condition_fields) {
        # An explicitly unknown contributor field remains unknown after
        # reopening, even when the backend supplies a run-wide scalar default.
        summarized <- field %in% names(values) || identical(status[[field]], "unknown")
        v <- if (summarized) values[[field]] else canonical[[field]]
        if (!summarized && field == "activationType" && all(is.na(v)))
            v <- row$x_mspurity_activation_type
        if (field %in% cv_fields) {
            key <- paste0(field, "Accession")
            id <- values[[key]]
            if (!key %in% names(values)) id <- row[[key]]
            if (all(is.na(id))) id <- v[grepl("^[a-zA-Z]+:[0-9]+$", v)]
            id <- toupper(trimws(as.character(id)))
            id[!grepl("^[A-Z]+:[0-9]+$", id)] <- NA_character_
            stage[[key]] <- unique(id)
            if (all(is.na(v))) v <- id
        }
        v <- .sm_normalize(v, field)
        if (identical(status[[field]], "mixed")) v <- c(v, NA)
        stage[[field]] <- unique(v)
    }
    stage
}

.sm_parent_columns <- c("msLevel", "run_id", "dataOrigin", "fileid",
    "spectrumRef", "spectrumId", "nativeId", "parentLineageScanId", "lineageScanId",
    "parentPid", "pid",
    "precScanNum", "precursorScanNum", "acquisitionNum", "scanIndex",
    "sourceSpectrumId", "spectrum_id_", "rtime")

.sm_parent_rows <- function(x, index, level) {
    cols <- intersect(.sm_parent_columns, Spectra::spectraVariables(x))
    d <- as.data.frame(Spectra::spectraData(x, columns = cols))
    run <- rep(NA_character_, nrow(d))
    for (key in intersect(c("run_id", "dataOrigin", "fileid"), cols)) {
        v <- as.character(d[[key]])
        use <- is.na(run) & !is.na(v) & nzchar(v)
        run[use] <- paste(key, v[use], sep = "=")
    }
    # A collection with no run metadata can still resolve a unique local
    # reference. Never fill a missing run with some other entry's run.
    if (all(is.na(run))) run[] <- "input"
    rows <- rep(NA_integer_, level)
    rows[level] <- index
    for (stage in seq.int(level, 3L)) {
        child <- rows[stage]
        candidates <- which(run == run[child] & d$msLevel == stage - 1L)
        if (!length(candidates)) break
        for (key in intersect(c("rtime", "scanIndex", "acquisitionNum"), cols)) {
            v <- d[[key]]
            if (!is.finite(v[child])) next
            before <- if (key == "rtime") v <= v[child] else v < v[child]
            candidates <- candidates[is.na(before[candidates]) | before[candidates]]
        }
        keys <- list(parentLineageScanId = "lineageScanId", parentPid = "pid",
            spectrumRef = c("spectrumId", "nativeId"),
            precScanNum = c("acquisitionNum", "scanIndex", "sourceSpectrumId"),
            precursorScanNum = c("acquisitionNum", "scanIndex", "sourceSpectrumId"))
        found <- integer()
        for (key in intersect(names(keys), cols)) {
            ref <- d[[key]][child]
            if (is.na(ref) || !nzchar(ref) || ref %in% c("0", "-1")) next
            ids <- intersect(keys[[key]], cols)
            # A dataset-assigned spectrum ID is independent of the raw scan number.
            if ("spectrum_id_" %in% cols) ids <- setdiff(ids, "sourceSpectrumId")
            found <- unique(unlist(lapply(ids, function(id)
                candidates[which(d[[id]][candidates] == ref)]), use.names = FALSE))
            break
        }
        if (length(found) != 1L) break
        rows[stage - 1L] <- found
    }
    rows
}

.sm_stages <- function(x, index, level) {
    stages <- vector("list", level)
    row <- .sm_metadata_row(x, index)
    values <- .sm_summary(row, "msPurityAcquisitionValues")
    stages[[level]] <- .sm_stage(row, values, .sm_summary(row, "msPurityAcquisitionStatus"))
    if (level == 2L) return(stages)
    path <- .sm_summary(row, "msPurityPrecursorPath")
    ancestors <- .sm_summary(row, "msPurityAncestorConditions")
    if (is.null(path) && is.null(ancestors) && is.null(values)) {
        parents <- .sm_parent_rows(x, index, level)
        for (i in which(!is.na(parents) & seq_along(parents) < level))
            stages[[i]] <- .sm_stage(.sm_metadata_row(x, parents[i]))
        return(stages)
    }
    for (i in seq.int(2L, level - 1L)) {
        j <- match(i, ancestors$msLevel)
        if (!is.na(j)) {
            # JSON encodes the stage table by columns and nulls for unknowns.
            ancestor <- lapply(ancestors, function(v) if (length(v)) v[j] else NA)
            stages[[i]] <- .sm_stage(ancestor)
        }
        j <- match(i, path$msLevel)
        if (!is.na(j)) stages[[i]]$precursorMz <- path$precursorMz[j]
    }
    stages
}

.sm_field_check <- function(a, b, tolerance = NULL) {
    known_a <- a[!is.na(a)]
    known_b <- b[!is.na(b)]
    if (!length(known_a) || !length(known_b)) return("missing")
    shared <- if (is.null(tolerance)) any(known_a %in% known_b) else
        any(abs(outer(known_a, known_b, "-")) <= tolerance)
    if (!shared) {
        # An unrecognized accession and a free-text label may describe the
        # same term. Only comparable representations establish a conflict.
        if (is.character(known_a) && is.character(known_b)) {
            cv <- "^[a-zA-Z]+:[0-9]+$"
            if (any(outer(grepl(cv, known_a), grepl(cv, known_b), "!=")))
                return("partial")
        }
        return("mismatch")
    }
    if (length(a) == 1L && length(b) == 1L) "match" else "partial"
}

.sm_condition_checks <- function(a, b, energy_tolerance) {
    out <- stats::setNames(character(length(.sm_condition_fields)), .sm_condition_fields)
    for (field in setdiff(names(out), "collisionEnergy")) {
        outcome <- .sm_field_check(a[[field]], b[[field]])
        key <- paste0(field, "Accession")
        # Uniform accessions resolve label synonyms. Preserve uncertainty when
        # a contributor set also contains unknowns or several conditions.
        if (length(a[[key]]) == 1L && length(b[[key]]) == 1L &&
                !anyNA(c(a[[key]], b[[key]])) &&
                length(a[[field]]) == 1L && length(b[[field]]) == 1L &&
                !anyNA(c(a[[field]], b[[field]])))
            outcome <- .sm_field_check(a[[key]], b[[key]])
        out[field] <- outcome
    }
    # With several scales the contributor summary does not identify which
    # energy belongs to which unit. Only compare values on one known scale.
    out["collisionEnergy"] <- if (out["collisionEnergyUnit"] == "match")
        .sm_field_check(a$collisionEnergy, b$collisionEnergy, energy_tolerance)
    else "unit_unverified"
    out
}
