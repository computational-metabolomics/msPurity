# Pids identify spectra within the supplied
# object; spectrum_id is a unique source key, not an acquisition number.

# Retain scalar source metadata and keep dataset keys separate from native
# acquisition identifiers. All values are read through public Spectra accessors.
.pa_msn_spectra_metadata <- function(object) {
    if (!methods::is(object, "Spectra"))
        stop("MSn acquisition preparation requires a Spectra object.")
    columns <- setdiff(Spectra::spectraVariables(object),
        Spectra::peaksVariables(object))
    metadata <- as.data.frame(Spectra::spectraData(object, columns = columns),
        optional = TRUE)
    n <- nrow(metadata)
    column <- function(name, default = NA_real_) {
        if (name %in% names(metadata)) metadata[[name]] else rep(default, n)
    }
    if (!"pid" %in% names(metadata)) {
        if (n) stop("MSn source spectra require retained pid identifiers.")
        metadata$pid <- integer()
    }
    run <- as.character(column("run_id", NA_character_))
    origin <- as.character(column("dataOrigin", NA_character_))
    missing <- is.na(run) | !nzchar(run)
    run[missing] <- origin[missing]
    missing <- is.na(run) | !nzchar(run)
    file <- as.character(column("fileid", NA_character_))
    run[missing] <- ifelse(is.na(file[missing]), "input",
        paste0("fileid:", file[missing]))
    source_key <- NULL
    if ("spectrum_id_" %in% names(metadata)) {
        source_key <- as.character(metadata$spectrum_id_)
        id <- source_key
    } else if ("spectrum_id" %in% names(metadata)) {
        id <- as.character(metadata$spectrum_id)
        source_key <- as.character(column("source_spectrum_id", NA_character_))
        source_key[is.na(source_key)] <- id[is.na(source_key)]
    } else {
        native <- as.character(column("spectrumId", NA_character_))
        if (all(!is.na(native) & nzchar(native))) {
            key <- native
        } else {
            index <- column("scanIndex")
            key <- if (all(is.finite(index))) as.character(index) else
                paste0("pid:", metadata$pid)
        }
        id <- paste(utils::URLencode(run, reserved = TRUE),
            utils::URLencode(key, reserved = TRUE), sep = "|")
        if (is.null(source_key)) source_key <- rep(NA_character_, n)
    }
    if (anyNA(id) || any(!nzchar(id)) || anyDuplicated(id))
        stop("MSn source spectra require unique known source keys.")
    normalized <- .pa_filter_metadata(metadata)
    for (name in names(normalized)) metadata[[name]] <- normalized[[name]]
    metadata$spectrum_id <- id
    metadata$source_spectrum_id <- source_key
    metadata$run_id <- run
    metadata$data_origin <- origin
    metadata$source_index <- seq_len(n)
    mapping <- c(ms_level = "msLevel", acquisition_num = "acquisitionNum",
        scan_index = "scanIndex", rtime = "rtime", precursor_mz = "precursorMz",
        precursor_charge = "precursorCharge",
        precursor_intensity = "precursorIntensity")
    for (name in names(mapping)) metadata[[name]] <- column(mapping[[name]])
    metadata$parent_reference_raw <- as.character(column("precScanNum"))
    metadata
}

.pa_msn_parent_candidates <- function(scans, convention) {
    reference <- if (convention == "spectrumId") scans$spectrumRef else
        scans$parent_reference_raw
    column <- switch(convention, spectrumId = "spectrumId",
        sourceSpectrumId = "source_spectrum_id",
        acquisitionNum = "acquisition_num", scanIndex = "scan_index")
    value <- as.character(scans[[column]])
    key <- paste(scans$run_id, value, sep = "\r")
    known <- !is.na(value) & nzchar(value)
    duplicate <- anyDuplicated(key[known]) > 0L
    rows <- which(known)
    parent <- if (duplicate) rep(NA_integer_, nrow(scans)) else
        rows[match(paste(scans$run_id, reference, sep = "\r"), key[rows])]
    parent[is.na(reference) | scans$ms_level == 1L] <- NA_integer_
    valid <- !is.na(parent)
    rows <- which(valid)
    p <- scans[parent[rows], , drop = FALSE]
    child <- scans[rows, , drop = FALSE]
    index_known <- is.finite(p$scan_index) & is.finite(child$scan_index)
    before <- ifelse(index_known, p$scan_index < child$scan_index,
        !is.finite(p$acquisition_num) | !is.finite(child$acquisition_num) |
            p$acquisition_num < child$acquisition_num)
    valid[rows] <- p$ms_level < child$ms_level & before &
        (!is.finite(p$rtime) | !is.finite(child$rtime) |
                p$rtime <= child$rtime) &
        (is.na(p$polarity) | is.na(child$polarity) |
                p$polarity == child$polarity)
    list(parent = parent, valid = valid, duplicate = duplicate)
}

# Required columns: pid, spectrum_id, run_id, ms_level. Optional numeric
# metadata are acquisition_num, scan_index, rtime, precursor_mz and polarity;
# parent_reference_raw holds the declared source reference. No peaks are read.
.pa_msn_ancestry <- function(scans, parentReference = "auto") {
    parentReference <- match.arg(parentReference,
        c("auto", "acquisitionNum", "scanIndex", "spectrumId",
            "sourceSpectrumId"))
    required <- c("pid", "spectrum_id", "run_id", "ms_level")
    if (!is.data.frame(scans) || !all(required %in% names(scans)))
        stop("MSn metadata require pid, spectrum_id, run_id and ms_level.")
    if (!is.numeric(scans$pid) || any(!is.finite(scans$pid)) ||
            anyDuplicated(scans$pid))
        stop("MSn spectrum pids must be unique finite numbers.")
    for (column in c("spectrum_id", "run_id")) {
        scans[[column]] <- as.character(scans[[column]])
        if (anyNA(scans[[column]]) || any(!nzchar(scans[[column]])))
            stop("MSn source identities and runs must be known.")
    }
    if (anyDuplicated(scans$spectrum_id))
        stop("MSn source spectrum identities must be unique.")
    if (!is.numeric(scans$ms_level) || anyNA(scans$ms_level) ||
            any(!scans$ms_level %in% 1:4))
        stop("This MSn increment supports known MS1-MS4 levels.")
    for (column in c("acquisition_num", "scan_index", "rtime",
            "precursor_mz", "polarity")) {
        if (!column %in% names(scans) || all(is.na(scans[[column]])))
            scans[[column]] <- rep(NA_real_, nrow(scans))
        if (!is.numeric(scans[[column]]) ||
                any(!is.na(scans[[column]]) & !is.finite(scans[[column]])))
            stop("MSn metadata ", column, " must be numeric and finite or NA.")
    }
    if (!"parent_reference_raw" %in% names(scans))
        scans$parent_reference_raw <- rep(NA_character_, nrow(scans))
    scans$parent_reference_raw <- as.character(scans$parent_reference_raw)
    scans$parent_reference_raw[scans$parent_reference_raw %in%
            c("", "0", "-1")] <- NA_character_
    for (column in c("spectrumId", "spectrumRef", "source_spectrum_id")) {
        if (!column %in% names(scans))
            scans[[column]] <- rep(NA_character_, nrow(scans))
        scans[[column]] <- as.character(scans[[column]])
        scans[[column]][!is.na(scans[[column]]) &
                !nzchar(scans[[column]])] <- NA_character_
    }
    conventions <- c("sourceSpectrumId", "acquisitionNum", "scanIndex",
        "spectrumId")
    candidates <- stats::setNames(lapply(conventions, function(convention) {
        .pa_msn_parent_candidates(scans, convention)
    }), conventions)
    native <- candidates$spectrumId
    has_native <- !is.na(scans$spectrumRef) & scans$ms_level > 1L
    convention <- parentReference
    if (convention == "auto") {
        for (i in seq_len(nrow(scans))) {
            resolved <- vapply(candidates, function(candidate) {
                if (candidate$valid[i]) candidate$parent[i] else NA_integer_
            }, integer(1L))
            if (length(unique(resolved[!is.na(resolved)])) > 1L)
                stop("Ambiguous parent references; specify parentReference.")
        }
        numeric_candidates <- candidates[conventions != "spectrumId"]
        counts <- vapply(numeric_candidates,
            function(candidate) sum(candidate$valid & !has_native), integer(1L))
        convention <- names(numeric_candidates)[which.max(counts)]
        covered <- Reduce(`|`, lapply(numeric_candidates, `[[`, "valid"))
        if (any(covered & !has_native & !candidates[[convention]]$valid))
            stop("Inconsistent parent conventions; specify parentReference.")
    }
    chosen <- candidates[[convention]]
    if (chosen$duplicate && any(!is.na(scans$parent_reference_raw) &
                (parentReference != "auto" | !has_native)))
        stop("Duplicate ", convention, " keys within a source run.")
    if (native$duplicate && any(has_native))
        stop("Duplicate spectrumId keys within a source run.")
    if (parentReference != "auto" && convention != "spectrumId" &&
            any(native$valid & chosen$valid &
                    native$parent != chosen$parent, na.rm = TRUE))
        stop("Native spectrumRef conflicts with ", convention, " references.")
    reference_type <- rep(convention, nrow(scans))
    if (parentReference == "auto" && any(has_native)) {
        chosen$parent[has_native] <- native$parent[has_native]
        chosen$valid[has_native] <- native$valid[has_native]
        reference_type[has_native] <- "spectrumId"
    }
    use_native <- reference_type == "spectrumId"
    scans$parent_reference_raw[use_native] <- scans$spectrumRef[use_native]
    scans$parentReferenceType <- reference_type
    referenced <- scans$ms_level > 1L & !is.na(scans$parent_reference_raw)
    scans$parent_status <- ifelse(scans$ms_level == 1L, "root",
        "missing_reference")
    scans$parent_status[referenced] <- ifelse(
        is.na(chosen$parent[referenced]), "unresolved_reference",
        "invalid_reference")
    scans$parent_status[chosen$valid] <- "valid"
    scans$parent_id <- rep(NA_character_, nrow(scans))
    scans$parent_id[chosen$valid] <- scans$spectrum_id[
        chosen$parent[chosen$valid]]
    parent_rows <- chosen$parent
    parent_rows[!chosen$valid] <- NA_integer_
    scans$parent_pid <- scans$pid[parent_rows]
    root_pid <- scans$pid[rep(NA_integer_, nrow(scans))]
    path_status <- ifelse(scans$ms_level == 2L, "root",
        ifelse(scans$ms_level == 1L, "not_fragmentation", scans$parent_status))
    paths <- rep(list(NULL), nrow(scans))
    for (i in which(scans$ms_level >= 2L)) {
        path <- i
        at <- i
        while (scans$ms_level[at] > 2L) {
            if (scans$parent_status[at] != "valid") {
                path_status[i] <- scans$parent_status[at]
                break
            }
            parent <- parent_rows[at]
            if (scans$ms_level[parent] != scans$ms_level[at] - 1L) {
                path_status[i] <- "missing_stage"
                break
            }
            path <- c(parent, path)
            at <- parent
        }
        if (scans$ms_level[at] != 2L)
            next
        if (scans$ms_level[i] > 2L &&
                any(!is.finite(scans$precursor_mz[path]) |
                        scans$precursor_mz[path] <= 0)) {
            path_status[i] <- "missing_precursor"
            next
        }
        root_pid[i] <- scans$pid[at]
        path_status[i] <- if (i == at) "root" else "valid"
        paths[i] <- list(scans$pid[path])
    }
    scans$root_pid <- root_pid
    scans$path_status <- path_status
    scans$path <- paths
    if (parentReference == "auto") {
        resolved <- unique(reference_type[chosen$valid])
        convention <- if (!length(resolved)) "unresolved" else
            if (length(resolved) == 1L) resolved else "mixed"
    }
    list(scans = scans, parentReference = convention)
}

# Cached ancestry is valid only for the same acquisition metadata. Processing
# flags and history may change independently and are not acquisition fields.
.pa_msn_validate_preparation <- function(pa) {
    preparation <- pa@params$msnLineage
    if (is.null(preparation)) return(invisible(TRUE))
    stored <- preparation$source_metadata
    if (!is.data.frame(stored) || !"pid" %in% names(stored))
        stop("Stored MSn ancestry lacks its source metadata. ",
            "Run frag4feature(includeMSn = TRUE) again.")
    current <- .pa_msn_spectra_metadata(pa@spectra)
    if (anyDuplicated(current$pid) || !setequal(current$pid, stored$pid))
        stop("Stored MSn ancestry no longer matches the source spectra. ",
            "Run frag4feature(includeMSn = TRUE) again.")
    fields <- unique(c(.PA_ACQUISITION_FIELDS, .MSP_HEADER_COLS,
            unname(.MSP_HEADER_RENAME), unlist(.PA_FILTER_AXES),
            "spectrum_id", "spectrum_id_", "source_spectrum_id", "run_id",
            "data_origin", "ms_level", "acquisition_num", "scan_index",
            "precursor_mz", "parent_reference_raw", "spectrumId", "spectrumRef",
            "isolationWindowLowerMz",
            "isolationWindowUpperMz"))
    fields <- intersect(fields, union(names(current), names(stored)))
    same_fields <- all(fields %in% names(current)) &&
        all(fields %in% names(stored))
    if (!same_fields || !isTRUE(all.equal(current[, fields, drop = FALSE],
                stored[match(current$pid, stored$pid), fields, drop = FALSE],
                check.attributes = FALSE, tolerance = 0)))
        stop("Stored MSn ancestry no longer matches acquisition metadata. ",
            "Run frag4feature(includeMSn = TRUE) again.")
    invisible(TRUE)
}

# Acquisition fields have one source: the current spectrum, not its root's
# feature link. Reuse the header and acquisition vocabularies already retained
# by msPurity and the same alias normalization as parameter selection.
.pa_msn_stage_metadata <- function(scans) {
    mapping <- c(acquisitionNum = "acquisition_num",
        scanIndex = "scan_index", seqNum = "scan_index", rtime = "rtime",
        retentionTime = "rtime", precursorMZ = "precursor_mz",
        precursorMz = "precursor_mz", precursor_mz = "precursor_mz",
        precursorIntensity = "precursor_intensity",
        precursorCharge = "precursor_charge", massAnalyzer = "mass_analyzer",
        activationType = "activation_type",
        collisionEnergy = "collision_energy",
        collisionEnergyUnit = "collision_energy_unit",
        dataOrigin = "data_origin")
    for (field in names(mapping)) {
        source <- mapping[[field]]
        if (source %in% names(scans)) scans[[field]] <- scans[[source]]
    }
    scans$msLevel <- scans$ms_level
    normalized <- .pa_filter_metadata(scans)
    for (field in names(normalized)) scans[[field]] <- normalized[[field]]
    scans
}

# links are pre-purity MS2 association occurrences, with pid and grpid.
# Descendants inherit feature membership, never the root's measured purity.
.pa_msn_feature_links <- function(links, ancestry, plim = NA_real_) {
    scans <- ancestry$scans
    if (!is.data.frame(links) || !all(c("pid", "grpid") %in% names(links)))
        stop("Root feature links require pid and grpid.")
    if (identical(plim, NA)) plim <- NA_real_
    if (!is.numeric(plim) || length(plim) != 1L ||
            (!is.na(plim) && !is.finite(plim)))
        stop("plim must be one finite number or NA.")
    root <- match(links$pid, scans$pid)
    if (anyNA(root) || any(scans$ms_level[root] != 2L))
        stop("Root feature links must refer to supplied MS2 spectra.")
    if (!"inPurity" %in% names(links))
        links$inPurity <- rep(NA_real_, nrow(links))
    if (!"featureMz" %in% names(links)) {
        column <- intersect(c("mzmed", "mz"), names(links))
        links$featureMz <- if (length(column)) links[[column[1L]]] else
            rep(NA_real_, nrow(links))
    }
    links$msLevel <- rep(2L, nrow(links))
    links$msPurityRootPid <- links$pid
    for (column in setdiff(c("seqNum", "retentionTime", "precurMtchID",
                    "precurMtchScan", "precurMtchRT", "precurMtchMZ",
                    "precurMtchPPM"), names(links)))
        links[[column]] <- rep(NA_real_, nrow(links))
    purity_columns <- intersect(c("precursorNearest", "aMz", "aPurity",
            "apkNm", "iMz", "iPurity", "ipkNm", "inPkNm", "inpkNm",
            "inPurity", "purity_pass_flag"), names(links))
    stage <- .pa_msn_stage_metadata(scans)
    stage_columns <- unique(c(.PA_ACQUISITION_FIELDS, .MSP_HEADER_COLS,
            unname(.MSP_HEADER_RENAME),
            unlist(.PA_FILTER_AXES, use.names = FALSE),
            "spectrum_id", "spectrum_id_", "source_spectrum_id", "source_index",
            "spectrum_index", "run_id", "parent_reference_raw",
            "acquisition_num", "scan_index", "precursor_mz",
            "isolationWindowLowerMz", "isolationWindowUpperMz"))
    feature_columns <- c("grpid", "sample", "cid", "featureMz", "mz", "mzmed",
        "mzmin", "mzmax", "rt", "rtmed", "rtmin", "rtmax", "into", "intb",
        "maxo", "sn", "npeaks", "peakidx", "index", "featureLinkIndex",
        "contributorIndex", "msPurityRootPid")
    source_columns <- setdiff(names(stage), c("parent_id", "parent_pid",
            "parent_status", "root_pid", "path_status", "path",
            feature_columns))
    stage_columns <- union(stage_columns, source_columns)
    stage_columns <- setdiff(stage_columns, purity_columns)
    stage_columns <- intersect(stage_columns, union(names(links), names(stage)))
    for (column in setdiff(stage_columns, names(links)))
        links[[column]] <- stage[[column]][root]
    linked_roots <- unique(root)
    rows_by_root <- split(seq_len(nrow(links)), match(root, linked_roots))
    root_rows <- match(scans$root_pid, scans$pid)
    root_groups <- match(root_rows, linked_roots)
    parent_rows <- match(scans$parent_id, scans$spectrum_id)
    descendants <- which(scans$ms_level > 2L & scans$path_status == "valid" &
            !is.na(root_groups))
    occurrences <- rows_by_root[root_groups[descendants]]
    child_rows <- rep(descendants, lengths(occurrences))
    child <- links[unlist(occurrences, use.names = FALSE), , drop = FALSE]
    if (length(child_rows)) {
        child$pid <- scans$pid[child_rows]
        child$msLevel <- scans$ms_level[child_rows]
        for (column in purity_columns) child[[column]][] <- NA
        for (column in stage_columns) {
            if (column %in% names(stage))
                child[[column]] <- stage[[column]][child_rows]
            else child[[column]][] <- NA
        }
        child$seqNum <- scans$scan_index[child_rows]
        child$retentionTime <- scans$rtime[child_rows]
        child$precurMtchID <- child$seqNum
        reference <- scans$parent_reference_raw[child_rows]
        child$precurMtchScan <- if (all(grepl("^[0-9]+$", reference)))
            as.numeric(reference) else reference
        child$precurMtchRT <- scans$rtime[parent_rows[child_rows]]
        child$precurMtchMZ <- scans$precursor_mz[child_rows]
        child$precurMtchPPM <- rep(NA_real_, length(child_rows))
        if ("precScanNum" %in% names(stage) &&
                "precursorScanNum" %in% names(child))
            child$precursorScanNum <- stage$precScanNum[child_rows]
        if ("parent_reference_raw" %in% names(child))
            child$parent_reference_raw <- reference
        if ("precursorRT" %in% names(child))
            child$precursorRT <- child$precurMtchRT
    }
    result <- rbind(links, child)
    if (!is.na(plim) && plim > 0)
        result <- result[result$msLevel != 2L |
                (!is.na(result$inPurity) & result$inPurity > plim), ,
            drop = FALSE]
    rownames(result) <- NULL
    result
}
