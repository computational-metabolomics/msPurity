#' Construct Mass Spectral Trees from measured acquisition metadata
#'
#' Represent measured spectra and their validated acquisition parents as
#' portable node and edge tables using spectrum metadata.
#'
#' @param object A \code{Spectra} object containing measured scans, or a
#'   \code{purityA} object whose source scans have retained acquisition metadata.
#' @param parentReference Interpretation of recorded parent references. One of
#'   \code{"auto"}, \code{"acquisitionNum"}, \code{"scanIndex"},
#'   \code{"spectrumId"} or \code{"sourceSpectrumId"}. Native
#'   \code{spectrumRef} values identify parents by \code{spectrumId}.
#'   Ambiguous automatic interpretations fail with an error.
#'
#' @details
#' Each supplied scan becomes one node. Only valid parent relationships resolved
#' within the supplied scans become edges. Missing, excluded and invalid parents
#' retain their resolution status in \code{parentStatus}.
#' The supported acquisition levels are MS1 through MS4.
#'
#' For \code{purityA} inputs, the source scans and their retained acquisition
#' preparation define the tree. Valid cached preparation is reused, and its
#' stored scan order determines node order. Changes to acquisition metadata
#' invalidate that preparation. Feature links and filtering flags do not remove
#' source scans from this tree. Supply a subset of the measured \code{Spectra}
#' instead when a smaller eligible set is required.
#'
#' For \code{Spectra} inputs, nodes follow input order. Each scan needs a known
#' run, supplied by \code{run_id}, \code{dataOrigin} or \code{fileid}. Existing
#' source keys or \code{lineageScanId} values are retained. Otherwise lineage
#' identifiers use the run and native \code{spectrumId}, or a complete set of
#' positive integer \code{scanIndex} values. With neither native field,
#' retained \code{pid} values preserve msPurity's existing lineage keys.
#' Generic inputs without \code{pid} can use per-run \code{acquisitionNum}.
#' Partial native identity fields require an explicit lineage key.
#' Native spectrum identifiers, source keys and object-local \code{pid} values
#' remain separate columns.
#'
#' Averaged and consensus spectra are outside this acquisition-tree interface.
#' Inputs carrying an averaging level or a \code{combineSpectra} processing
#' record are rejected. Other derived inputs must not be supplied as measured
#' scans. Peak filtering of measured scans is allowed and does not change
#' acquisition relationships.
#'
#' @return A \code{MassSpectralTrees} list containing:
#' \itemize{
#'   \item \code{nodes}: one row per scan, keyed by \code{lineageScanId}, with
#'     original identity fields, acquisition metadata and parent status. A stored
#'     \code{inPurity} score is included when supplied. The
#'     \code{pid} column is \code{NA} when the input did not provide it.
#'   \item \code{edges}: directed acquisition relationships from
#'     \code{parentLineageScanId} to \code{lineageScanId}, with the recorded
#'     parent reference and its interpretation.
#'   \item \code{metadata}: the structure, input type, requested parent reference
#'     and resolved \code{parentReference} convention.
#' }
#' Empty inputs return the same columns with zero rows.
#'
#' @examples
#' metadata <- S4Vectors::DataFrame(
#'     spectrumId = c("scan=10", "scan=20", "scan=30"),
#'     dataOrigin = "example.mzML", msLevel = 1:3,
#'     acquisitionNum = c(10L, 20L, 30L), scanIndex = 1:3,
#'     precScanNum = c(0L, 10L, 20L), rtime = c(1, 2, 3),
#'     precursorMz = c(NA_real_, 300, 150))
#' metadata$mz <- list(300, 150, 75)
#' metadata$intensity <- list(100, 80, 40)
#' spectra <- Spectra::Spectra(metadata)
#' trees <- massSpectralTrees(spectra, parentReference = "acquisitionNum")
#' trees$nodes[, c("lineageScanId", "msLevel", "parentStatus")]
#' trees$edges
#' @export
massSpectralTrees <- function(object,
        parentReference = c("auto", "acquisitionNum", "scanIndex",
            "spectrumId", "sourceSpectrumId")) {
    parentReference <- match.arg(parentReference)
    if (methods::is(object, "purityA")) {
        .mst_measured_input(object@spectra, Spectra::spectraVariables(object@spectra))
        prepared <- .pa_msn_preparation(object, parentReference)
        supplied_pid <- TRUE
        source_type <- "purityA"
    } else if (methods::is(object, "Spectra")) {
        variables <- Spectra::spectraVariables(object)
        .mst_measured_input(object, variables)
        supplied_pid <- "pid" %in% variables
        scans <- .pa_msn_spectra_metadata(object, requirePid = FALSE)
        scans <- .mst_source_identities(object, scans, variables)
        prepared <- .pa_msn_ancestry(scans, parentReference)
        source_type <- "Spectra"
    } else {
        stop("massSpectralTrees requires a Spectra or purityA object.",
            call. = FALSE)
    }
    tree <- .mst_tables(prepared, supplied_pid, source_type)
    tree$metadata$requestedParentReference <- if (source_type == "purityA" &&
            parentReference == "auto" && prepared$parentReference %in%
            c("acquisitionNum", "scanIndex", "sourceSpectrumId")) prepared$parentReference
        else parentReference
    # Cached ancestry can predate an assessment or use a different row order.
    if (source_type == "purityA") {
        source <- object@spectra
        rows <- match(tree$nodes$pid, source$pid)
    } else {
        source <- object
        rows <- seq_along(source)
    }
    if ("inPurity" %in% Spectra::spectraVariables(source))
        tree$nodes$inPurity <- source$inPurity[rows]
    tree
}

.mst_measured_input <- function(object, variables) {
    if (!length(object)) return(invisible(NULL))
    averaged <- "av_level" %in% variables &&
        any(!is.na(object$av_level) & nzchar(as.character(object$av_level)))
    combined <- any(grepl("combineSpectra", Spectra::processingLog(object),
        fixed = TRUE))
    if (averaged || combined)
        stop("Averaged or consensus spectra are not supported; supply ",
            "measured acquisition scans.", call. = FALSE)
    invisible(NULL)
}

.mst_source_identities <- function(object, scans, variables) {
    if (!nrow(scans)) return(scans)
    known_run <- rep(FALSE, nrow(scans))
    for (field in intersect(c("run_id", "dataOrigin", "fileid"), variables)) {
        value <- as.character(object[[field]])
        known_run <- known_run | (!is.na(value) & nzchar(value))
    }
    if (any(!known_run))
        stop("Measured spectra require a known run_id, dataOrigin or fileid ",
            "for stable acquisition identities.", call. = FALSE)
    if (any(c("spectrum_id_", "lineageScanId") %in% variables)) return(scans)

    native <- if ("spectrumId" %in% variables)
        as.character(object$spectrumId) else rep(NA_character_, nrow(scans))
    known_native <- !is.na(native) & nzchar(native)
    if (any(known_native)) {
        if (!all(known_native))
            stop("Incomplete spectrumId values cannot define stable ",
                "acquisition identities; supply retained lineageScanId ",
                "values for every scan.", call. = FALSE)
        return(scans)
    }
    for (field in c("scan_index", "acquisition_num")) {
        if (field == "acquisition_num" && "pid" %in% variables) return(scans)
        value <- scans[[field]]
        if (all(is.na(value))) next
        if (!is.numeric(value) || anyNA(value) ||
                any(!is.finite(value) | value < 1 | value != floor(value)))
            stop("Complete positive integer scanIndex or acquisitionNum ",
                "values are required for stable acquisition identities.",
                call. = FALSE)
        if (field == "acquisition_num") {
            scans$lineageScanId <- paste(
                utils::URLencode(scans$run_id, reserved = TRUE),
                paste0("acquisitionNum:", as.character(value)), sep = "|")
        }
        return(scans)
    }
    stop("Measured spectra need stable source keys, native spectrumId, ",
        "scanIndex or acquisitionNum values; pid or input position alone ",
        "cannot define stable acquisition identities.", call. = FALSE)
}

.mst_tables <- function(prepared, supplied_pid, source_type) {
    scans <- prepared$scans
    nodes <- data.frame(
        lineageScanId = as.character(scans$lineageScanId),
        parentLineageScanId = as.character(scans$parentLineageScanId),
        parentStatus = as.character(scans$parentStatus),
        parentReferenceType = as.character(scans$parentReferenceType),
        rawParentReference = as.character(scans$rawParentReference),
        pid = if (supplied_pid) as.numeric(scans$pid) else
            rep(NA_real_, nrow(scans)),
        spectrumId = as.character(scans$spectrumId),
        sourceSpectrumId = as.character(scans$sourceSpectrumId),
        run_id = as.character(scans$run_id),
        dataOrigin = as.character(scans$data_origin),
        acquisitionNum = as.numeric(scans$acquisition_num),
        scanIndex = as.numeric(scans$scan_index),
        msLevel = as.integer(scans$ms_level),
        rtime = as.numeric(scans$rtime),
        precursorMz = as.numeric(scans$precursor_mz),
        stringsAsFactors = FALSE)
    linked <- !is.na(nodes$parentLineageScanId)
    edges <- nodes[linked, c("parentLineageScanId", "lineageScanId",
        "parentStatus", "parentReferenceType", "rawParentReference"), drop = FALSE]
    rownames(edges) <- NULL
    edges$relationship <- rep("acquisition", nrow(edges))
    structure(list(nodes = nodes, edges = edges,
        metadata = list(structure = "MST", sourceType = source_type,
            parentReference = prepared$parentReference)),
        class = "MassSpectralTrees")
}
