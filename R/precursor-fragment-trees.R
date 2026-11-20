#' Construct a precursor-fragment tree from measured spectra
#'
#' Link each fragmentation scan's observed products to a measured peak in
#' its acquisition parent. A `purityA` input reuses recorded parent peak
#' assignments from MS3/MS4 purity assessments.
#'
#' @inheritParams massSpectralTrees
#' @param object Measured `Spectra`, a `purityA` result (its retained source
#'   spectra), or a scan tree from [massSpectralTrees()].
#' @param ppm,tolerance For scans without a recorded assignment, the
#'   nonnegative relative tolerance in ppm and absolute
#'   tolerance in Da, defaulting to 5 ppm and 0 Da. The inclusive half-window is
#'   `tolerance + abs(precursorMz) * ppm / 1e6`.
#' @param spectra For a stored scan tree, supply the corresponding measured
#'   spectra. Acquisition metadata and directed links must agree. The
#'   supplied peak arrays define the measurements for this call.
#' @param passingOnly Use only products whose stored `pass_flag` equals one.
#'   Supply flagged spectra or the result of [filterFragSpectra()]. A parent
#'   peak retained in the source can supply precursor context even when
#'   filtering removed it from the products.
#' @details With a `purityA` object, the active assessment for each scan
#'   supplies its recorded parent, peak position and measured m/z. These
#'   must agree with the available parent measurement. A recorded unresolved
#'   selection keeps its status. If the recorded peak is no longer available,
#'   its products remain linked to an unmeasured selection. MS2 scans and
#'   older results without a peak observation use the mass window when the
#'   parent spectrum is present.
#'
#'   With measured `Spectra`, exactly one parent peak in the window supplies a measured
#'   precursor. Multiple candidates, no candidate or an unresolved parent
#'   produce an unmeasured selection with an explicit status. Purity scores
#'   describe the original isolation assessment. Known profile data and multiple isolated precursors
#'   per scan are unsupported. Subset large runs with
#'   [Spectra::filterPrecursorScan()] before reading peak arrays.
#'
#'   Peak IDs combine the scan lineage and the one-based peak position.
#'   Stored original positions take precedence. For queued processing they
#'   must be retained or recoverable. In-memory input without earlier
#'   position information uses the positions in the supplied measurement.
#' @return A `PrecursorFragmentTrees` list. Its `mst` holds the scan tree,
#'   including empty scans. Its `pft` holds peak nodes and fragmentation
#'   edges. Its `selections` table has one row per fragmentation scan,
#'   including scans with zero products. Its `assignment_source` column
#'   distinguishes `purityA` observations, `mz_window` matching and
#'   `unavailable` parents. Use [fragmentLinks()] for a table
#'   of products with measured masses and [subsetTree()] to follow a branch.
#' @examples
#' f <- system.file("extdata", "msn-standards",
#'     "msPurity_L_kynurenine_161166_MSn.mzML.gz", package = "msPurity")
#' s <- Spectra::Spectra(f, source = Spectra::MsBackendMzR(),
#'     BPPARAM = BiocParallel::SerialParam())
#' four <- s[match(c(966L, 977L, 985L, 987L), Spectra::acquisitionNum(s))]
#' trees <- precursorFragmentTrees(four)
#' products <- subsetTree(trees, scan = 985L, peak = 44L, mode = "out", depth = 1L)
#' fragmentLinks(products)
#' @export
precursorFragmentTrees <- function(object, ppm = 5, tolerance = 0,
        parentReference = "auto", spectra = NULL, passingOnly = FALSE) {
    for (value in list(ppm, tolerance))
        if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 1L ||
                !is.finite(value) || value < 0)
            stop("ppm and tolerance must be finite nonnegative numbers.", call. = FALSE)
    if (!is.logical(passingOnly) || length(passingOnly) != 1L || is.na(passingOnly))
        stop("passingOnly must be TRUE or FALSE.", call. = FALSE)

    input <- .pft_scan_input(object, spectra, if (missing(parentReference)) NULL else parentReference)
    mst <- input$mst
    spectra <- input$spectra
    scans <- mst$nodes
    events <- which(scans$msLevel > 1L)
    parents <- match(scans$parentLineageScanId[events], scans$lineageScanId)
    needed <- sort(unique(c(events, parents[!is.na(parents)])))
    source <- spectra[needed]
    vars <- Spectra::spectraVariables(source)
    if (any(Spectra::centroided(source) %in% FALSE) ||
            ("spectrum_representation" %in% vars &&
                any(source$spectrum_representation %in% "MS:1000128")))
        stop("Supply centroided spectra to construct a peak tree.", call. = FALSE)
    for (field in intersect(vars, c("n_precursors", "n_selected_ions",
            "nPrecursors", "nSelectedIons")))
        if (any(source[[field]] > 1, na.rm = TRUE))
            stop("Multiple isolated precursors per scan are unsupported.", call. = FALSE)

    peak_vars <- Spectra::peaksVariables(spectra)
    use_pa_flags <- passingOnly && methods::is(object, "purityA") &&
        length(object@filter_frag_params) > 0L
    if (passingOnly && !use_pa_flags && !"pass_flag" %in% peak_vars)
        stop("passingOnly requires a stored pass_flag or a result from filterFragSpectra().",
            call. = FALSE)
    columns <- c("mz", "intensity", intersect(c("original_peak_index", "pass_flag"), peak_vars))
    arrays <- if (length(needed)) as.list(Spectra::peaksData(source, columns = columns,
        BPPARAM = .msp_bpparam())) else list()
    positions <- if (length(needed)) .pa_filter_source_peaks(source, arrays)$original_positions
        else list()
    counts <- vapply(arrays, nrow, integer(1))
    scan_row <- rep(needed, counts)
    flat <- function(x) unlist(x, use.names = FALSE)
    mz <- as.numeric(flat(lapply(arrays, function(a) a[, "mz"])))
    intensity <- as.numeric(flat(lapply(arrays, function(a) a[, "intensity"])))
    if (any(!is.finite(mz) | mz <= 0) || any(!is.finite(intensity) | intensity < 0))
        stop("Peak m/z must be positive and intensities nonnegative and finite.", call. = FALSE)
    index <- as.integer(flat(positions))
    peak_id <- paste0("peak:", scans$lineageScanId[scan_row], ":", index, recycle0 = TRUE)
    products <- scans$msLevel[scan_row] > 1L
    flags <- if ("pass_flag" %in% peak_vars)
        as.numeric(flat(lapply(arrays, function(a) a[, "pass_flag"]))) else NULL
    if (use_pa_flags)
        flags <- as.numeric(flat(.pft_pa_flags(object, source, arrays, positions)))
    if (passingOnly) products <- products & flags %in% 1

    targets <- scans$precursorMz[events]
    recorded <- NULL
    if (methods::is(object, "purityA"))
        recorded <- .pft_pa_observations(object, spectra, scans, events)
    assignment <- .pft_assign_precursors(arrays, positions, counts, match(parents, needed),
        scans[events, , drop = FALSE], recorded, ppm, tolerance)
    selected <- assignment$row
    status <- assignment$status
    matched <- !is.na(selected)
    source_id <- peak_id[selected]
    source_id[!matched] <- paste0("selected:", scans$lineageScanId[events[!matched]],
        recycle0 = TRUE)
    selections <- data.frame(lineageScanId = scans$lineageScanId[events],
        run_id = scans$run_id[events],
        acquisitionNum = scans$acquisitionNum[events], parentLineageScanId = scans$parentLineageScanId[events],
        rawParentReference = scans$rawParentReference[events], isolationMz = targets,
        source = source_id, parent_peak_index = index[selected], measured_mz = mz[selected],
        match_status = status, assignment_source = assignment$source, stringsAsFactors = FALSE)

    # Include every eligible product and the parent peaks used by those acquisitions.
    keep <- products
    keep[selected[matched]] <- TRUE
    unknown <- events[!matched]
    rows <- c(scan_row[keep], unknown)
    nodes <- data.frame(node_id = c(peak_id[keep], source_id[!matched]),
        lineageScanId = scans$lineageScanId[rows], acquisitionNum = scans$acquisitionNum[rows],
        run_id = scans$run_id[rows], msLevel = scans$msLevel[rows], rtime = scans$rtime[rows],
        peak_index = c(index[keep], rep(NA_integer_, length(unknown))),
        mz = c(mz[keep], rep(NA_real_, length(unknown))),
        intensity = c(intensity[keep], rep(NA_real_, length(unknown))),
        measured = c(rep(TRUE, sum(keep)), rep(FALSE, length(unknown))),
        context = c(!products[keep], rep(TRUE, length(unknown))),
        match_status = c(rep("measured", sum(keep)), status[!matched]),
        stringsAsFactors = FALSE)
    if (length(unknown)) {
        nodes$label <- rep(NA_character_, nrow(nodes))
        nodes$label[!nodes$measured] <- paste0("Scan ", scans$acquisitionNum[unknown],
            "\nSelected m/z ", format(targets[!matched], digits = 7),
            "\n", status[!matched])
    }
    if (!is.null(flags)) nodes$pass_flag <- c(flags[keep], rep(NA_real_, length(unknown)))
    child <- match(scan_row[products], events)
    observation <- paste0("observation:", scans$lineageScanId[scan_row[products]],
        ":", index[products], recycle0 = TRUE)
    edges <- data.frame(source = source_id[child], target = peak_id[products],
        relationship = rep("fragmentation", length(child)), observation_id = observation,
        match_status = status[child], assignment_source = assignment$source[child],
        isolationMz = targets[child], stringsAsFactors = FALSE)
    metadata <- list(structure = "PFT", ppm = ppm, tolerance = tolerance,
        passingOnly = passingOnly, parentReference = mst$metadata$parentReference)
    pft <- structure(list(nodes = nodes, edges = edges, metadata = metadata),
        class = c("PFT", "SpectralTree"))
    structure(list(mst = mst, pft = pft, selections = selections, metadata = metadata),
        class = "PrecursorFragmentTrees")
}

.pft_scan_input <- function(object, spectra, parentReference) {
    if (inherits(object, "MassSpectralTrees")) {
        if (!methods::is(spectra, "Spectra"))
            stop("Supply the measured spectra for this scan tree.", call. = FALSE)
        mst <- object
        if (is.null(parentReference)) parentReference <- mst$metadata$requestedParentReference
        if (is.null(parentReference) && mst$metadata$parentReference %in%
                c("acquisitionNum", "scanIndex", "spectrumId", "sourceSpectrumId"))
            parentReference <- mst$metadata$parentReference
        if (is.null(parentReference)) parentReference <- "auto"
        current <- massSpectralTrees(spectra, parentReference)
        order <- match(mst$nodes$lineageScanId, current$nodes$lineageScanId)
        fields <- setdiff(names(current$nodes), "inPurity")
        same <- vapply(fields, function(field)
            identical(mst$nodes[[field]], current$nodes[[field]][order]), logical(1))
        if (length(spectra) != length(order) || anyNA(order) || !all(same))
            stop("The scan tree and supplied acquisition metadata disagree.", call. = FALSE)
        .tree_table_index(mst)
        expected <- current$edges[match(mst$edges$lineageScanId,
            current$edges$lineageScanId), , drop = FALSE]
        if (anyDuplicated(mst$edges$lineageScanId) || nrow(expected) != nrow(current$edges) ||
                !identical(mst$edges$parentLineageScanId, expected$parentLineageScanId) ||
                !identical(mst$edges$lineageScanId, expected$lineageScanId))
            stop("The scan tree and supplied acquisition edges disagree.", call. = FALSE)
        spectra <- spectra[order]
    } else {
        if (!is.null(spectra))
            stop("Supply spectra separately only with a stored scan tree.", call. = FALSE)
        if (is.null(parentReference)) parentReference <- "auto"
        mst <- massSpectralTrees(object, parentReference)
        spectra <- if (methods::is(object, "purityA"))
            object@spectra[match(mst$nodes$pid, object@spectra$pid)] else object
    }
    list(mst = mst, spectra = spectra)
}

# Product flags belong to the filtered output. Keep the source arrays for
# parent measurements, including peaks removed from that output.
.pft_pa_flags <- function(pa, source, arrays, positions) {
    filtered <- if (isTRUE(pa@filter_frag_params$allfrag))
        allFragSpectra(pa, legacy = FALSE) else groupedSpectra(pa, legacy = FALSE)
    flags <- lapply(positions, function(pos) rep(NA_real_, length(pos)))
    if (!length(filtered)) return(flags)
    vars <- c("original_peak_index", "pass_flag")
    if (!all(vars %in% Spectra::peaksVariables(filtered)))
        stop("Filtered spectra need stored original peak positions and pass_flag.", call. = FALSE)
    if (anyNA(filtered$pid) || anyDuplicated(filtered$pid))
        stop("Filtered spectra need unique scan pids.", call. = FALSE)
    source_row <- match(filtered$pid, source$pid)
    if (anyNA(source_row))
        stop("A filtered scan is missing from the retained source.", call. = FALSE)
    values <- Spectra::peaksData(filtered, columns = c("mz", "intensity", vars),
        BPPARAM = .msp_bpparam())
    for (i in seq_along(values)) {
        j <- source_row[i]
        m <- values[[i]]
        row <- match(m[, "original_peak_index"], positions[[j]])
        if (anyNA(row) || anyDuplicated(row) || !isTRUE(all(
                m[, "mz"] == arrays[[j]][row, "mz"] &
                m[, "intensity"] == arrays[[j]][row, "intensity"])))
            stop("Filtered peak positions and measurements must agree with the source.", call. = FALSE)
        flags[[j]][row] <- m[, "pass_flag"]
    }
    flags
}

# A later purityA run can assess only some stages. Each child names the
# processing step that produced its own assessment.
.pft_pa_observations <- function(pa, spectra, scans, events) {
    if (!length(events) || !"msPurityProcessingId" %in% Spectra::spectraVariables(spectra))
        return(NULL)
    step <- spectra$msPurityProcessingId[events]
    records <- processingParameters(pa, "purityA")
    selected <- lapply(intersect(unique(step), names(records)), function(id) {
        obs <- records[[id]]$msn_purity$observations
        if (!NROW(obs) || !"parent_peak_index" %in% names(obs))
            return(NULL)
        if (anyDuplicated(obs$childLineageScanId))
            stop("A purityA assessment contains duplicate child scan identities.", call. = FALSE)
        event <- which(step == id & scans$msLevel[events] > 2L)
        row <- match(scans$lineageScanId[events[event]], obs$childLineageScanId)
        keep <- !is.na(row)
        obs <- obs[row[keep], c("parentLineageScanId", "parent_peak_index", "matched_mz",
            "status"), drop = FALSE]
        obs$event <- event[keep]
        obs
    })
    do.call(rbind, selected)
}

.pft_assign_precursors <- function(arrays, positions, counts, parent_array, events,
        recorded, ppm, tolerance) {
    offsets <- cumsum(as.double(c(0L, counts)))
    selected <- rep(NA_real_, nrow(events))
    status <- events$parentStatus
    source <- rep("mz_window", nrow(events))
    source[is.na(parent_array)] <- "unavailable"

    if (NROW(recorded)) {
        source[recorded$event] <- "purityA"
        status[recorded$event] <- recorded$status
        for (k in which(!is.na(recorded$parent_peak_index))) {
            i <- recorded$event[k]
            j <- parent_array[i]
            status[i] <- "recorded_peak_unavailable"
            if (is.na(j) || !isTRUE(recorded$parentLineageScanId[k] == events$parentLineageScanId[i]))
                next
            # The recorded row is in the assessment's peak array. Filtering
            # can move that row while retaining its original position.
            pos <- recorded$parent_peak_index[k]
            rows <- unique(c(pos, match(pos, positions[[j]])))
            rows <- rows[!is.na(rows) & rows >= 1L & rows <= counts[j]]
            hits <- rows[which(arrays[[j]][rows, "mz"] == recorded$matched_mz[k])]
            if (length(hits) == 1L) {
                selected[i] <- offsets[j] + hits
                status[i] <- "matched"
            }
        }
    }

    for (i in which(source == "mz_window")) {
        j <- parent_array[i]
        target <- events$precursorMz[i]
        if (!is.finite(target) || target <= 0) {
            status[i] <- "missing_precursor_mz"
            next
        }
        window <- tolerance + abs(target) * ppm / 1e6
        hits <- which(abs(arrays[[j]][, "mz"] - target) <= window)
        if (length(hits) == 1L) selected[i] <- offsets[j] + hits
        status[i] <- if (!length(hits)) "unmatched" else if (length(hits) == 1L)
            "matched" else "ambiguous"
    }
    list(row = selected, status = status, source = source)
}
