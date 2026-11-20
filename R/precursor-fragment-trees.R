#' Construct a precursor-fragment tree from measured spectra
#'
#' Match each fragmentation scan's recorded precursor m/z against peaks in
#' its acquisition parent and link the selection to its observed products.
#'
#' @inheritParams massSpectralTrees
#' @param object Measured `Spectra`, a `purityA` result (its retained source
#'   spectra), or a scan tree from [massSpectralTrees()].
#' @param ppm,tolerance Nonnegative relative tolerance in ppm and absolute
#'   tolerance in Da, defaulting to 5 ppm and 0 Da. The inclusive half-window is
#'   `tolerance + abs(precursorMz) * ppm / 1e6`.
#' @param spectra For a stored scan tree, supply the corresponding measured
#'   spectra. Acquisition metadata and directed links must agree. The
#'   supplied peak arrays define the measurements for this call.
#' @param passingOnly Use only products whose stored `pass_flag` equals one.
#'   Requires flagged spectra, e.g. `allFragSpectra(pa, legacy = FALSE)`
#'   after filtering. A retained parent peak may still supply measured
#'   precursor context when its own product flag failed.
#' @details Exactly one parent peak in the window supplies a measured
#'   precursor. Multiple candidates, no candidate or an unresolved parent
#'   produce an unmeasured selection with an explicit status. Purity scores
#'   describe the original isolation assessment. The peak match here uses the
#'   stated mass window. Known profile data and multiple isolated precursors
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
#'   including scans with zero products. Use [fragmentLinks()] for a table
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
    if (passingOnly && !"pass_flag" %in% peak_vars)
        stop("passingOnly requires a stored pass_flag peak variable.", call. = FALSE)
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
    if (passingOnly) products <- products & flags %in% 1

    targets <- scans$precursorMz[events]
    assignment <- .pft_assign_precursors(arrays, counts, match(parents, needed),
        targets, scans$parentStatus[events], ppm, tolerance)
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
        match_status = status, stringsAsFactors = FALSE)

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
        match_status = status[child], isolationMz = targets[child], stringsAsFactors = FALSE)
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

.pft_assign_precursors <- function(arrays, counts, parent_array, targets, status, ppm, tolerance) {
    # Array offsets give direct access to each resolved parent's peaks.
    offsets <- cumsum(as.double(c(0L, counts)))
    selected <- rep(NA_real_, length(targets))
    for (i in which(!is.na(parent_array))) {
        j <- parent_array[i]
        if (!is.finite(targets[i]) || targets[i] <= 0) {
            status[i] <- "missing_precursor_mz"
            next
        }
        window <- tolerance + abs(targets[i]) * ppm / 1e6
        hits <- which(abs(arrays[[j]][, "mz"] - targets[i]) <= window)
        if (length(hits) == 1L) selected[i] <- offsets[j] + hits
        status[i] <- if (!length(hits)) "unmatched" else if (length(hits) == 1L)
            "matched" else "ambiguous"
    }
    list(row = selected, status = status)
}
