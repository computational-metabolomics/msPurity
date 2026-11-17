# Assess selected MS2 events with their resolved settings, retaining the full
# file for reference lookup and interpolation. Read its peak arrays only once.
.pa_assess_purity_configured <- function(filepth, parameter_ids = NULL,
    parameter_sets = NULL, ...) {
    args <- list(...)
    if (!length(parameter_ids))
        return(do.call(.pa_assess_purity_single,
                c(list(filepth = filepth), args)))
    ids <- unique(parameter_ids)
    scans <- if (length(ids) > 1L && any(Spectra::msLevel(filepth) == 1L))
        getscans(filepth, args$mzRback) else NULL
    results <- lapply(ids, function(id) {
        settings <- parameter_sets[[id]]
        settings$not_applicable <- NULL
        effective <- utils::modifyList(args, settings, keep.null = TRUE)
        targets <- list(filepth = filepth,
            target_positions = which(parameter_ids == id), scans = scans)
        do.call(.pa_assess_purity_single, c(targets, effective))
    })
    names(results) <- ids
    if (length(results) == 1L) {
        result <- results[[1L]]
        result$resolved$parameter_id <- ids
        return(result)
    }
    # Each worker reports its actual windows and used MS1 measurements.
    # Do not imply that one scalar window or weighting applies to every event.
    combine <- function(tables, positions) {
        table <- do.call(rbind, tables)
        if (is.null(table)) return(NULL)
        table <- table[match(positions, table$seqNum), , drop = FALSE]
        rownames(table) <- NULL
        table
    }
    positions <- Spectra::scanIndex(filepth)[Spectra::msLevel(filepth) == 2L]
    purity <- combine(lapply(results, `[[`, "purity"), positions)
    resolved <- lapply(results, `[[`, "resolved")
    list(purity = purity, resolved = list(
        dataOrigin = Spectra::dataOrigin(filepth)[1L],
        status = "condition_specific", groups = resolved,
        windows = combine(lapply(resolved, `[[`, "windows"), positions),
        observations = combine(lapply(resolved, `[[`, "observations"),
            positions)))
}
