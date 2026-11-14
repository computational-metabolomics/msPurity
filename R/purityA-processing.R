# Structured processing history shared by purityA and extracted Spectra.
# Records store references and settings rather than explicit peak-array copies.

.pa_processing_empty <- function() {
    list(version = 2L, records = list(), active = list())
}

.pa_processing_steps <- function() {
    c(
        "purityA", "frag4feature", "filterFragSpectra",
        "averageIntraFragSpectra", "averageInterFragSpectra",
        "averageAllFragSpectra"
    )
}

.pa_processing_validate <- function(history) {
    if (is.list(history) && identical(history$version, 1L))
        stop("Processing history version 1 is from an earlier development ",
            "trial and is not supported. Use a current object or reprocess ",
            "the original inputs.")
    if (!is.list(history) || !identical(history$version, 2L) ||
            !is.list(history$records) || !is.list(history$active))
        stop("Invalid msPurity processing history.")
    records <- history$records
    if (length(records) && (is.null(names(records)) ||
                anyNA(names(records)) || any(!nzchar(names(records))) ||
                anyDuplicated(names(records))))
        stop("Processing records must have unique step identifiers.")
    earlier <- character()
    for (id in names(records)) {
        record <- records[[id]]
        if (!is.list(record) || !identical(record$step_id, id) ||
                !is.character(record$step) || length(record$step) != 1L ||
                !record$step %in% .pa_processing_steps() ||
                !is.character(record$status) || length(record$status) != 1L ||
                !record$status %in% c("completed", "imported") ||
                !is.logical(record$complete) || length(record$complete) != 1L ||
                anyNA(record$complete) ||
                !is.list(record$software) || !is.list(record$time) ||
                !is.list(record$input) || !is.list(record$output) ||
                !is.list(record$parameters) ||
                (length(record$parameters) &&
                        (is.null(names(record$parameters)) ||
                                anyNA(names(record$parameters)) ||
                                any(!nzchar(names(record$parameters))) ||
                                anyDuplicated(names(record$parameters)))))
            stop("Invalid processing record: ", id)
        if (!is.character(record$input$previous_steps) ||
                !all(record$input$previous_steps %in% earlier))
            stop("Processing history refers to an unknown or later input step.")
        earlier <- c(earlier, id)
    }
    active <- history$active
    if (length(active) && (is.null(names(active)) || anyNA(names(active)) ||
                any(!names(active) %in% .pa_processing_steps()) ||
                anyDuplicated(names(active))))
        stop("Active processing settings must name processing steps.")
    for (step in names(active)) {
        id <- active[[step]]
        if (!is.character(id) || length(id) != 1L || is.na(id) ||
                !id %in% names(records) || records[[id]]$step != step)
            stop("Active processing settings refer to an invalid step.")
    }
    history
}

.pa_processing_refs <- function(sp) {
    cols <- intersect(
        c(
            "pid", "fileid", "seqNum", "spectrumId",
            "acquisitionNum", "scanIndex", "dataOrigin",
            "grpid", "av_level", "sample", "msPurityEffectiveParameterId"
        ),
        Spectra::spectraVariables(sp)
    )
    as.data.frame(Spectra::spectraData(sp, columns = cols))
}

.pa_processing_tips <- function(sp) {
    if (!"msPurityProcessingId" %in% Spectra::spectraVariables(sp))
        return(character())
    ids <- as.character(sp$msPurityProcessingId)
    unique(ids[!is.na(ids) & nzchar(ids)])
}

.pa_processing_input <- function(pa, step, allfrag = FALSE) {
    sp <- switch(step,
        frag4feature = pa@spectra,
        filterFragSpectra = pa@fragSpectra,
        averageInterFragSpectra = if (length(pa@avSpectra))
            pa@avSpectra[pa@avSpectra$av_level %in% "av_intra"] else
            pa@avSpectra,
        pa@fragSpectra
    )
    if (step == "filterFragSpectra" && allfrag)
        sp <- pa@spectra
    out <- list(
        previous_steps = .pa_processing_tips(sp),
        spectra = .pa_processing_refs(sp)
    )
    if (step == "filterFragSpectra" && allfrag)
        out$previous_steps <- union(
            out$previous_steps,
            .pa_processing_tips(pa@fragSpectra)
        )
    if (step != "purityA" && nrow(pa@grped_df))
        out$feature_links <- pa@grped_df[, intersect(
            c("grpid", "pid", "sample", "cid"), names(pa@grped_df)
        ),
        drop = FALSE
        ]
    out
}

.pa_processing_attach <- function(pa) {
    history <- .pa_processing_validate(pa@params$processing)
    for (slot in c("spectra", "fragSpectra", "avSpectra")) {
        sp <- methods::slot(pa, slot)
        meta <- .pa_spectra_metadata(sp)
        meta$msPurityProcessing <- history
        .pa_spectra_metadata(sp) <- meta
        methods::slot(pa, slot) <- sp
    }
    pa
}

.pa_processing_set_tip <- function(pa, step, id) {
    slot <- switch(step,
        purityA = "spectra",
        frag4feature = "fragSpectra",
        filterFragSpectra = "fragSpectra",
        "avSpectra"
    )
    sp <- methods::slot(pa, slot)
    if (!length(sp))
        return(pa)
    idx <- seq_along(sp)
    if (slot == "avSpectra") {
        level <- switch(step,
            averageIntraFragSpectra = "av_intra",
            averageInterFragSpectra = "av_inter",
            averageAllFragSpectra = "av_all"
        )
        idx <- which(sp$av_level == level)
    }
    ids <- if ("msPurityProcessingId" %in% Spectra::spectraVariables(sp))
        as.character(sp$msPurityProcessingId) else
        rep(NA_character_, length(sp))
    ids[idx] <- id
    sp$msPurityProcessingId <- ids
    methods::slot(pa, slot) <- sp
    pa
}

.pa_processing_import_id <- function(pa, step, parameters) {
    # An imported record identifies a known snapshot, not a new execution.
    # Keep the snapshot field name stable for imported-setting identities.
    snapshot <- list(stage = step, parameters = parameters,
        files = pa@fileList, spectra = .pa_processing_refs(pa@spectra))
    path <- tempfile("mspurity-history-")
    on.exit(unlink(path))
    saveRDS(snapshot, path, compress = FALSE, version = 2)
    paste0("imported-", step, "-", unname(tools::md5sum(path)))
}

.pa_processing_record <- function(
    pa, step, parameters, started, input,
    imported = FALSE
) {
    if (is.null(pa@params$processing))
        pa@params$processing <- .pa_processing_empty()
    id <- if (imported) .pa_processing_import_id(pa, step, parameters)
    else paste0("processing-", .mzs_uid())
    finished <- if (imported) as.POSIXct(NA) else Sys.time()
    if (imported)
        started <- as.POSIXct(NA)
    pa <- .pa_processing_set_tip(pa, step, id)
    slot <- switch(step,
        purityA = "spectra",
        frag4feature = "fragSpectra",
        filterFragSpectra = "fragSpectra",
        "avSpectra"
    )
    output <- .pa_processing_input(pa, step)
    output$previous_steps <- NULL
    sp <- methods::slot(pa, slot)
    selected <- if (length(sp)) sp[sp$msPurityProcessingId %in% id] else sp
    output$spectra <- .pa_processing_refs(selected)
    pa@params$processing$records[[id]] <- list(
        step_id = id, step = step, parameters = parameters,
        status = if (imported) "imported" else "completed",
        complete = !imported,
        software = list(package = "msPurity", version = if (imported)
                NA_character_ else
                as.character(utils::packageVersion("msPurity"))),
        time = list(
            started = started, finished = finished,
            elapsed_seconds = if (imported) NA_real_ else
                as.numeric(difftime(finished, started, units = "secs"))
        ),
        input = input, output = output
    )
    pa@params$processing$active[[step]] <- id
    .pa_processing_attach(pa)
}

.pa_processing_init <- function(pa) {
    if (!is.null(pa@params$processing)) {
        .pa_processing_validate(pa@params$processing)
        return(pa)
    }
    pa@params$processing <- .pa_processing_empty()
    known <- list(
        purityA = pa@params$purityA,
        frag4feature = pa@params$frag4feature,
        filterFragSpectra = pa@filter_frag_params,
        averageIntraFragSpectra = pa@av_intra_params,
        averageInterFragSpectra = pa@av_inter_params,
        averageAllFragSpectra = pa@av_all_params
    )
    for (step in names(known)) {
        if (!length(known[[step]]))
            next
        # Known settings are a snapshot, not evidence of an execution date
        # or a complete chain of dependencies in an older object.
        pa <- .pa_processing_record(pa, step, known[[step]],
            started = as.POSIXct(NA),
            input = list(previous_steps = character(), state = "unknown"),
            imported = TRUE
        )
    }
    .pa_processing_attach(pa)
}

.pa_processing_invalidate <- function(pa, levels = c("intra", "inter", "all")) {
    steps <- c(
        intra = "averageIntraFragSpectra",
        inter = "averageInterFragSpectra",
        all = "averageAllFragSpectra"
    )
    if (length(pa@avSpectra))
        pa@avSpectra <- pa@avSpectra[
            !pa@avSpectra$av_level %in% paste0("av_", levels)
        ]
    for (level in levels) {
        methods::slot(pa, paste0("av_", level, "_params")) <- list()
        pa@params$processing$active[[steps[[level]]]] <- NULL
        for (group in names(pa@av_spectra))
            pa@av_spectra[[group]][[paste0("av_", level)]] <- NULL
    }
    if (all(c("intra", "inter", "all") %in% levels)) {
        pa@av_spectra <- list()
        meta <- .pa_spectra_metadata(pa@avSpectra)
        meta$groups <- NULL
        meta$columns <- NULL
        .pa_spectra_metadata(pa@avSpectra) <- meta
    }
    pa
}

.pa_processing_relink <- function(pa) {
    pa <- .pa_processing_invalidate(pa)
    pa@filter_frag_params <- list()
    pa@all_frag_scans <- .pa_empty_frame()
    pa@params$processing$active$filterFragSpectra <- NULL
    if ("purity_pass_flag" %in% Spectra::spectraVariables(pa@spectra))
        pa@spectra$purity_pass_flag <- NULL
    meta <- .pa_spectra_metadata(pa@spectra)
    meta$puritydfColumns <- setdiff(meta$puritydfColumns, "purity_pass_flag")
    .pa_spectra_metadata(pa@spectra) <- meta
    pa@puritydf$purity_pass_flag <- NULL
    pa@grped_df$purity_pass_flag <- NULL
    pa
}

.pa_processing_select <- function(history, step) {
    if (!is.null(step) && (!is.character(step) || anyNA(step) ||
                any(!step %in% .pa_processing_steps())))
        stop("step must be NULL or msPurity processing function names.")
    records <- history$records
    if (!is.null(step))
        records <- records[vapply(
            records, function(x) x$step %in% step,
            logical(1)
        )]
    records
}

#' Processing history and effective settings
#'
#' Retrieve processing records from a purityA object or an
#' extracted Spectra object. Records are named by unique processing step
#' identifiers. Repeated calls retain earlier records. For Spectra,
#' only records reachable from the retained spectra are returned, including
#' their input steps. Subsetting therefore selects the corresponding history.
#' Concatenating subsets with the same processing history is supported;
#' histories from independently processed Spectra objects are not merged.
#' Retrieve their histories separately before combining spectra.
#'
#' Each record contains the processing operation as \code{step}, its unique
#' execution identifier as \code{step_id}, effective parameters (including
#' defaults), software version, execution times, and input/output references.
#' Records store settings and input/output references rather than explicit
#' copies of peak arrays. Function-valued settings remain functions and survive
#' RDS serialization. Their R environments are retained, not deep-copied, and
#' can contain other objects. Avoid capturing large or mutable analysis objects
#' in parameter functions. Older objects are imported with known settings
#' marked incomplete; dates and missing settings stay unknown.
#' Where worker settings differ, cores records the requested value and
#' effective_cores records the execution used. Imported settings do not
#' invent missing execution information.
#'
#' For new purity assessments, `parameters$resolved` contains one settings
#' list per input file, keyed by file position. Each list identifies the file,
#' the isolation offsets actually selected, their source and selection method,
#' and the `iwNormFun` function when normalization was applied. With
#' `iwNorm = FALSE`, the resolved function is `NULL` and its source is
#' `not_used`. Call arguments remain available at the parameter
#' list's top level and in the original phase-specific settings.
#' Sources distinguish caller overrides, recorded metadata and approximations.
#' When assessment cannot use these settings because MS1 or MS2 measurements
#' are absent, their source is `not_used` and the status explains why.
#'
#' Averaged Spectra retain `msPurityContributors` for the original feature-link
#' occurrences supporting their retained peaks. The
#' `msPurityClusterContributors` list maps averaged peak clusters to those
#' occurrences' `contributorIndex` values within the feature. These mappings
#' store occurrence references, not copies of peak arrays. Historical averaged
#' inputs without cluster mappings retain unknown fine-grained contributors.
#'
#' @param object A purityA or Spectra object.
#' @param step NULL for all steps, or processing function names such as
#'   "filterFragSpectra" and "averageAllFragSpectra".
#' @return A named list of records for processingHistory(), or a named list
#'   of their effective parameter lists for processingParameters().
#' @examples
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'     "2_frag4feature_pa.rds", package = "msPurity"))
#' pa <- filterFragSpectra(pa, plim = 0.7)
#' records <- processingHistory(pa, step = "filterFragSpectra")
#' records[[length(records)]]$parameters
#' processingParameters(groupedSpectra(pa, legacy = FALSE),
#'     step = "filterFragSpectra")
#' @name processingHistory
#' @aliases processingParameters
#' @aliases processingHistory,purityA-method
#' @aliases processingHistory,Spectra-method
#' @aliases processingParameters,purityA-method
#' @aliases processingParameters,Spectra-method
#' @export processingHistory processingParameters
NULL

setMethod("processingHistory", "purityA", function(object, step = NULL) {
    object <- .pa_update(object)
    .pa_processing_select(object@params$processing, step)
})

setMethod("processingHistory", "Spectra", function(object, step = NULL) {
    history <- .pa_spectra_metadata(object)$msPurityProcessing
    if (is.null(history))
        return(.pa_processing_select(.pa_processing_empty(), step))
    history <- .pa_processing_validate(history)
    keep <- .pa_processing_tips(object)
    if (!all(keep %in% names(history$records)))
        stop("Spectrum processing identifier is absent from its history.")
    repeat {
        parents <- unique(unlist(lapply(
            history$records[keep],
            function(x) x$input$previous_steps
        ), use.names = FALSE))
        expanded <- union(keep, parents)
        if (length(expanded) == length(keep))
            break
        keep <- expanded
    }
    history$records <- history$records[names(history$records) %in% keep]
    .pa_processing_select(history, step)
})

for (.pa_history_class in c("purityA", "Spectra")) {
    setMethod(
        "processingParameters", .pa_history_class,
        function(object, step = NULL) {
            lapply(processingHistory(object, step), function(x) x$parameters)
        }
    )
}
rm(.pa_history_class)
