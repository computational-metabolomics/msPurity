#' @include purityX-class.R
NULL

.px_fingerprint <- function(value) {
    path <- tempfile("mspurity-prediction-")
    on.exit(unlink(path))
    saveRDS(value, path, compress = FALSE, version = 2)
    list(version = 1L, algorithm = "sha256", value = .mzs_sha256(path))
}

.px_record <- function(object, prepared, calculations, settings, controls,
    started) {
    parameters <- c(settings, controls)
    parameters$resolved <- list(isotope_model_source =
            if (!settings$isotopes) "not_used" else if (is.null(settings$im))
                "pcalc_default" else "caller")
    input <- list(previous_steps = character(), reference = prepared$inputInfo,
        fingerprint = .px_fingerprint(list(reference = prepared$inputInfo,
                targets = prepared$targets, groups = prepared$groups,
                observations = calculations$observations)),
        observations = calculations$observations,
        excluded_peaks = calculations$excluded_peaks,
        exclusions = prepared$exclusions)
    key <- if ("grpid" %in% names(object@predictions)) "grpid" else "id"
    references <- data.frame(row = seq_len(nrow(object@predictions)),
        output_id = object@predictions[[key]])
    if (key == "grpid") {
        references$feature_id <- prepared$groups$feature_id[
            match(references$output_id, prepared$groups$grpid)]
    } else {
        references$peak_id <- prepared$targets$peak_id[
            match(references$output_id, prepared$targets$id)]
    }
    output <- list(prediction_rows = references,
        predictions = calculations$predictions)
    id <- paste0("processing-", .mzs_uid())
    record <- .msp_processing_record("purityX", id, parameters, started,
        input, output)
    record$software$implementation <- "msPurity::pcalc"
    record$software$dependencies <- lapply(c("xcms", "Spectra"), function(pkg) {
        list(package = pkg, version = as.character(utils::packageVersion(pkg)))
    })
    history <- .pa_processing_empty()
    history$records[[id]] <- record
    history$active$purityX <- id
    object@params <- list(purityX = parameters,
        processing = .pa_processing_validate(history))
    object
}

.px_import_history <- function(object) {
    history <- .pa_processing_empty()
    known <- list()
    for (name in c("purityType", "offsets", "cores", "fileignore")) {
        value <- methods::slot(object, name)
        if (length(value)) known[name] <- list(value)
    }
    if (!length(known) && !nrow(object@predictions)) return(history)
    identity <- .px_fingerprint(list(parameters = known,
            predictions = object@predictions))
    id <- paste0("imported-purityX-", identity$value)
    record <- .msp_processing_record("purityX", id, known, as.POSIXct(NA),
        input = list(previous_steps = character(), reference = NULL),
        output = list(prediction_rows = seq_len(nrow(object@predictions))),
        imported = TRUE)
    history$records[[id]] <- record
    history$active$purityX <- id
    history
}

#' Update saved feature-purity predictions
#'
#' Add processing-history storage to older purityX objects without recalculating
#' their predictions. Only settings stored in the saved object are imported.
#' Execution times, software versions, missing settings and input references
#' remain unknown, and the imported record is marked incomplete.
#'
#' @param object A purityX object.
#' @param ... Unused.
#' @param verbose Logical; unused.
#' @return An updated purityX object retaining its saved predictions.
#' @aliases updateObject,purityX-method
#' @seealso \code{\link{processingHistory}}
#' @export
methods::setMethod("updateObject", "purityX",
    function(object, ..., verbose = FALSE) {
        if (!methods::.hasSlot(object, "params")) {
            updated <- methods::new("purityX")
            for (name in setdiff(methods::slotNames("purityX"), "params")) {
                if (methods::.hasSlot(object, name))
                    methods::slot(updated, name) <- methods::slot(object, name)
            }
            object <- updated
        }
        if (is.null(object@params$processing)) {
            object@params$processing <- .px_import_history(object)
        } else {
            object@params$processing <- .pa_processing_validate(
                object@params$processing)
        }
        methods::validObject(object)
        object
    }
)

#' @rdname processingHistory
#' @details For purityX, each prediction call creates a new result with its own
#'   execution record. Its input references identify the MS1 observations and
#'   chromatographic peaks actually used. Predictions are feature summaries;
#'   they are not converted to Spectra. Retrieve histories separately when
#'   comparing results from independent calls. Saved older results expose only
#'   their known settings as incomplete historical records.
#' @aliases processingHistory,purityX-method
#' @export
methods::setMethod("processingHistory", "purityX",
    function(object, step = NULL) {
        object <- updateObject(object)
        .pa_processing_select(object@params$processing, step)
    }
)

#' @rdname processingHistory
#' @aliases processingParameters,purityX-method
#' @export
methods::setMethod("processingParameters", "purityX",
    function(object, step = NULL) {
        lapply(processingHistory(object, step), function(record) {
            record$parameters
        })
    }
)
