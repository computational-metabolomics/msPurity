# Complete filtering parameters and acquisition-metadata selection.

.PA_FILTER_AXES <- list(
    msLevel = c("msLevel", "ms_level"),
    massAnalyzer = c("massAnalyzer", "analyzer", "analyzerType",
        "massAnalyzerType", "analyzerName"),
    scanFilterAnalyzer = "scanFilterAnalyzer",
    activationType = c("activationType", "activationMethod"),
    collisionEnergy = c("collisionEnergy", "collision_energy"),
    collisionEnergyUnit = c("collisionEnergyUnit", "collisionEnergyUnits",
        "collisionEnergyUnitName"),
    polarity = "polarity"
)

.PA_ACTIVATION_TYPES <- c(
    CID = "CID", "collision-induced dissociation" = "CID",
    "MS:1000133" = "CID",
    HCD = "HCD", "beam-type collision-induced dissociation" = "HCD",
    "MS:1000422" = "HCD",
    ETD = "ETD", "electron transfer dissociation" = "ETD",
    "MS:1000598" = "ETD")

# Normalize only recognized single methods, never parts of combined methods.
.pa_activation_type <- function(values) {
    values <- as.character(values)
    index <- match(tolower(trimws(values)),
        tolower(names(.PA_ACTIVATION_TYPES)))
    known <- !is.na(index)
    values[known] <- unname(.PA_ACTIVATION_TYPES[index[known]])
    values
}

.pa_filter_named <- function(x, label, empty = TRUE) {
    if (!is.list(x) || (!empty && !length(x)) ||
            (length(x) && (is.null(names(x)) || anyNA(names(x)) ||
                        any(!nzchar(names(x))) || anyDuplicated(names(x)))))
        stop(label, " must be a list with unique, nonempty names.")
    x
}

.pa_filter_select <- function(select) {
    select <- .pa_filter_named(select, "select")
    if (any(!names(select) %in% names(.PA_FILTER_AXES)))
        stop("select contains unsupported acquisition dimensions.")
    for (name in names(select)) {
        value <- select[[name]]
        numeric_axis <- name %in% c("msLevel", "collisionEnergy", "polarity")
        unknown <- is.logical(value) && length(value) == 1L && is.na(value)
        if (!length(value) || (!unknown && if (numeric_axis)
                    !is.numeric(value) else !is.character(value)))
            stop("select$", name, " has an invalid selector type.")
        if (numeric_axis && (any(is.nan(value)) ||
                    any(!is.finite(value[!is.na(value)]))))
            stop("select$", name, " must contain finite values or NA.")
        known <- value[!is.na(value)]
        if (name == "msLevel" && any(known < 1 | known != trunc(known)))
            stop("select$msLevel must contain positive integers or NA.")
        if (name == "collisionEnergy" && any(known < 0))
            stop("select$collisionEnergy must be nonnegative or NA.")
        if (name == "polarity" && any(!known %in% c(0, 1)))
            stop("select$polarity must contain 0, 1 or NA.")
        if (is.character(value) && any(!nzchar(known)))
            stop("select$", name, " must use NA, not empty text, for unknown.")
    }
    select
}

.pa_filter_values <- function(values, allow_skip = TRUE) {
    for (name in c("ilim", "plim", "ra", "snr")) {
        value <- values[[name]]
        skip <- name == "plim" && allow_skip && is.numeric(value) &&
            length(value) == 1L && is.na(value) && !is.nan(value)
        if (!skip && (!is.numeric(value) || length(value) != 1L ||
                    !is.finite(value) || value < 0))
            stop(name, " must be one finite, nonnegative number.")
        if (name == "plim" && !skip && value > 1)
            stop("plim must be between zero and one, ",
                "or NA in a parameter object.")
    }
    if (!is.character(values$snmeth) || length(values$snmeth) != 1L ||
            anyNA(values$snmeth) || !values$snmeth %in% c("median", "mean"))
        stop("snmeth must be 'median' or 'mean'.")
    values
}

.pa_filter_param_values <- function(object) {
    list(select = object@select, ilim = object@ilim, plim = object@plim,
        ra = object@ra, snr = object@snr, snmeth = object@snmeth)
}

#' Complete parameters for fragmentation-spectrum filtering
#'
#' A complete filtering parameter set extending \code{ProtGenerics::Param}.
#' Defaults are fixed when the object is constructed. Supply one object or a
#' named list of objects to \code{filterFragSpectra(param = )}. Scalar calls
#' remain available. Call-wide \code{rmp}, \code{allfrag} and \code{cores}
#' are not fields of this object; parameter-set filtering requires one core.
#'
#' \code{select} assigns settings to spectra; it does not discard spectra.
#' Fields are combined with AND and values within a field with OR. Supported
#' fields are \code{msLevel}, \code{massAnalyzer}, \code{scanFilterAnalyzer},
#' \code{activationType}, \code{collisionEnergy}, \code{collisionEnergyUnit}
#' and \code{polarity}. CID, HCD and ETD can be selected by abbreviation,
#' recognized CV name or accession. Source labels and accessions are retained.
#' Available structured activation metadata take precedence over scan filters;
#' contradictory known labels and accessions are rejected. Combined or
#' unrecognized methods are not reduced to individual methods.
#' When metadata are unknown, ordinary serial scan filters
#' can supply the final CID/HCD activation type and an \code{FTMS} or
#' \code{ITMS} analyzer category. This category is stored as
#' \code{scanFilterAnalyzer}; it does not identify the exact
#' \code{massAnalyzer}. Existing known metadata take precedence. Unsupported
#' or ambiguous filters leave the activation type unknown. Collision energies
#' and their units are not inferred from scan filters.
#' Use \code{NA} to select unknown metadata. Spectra matched by no parameter
#' set receive the filtering call's scalar settings. Overlapping sets must
#' have the same complete effective settings.
#'
#' @param select Named list of exact acquisition-metadata selectors.
#' @param ilim Minimum peak intensity in recorded units.
#' @param plim Minimum MS2 precursor purity, between zero and one. \code{NA}
#'   explicitly skips purity filtering for spectra matched by this object.
#' @param ra Minimum relative abundance on the existing 0 to 100 scale.
#' @param snr Minimum signal-to-noise ratio, using the existing statistic.
#' @param snmeth Noise statistic: \code{"median"} or \code{"mean"}.
#' @return A \code{FilterFragSpectraParam} object. \code{as.list()} returns
#'   its complete settings, including the selectors.
#' @name FilterFragSpectraParam
#' @aliases FilterFragSpectraParam-class
#' @importClassesFrom ProtGenerics Param
#' @exportClass FilterFragSpectraParam
#' @seealso \code{\link{filterFragSpectra}}, \code{\link{processingHistory}}
#' @examples
#' filter_param <- FilterFragSpectraParam(plim = 0.7, snr = 3)
#' as.list(filter_param)
#' hcd_param <- FilterFragSpectraParam(
#'     select = list(activationType = "HCD"), snr = 2)
#' hcd_param
NULL

methods::setClass("FilterFragSpectraParam",
    contains = "Param",
    slots = c(select = "list", ilim = "numeric", plim = "numeric",
        ra = "numeric", snr = "numeric", snmeth = "character"),
    prototype = list(select = list(), ilim = 0, plim = 0.8, ra = 0,
        snr = 3, snmeth = "median"),
    validity = function(object) {
        tryCatch({
            .pa_filter_select(object@select)
            .pa_filter_values(.pa_filter_param_values(object))
            TRUE
        }, error = function(error) conditionMessage(error))
    }
)

#' @rdname FilterFragSpectraParam
#' @export
FilterFragSpectraParam <- function(select = list(), ilim = 0, plim = 0.8,
    ra = 0, snr = 3, snmeth = "median") {
    if (identical(plim, NA)) plim <- NA_real_
    methods::new("FilterFragSpectraParam", select = select, ilim = ilim,
        plim = plim, ra = ra, snr = snr, snmeth = snmeth)
}

#' @rdname FilterFragSpectraParam
#' @param x A filtering parameter object.
#' @param ... Ignored.
#' @export
methods::setMethod("as.list", "FilterFragSpectraParam", function(x, ...) {
    .pa_filter_param_values(x)
})

.pa_filter_inputs <- function(param) {
    if (methods::is(param, "FilterFragSpectraParam"))
        param <- list(default = param)
    param <- .pa_filter_named(param, "param", empty = FALSE)
    for (name in names(param)) {
        if (!methods::is(param[[name]], "FilterFragSpectraParam"))
            stop("Parameter set ", name, " requires FilterFragSpectraParam.")
        methods::validObject(param[[name]])
    }
    param
}

# Only ordinary serial filters are interpreted. A later unsupported event must
# not expose the CID/HCD method of an earlier precursor stage.
.pa_filter_string_metadata <- function(strings, ms_level) {
    result <- data.frame(
        scanFilterAnalyzer = rep(NA_character_, length(strings)),
        activationType = rep(NA_character_, length(strings)))
    number <- "(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
    event <- paste0("^", number, "@([[:alpha:]]+)", number, "$")
    range <- paste0("\\[\\s*", number, "\\s*-\\s*", number,
        "(?:\\s*,\\s*", number, "\\s*-\\s*", number, ")*\\s*\\]")
    for (i in seq_along(strings)) {
        text <- strings[i]
        if (is.na(text)) next
        prefix <- regmatches(text, regexpr("^(FTMS|ITMS)\\b", text,
                ignore.case = TRUE, perl = TRUE))
        if (length(prefix)) result$scanFilterAnalyzer[i] <- toupper(prefix)
        marker <- regexpr("\\bms[0-9]*\\b", text, ignore.case = TRUE,
            perl = TRUE)
        if (marker < 0L || is.na(ms_level[i])) next
        level <- sub("^ms", "", tolower(regmatches(text, marker)))
        level <- if (nzchar(level)) as.numeric(level) else 1
        if (level < 2L || level != ms_level[i]) next
        tail <- substring(text, marker + attr(marker, "match.length"))
        tail <- trimws(sub(paste0("\\s*", range, "\\s*$"), "", tail,
                perl = TRUE))
        events <- strsplit(tail, "[[:space:]]+")[[1L]]
        if (length(events) != level - 1L ||
                !all(grepl(event, events, ignore.case = TRUE, perl = TRUE)))
            next
        method <- toupper(sub(event, "\\1", events[length(events)],
                ignore.case = TRUE, perl = TRUE))
        if (method %in% c("CID", "HCD"))
            result$activationType[i] <- method
    }
    result
}

.pa_filter_metadata <- function(source) {
    n <- if (is.data.frame(source)) nrow(source) else length(source)
    available <- if (is.data.frame(source)) names(source) else
        Spectra::spectraVariables(source)
    values <- lapply(names(.PA_FILTER_AXES), function(axis) {
        numeric_axis <- axis %in% c("msLevel", "collisionEnergy", "polarity")
        out <- if (numeric_axis) rep(NA_real_, n) else rep(NA_character_, n)
        for (name in intersect(.PA_FILTER_AXES[[axis]], available)) {
            value <- source[[name]]
            if (numeric_axis) {
                if (!is.numeric(value) && !all(is.na(value)))
                    stop("Acquisition metadata ", name, " must be numeric.")
                value <- as.numeric(value)
                known <- value[!is.na(value)]
                if (any(!is.finite(known)))
                    stop("Acquisition metadata ", name,
                        " must be finite or NA.")
                if (axis == "msLevel" && any(known < 1 |
                            known != trunc(known)))
                    stop("Acquisition msLevel must be positive integer or NA.")
                if (axis == "collisionEnergy" && any(known < 0))
                    stop("Acquisition collisionEnergy must be nonnegative.")
                if (axis == "polarity" && any(!known %in% c(0, 1)))
                    stop("Acquisition polarity must be 0, 1 or NA.")
            } else {
                value <- as.character(value)
                value[!is.na(value) & (!nzchar(value) |
                            tolower(value) == "unknown")] <- NA_character_
            }
            compared_out <- if (axis == "activationType")
                .pa_activation_type(out) else out
            compared_value <- if (axis == "activationType")
                .pa_activation_type(value) else value
            if (any(!is.na(out) & !is.na(value) &
                        compared_out != compared_value))
                stop("Conflicting known acquisition aliases for ", axis, ".")
            take <- is.na(out) & !is.na(value)
            out[take] <- value[take]
        }
        if (axis %in% c("msLevel", "polarity")) as.integer(out) else out
    })
    names(values) <- names(.PA_FILTER_AXES)
    values <- as.data.frame(values, stringsAsFactors = FALSE)
    structured_activation <- rep(FALSE, n)
    if ("activationTypeAccession" %in% available) {
        accession <- as.character(source$activationTypeAccession)
        structured_activation <- !is.na(accession) & nzchar(trimws(accession))
        mapped <- .pa_activation_type(accession)
        label <- .pa_activation_type(values$activationType)
        methods <- unique(unname(.PA_ACTIVATION_TYPES))
        conflict <- label %in% methods & mapped %in% methods & label != mapped
        if (any(conflict))
            stop("Conflicting activation label and activationTypeAccession ",
                "at spectrum ", which(conflict)[1L], ".")
        take <- is.na(values$activationType) & mapped %in% methods
        values$activationType[take] <- mapped[take]
    }
    if ("filterString" %in% available) {
        derived <- .pa_filter_string_metadata(as.character(source$filterString),
            values$msLevel)
        for (field in names(derived)) {
            take <- is.na(values[[field]]) & !is.na(derived[[field]])
            if (field == "activationType")
                take <- take & !structured_activation
            values[[field]][take] <- derived[[field]][take]
        }
    }
    values
}

.pa_filter_call_controls <- function(fallback) {
    .pa_filter_values(fallback, allow_skip = FALSE)
    for (name in c("rmp", "allfrag")) {
        if (!is.logical(fallback[[name]]) || length(fallback[[name]]) != 1L ||
                anyNA(fallback[[name]]))
            stop(name, " must be TRUE or FALSE.")
    }
    if (!is.numeric(fallback$cores) || length(fallback$cores) != 1L ||
            !is.finite(fallback$cores) || fallback$cores != 1)
        stop("Parameter-set filtering currently supports cores = 1 only.")
    invisible(fallback)
}

.pa_filter_resolve <- function(source, param, fallback) {
    objects <- .pa_filter_inputs(param)
    .pa_filter_call_controls(fallback)
    metadata <- .pa_filter_metadata(source)
    specs <- lapply(objects, as.list)
    matched <- lapply(specs, function(spec) {
        keep <- rep(TRUE, length(source))
        for (axis in names(spec$select)) {
            wanted <- spec$select[[axis]]
            value <- metadata[[axis]]
            if (axis == "activationType") {
                wanted <- .pa_activation_type(wanted)
                value <- .pa_activation_type(value)
            }
            keep <- keep & ifelse(is.na(value), anyNA(wanted),
                value %in% wanted[!is.na(wanted)])
        }
        keep
    })
    settings <- sets <- list()
    assignment <- data.frame(spectrum_index = seq_along(source),
        parameter_id = rep(NA_character_, length(source)),
        fallback_used = rep(FALSE, length(source)), stringsAsFactors = FALSE)
    if ("pid" %in% Spectra::spectraVariables(source))
        assignment$pid <- source$pid
    assignment$parameter_set_names <- rep(list(character()), length(source))
    equal <- function(a, b) {
        isTRUE(all.equal(a, b, check.attributes = FALSE, tolerance = 0))
    }
    fields <- c("ilim", "plim", "ra", "snr", "snmeth")
    for (i in seq_along(source)) {
        ids <- names(matched)[vapply(matched, `[`, logical(1), i)]
        current <- if (length(ids)) specs[[ids[1L]]][fields]
        else fallback[fields]
        if (length(ids) > 1L && any(!vapply(specs[ids], function(spec) {
            equal(spec[fields], current)
        }, logical(1))))
            stop("Conflicting overlapping parameter sets at spectrum ", i, ".")
        same <- which(vapply(sets, equal, logical(1), current))
        if (!length(same)) {
            id <- sprintf("parameter-%04d", length(sets) + 1L)
            sets[[id]] <- current
        } else {
            id <- names(sets)[same[1L]]
        }
        settings[[i]] <- current
        assignment$parameter_id[i] <- id
        assignment$fallback_used[i] <- !length(ids)
        assignment$parameter_set_names[i] <- list(ids)
    }
    list(parameters = settings, assignments = assignment, metadata = metadata,
        parameter_sets = sets, parameter_objects = objects,
        unused_parameter_sets = names(matched)[
            !vapply(matched, any, logical(1))])
}
