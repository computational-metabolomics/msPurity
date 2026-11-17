.PA_PURITY_FIELDS <- c("mostIntense", "nearest", "offsets", "interpol",
    "iwNorm", "iwNormFun", "ilim", "isotopes", "im", "ppmInterp")

.pa_purity_values <- function(values) {
    for (name in c("mostIntense", "nearest", "iwNorm", "isotopes")) {
        value <- values[[name]]
        if (!is.logical(value) || length(value) != 1L || anyNA(value))
            stop(name, " must be TRUE or FALSE.")
    }
    .pa_automatic_offsets(values$offsets)
    if (!is.character(values$interpol) || length(values$interpol) != 1L ||
            anyNA(values$interpol) ||
            !values$interpol %in% c("none", "linear", "spline"))
        stop("interpol must be 'none', 'linear' or 'spline'.")
    for (name in c("ilim", "ppmInterp")) {
        value <- values[[name]]
        if (!is.numeric(value) || length(value) != 1L ||
                !is.finite(value) || value < 0)
            stop(name, " must be one finite, nonnegative number.")
    }
    if (values$ilim > 1) stop("ilim must be between zero and one.")
    if (!is.null(values$iwNormFun) && !is.function(values$iwNormFun))
        stop("iwNormFun must be a function or NULL.")
    im <- values$im
    if (!is.null(im) && (!is.matrix(im) || !is.numeric(im) ||
                nrow(im) < 1L || ncol(im) != 8L || any(!is.finite(im))))
        stop("im must be NULL or a finite numeric matrix with eight columns.")
    invisible(values)
}

# Preserve NULL optional fields that ProtGenerics' Param list conversion drops.
.pa_purity_param_values <- function(object) {
    list(select = object@select, mostIntense = object@mostIntense,
        nearest = object@nearest, offsets = object@offsets,
        interpol = object@interpol, iwNorm = object@iwNorm,
        iwNormFun = object@iwNormFun, ilim = object@ilim,
        isotopes = object@isotopes, im = object@im,
        ppmInterp = object@ppmInterp)
}

#' Complete parameters for precursor-purity assessment
#'
#' Supply one object or a named list to \code{purityA(param = )}. A matching
#' object supplies complete settings; unmatched MS2 events use the scalar
#' call settings. Unmatched MS3/MS4 events use the scalar \code{mostIntense}
#' and \code{offsets}, with no isotope or low-intensity exclusion.
#'
#' The object stores reusable scientific settings. The call's \code{msLevel}
#' selects the spectra to assess; \code{select} only assigns settings within
#' that selection. Worker counts, plotting and output destinations remain
#' call-wide arguments. Processing history records these controls alongside
#' the resolved settings; see \code{\link{processingParameters}}.
#'
#' Selectors assign settings to requested child events without discarding
#' their parent measurements. They use the acquisition fields and matching
#' rules described for \code{\link{FilterFragSpectraParam}}. Overlapping sets
#' must agree on the settings applicable to the matched MS level.
#'
#' MS3/MS4 use one validated parent and no interpolation or isolation-efficiency
#' weighting. Their effective settings mark MS2-only controls as not applicable.
#' A parameter object assigned to these levels must have \code{iwNorm = FALSE}
#' and \code{iwNormFun = NULL}. Other constructor defaults are fixed regardless
#' of selectors: an object assigned to higher levels requests its \code{ilim}
#' and \code{isotopes} settings. Isotope exclusion needs a model appropriate
#' for the measured parent spectrum.
#'
#' @param select Named list of exact acquisition-metadata selectors.
#' @param mostIntense Select the strongest in-window peak instead of the
#'   recorded-intensity/closest-mass target.
#' @param nearest Select the nearest permitted MS1 measurement; MS2 only.
#' @param offsets Lower and upper isolation offsets in Da, or \code{NA} to
#'   use each child's recorded bounds.
#' @param interpol MS2 interpolation: \code{"linear"}, \code{"none"}, or
#'   the legacy \code{"spline"} option. Linear is recommended for interpolation.
#' @param iwNorm Enable MS2 isolation-efficiency weighting.
#' @param iwNormFun MS2 weighting function receiving offsets from the window
#'   midpoint, or \code{NULL} for the default Gaussian.
#' @param ilim Fraction of the selected precursor intensity used as the
#'   low-intensity cutoff, between zero and one.
#' @param isotopes Exclude isotope peaks using the selected model.
#' @param im Isotope-model matrix accepted by \code{\link{pcalc}}, or
#'   \code{NULL} for its default model.
#' @param ppmInterp MS2 mass tolerance in ppm for tracking the selected peak
#'   across neighbouring MS1 measurements.
#' @return A \code{PurityAParam} extending \code{ProtGenerics::Param}.
#'   \code{as.list()} returns all settings, including selectors and functions.
#' @name PurityAParam
#' @aliases PurityAParam-class
#' @importClassesFrom ProtGenerics Param
#' @exportClass PurityAParam
#' @seealso \code{\link{purityA}}, \code{\link{processingParameters}}
#' @examples
#' ms2 <- PurityAParam(select = list(msLevel = 2L), interpol = "linear")
#' higher <- PurityAParam(select = list(msLevel = 3:4),
#'     ilim = 0, isotopes = FALSE)
#' as.list(higher)
#' ms2
NULL

methods::setClass("PurityAParam",
    contains = "Param",
    slots = c(select = "list", mostIntense = "logical", nearest = "logical",
        offsets = "numeric", interpol = "character", iwNorm = "logical",
        iwNormFun = "ANY", ilim = "numeric", isotopes = "logical",
        im = "ANY", ppmInterp = "numeric"),
    prototype = list(select = list(), mostIntense = FALSE, nearest = TRUE,
        offsets = NA_real_, interpol = "linear", iwNorm = FALSE,
        iwNormFun = NULL, ilim = 0.05, isotopes = TRUE, im = NULL,
        ppmInterp = 7),
    validity = function(object) {
        tryCatch({
            .pa_filter_select(object@select)
            .pa_purity_values(.pa_purity_param_values(object))
            TRUE
        }, error = function(error) conditionMessage(error))
    }
)

#' @rdname PurityAParam
#' @export
PurityAParam <- function(select = list(), mostIntense = FALSE, nearest = TRUE,
    offsets = NA_real_, interpol = "linear", iwNorm = FALSE, iwNormFun = NULL,
    ilim = 0.05, isotopes = TRUE, im = NULL, ppmInterp = 7) {
    if (is.logical(offsets) && all(is.na(offsets)))
        offsets <- as.numeric(offsets)
    methods::new("PurityAParam", select = select, mostIntense = mostIntense,
        nearest = nearest, offsets = offsets, interpol = interpol,
        iwNorm = iwNorm, iwNormFun = iwNormFun, ilim = ilim,
        isotopes = isotopes, im = im, ppmInterp = ppmInterp)
}

#' @rdname PurityAParam
#' @param x A purity parameter object.
#' @param ... Ignored.
#' @export
methods::setMethod("as.list", "PurityAParam", function(x, ...) {
    .pa_purity_param_values(x)
})

.pa_purity_inputs <- function(param) {
    if (methods::is(param, "PurityAParam")) param <- list(default = param)
    param <- .pa_filter_named(param, "param", empty = FALSE)
    for (name in names(param)) {
        if (!methods::is(param[[name]], "PurityAParam"))
            stop("Parameter set ", name, " requires PurityAParam.")
        methods::validObject(param[[name]])
    }
    param
}

.pa_purity_higher_values <- function(values) {
    fields <- c("nearest", "interpol", "iwNorm", "iwNormFun", "ppmInterp")
    values[fields] <- list(NA, "none", FALSE, NULL, NA_real_)
    values$not_applicable <- fields
    values
}

.pa_purity_resolve <- function(source, param, fallback) {
    objects <- .pa_purity_inputs(param)
    .pa_purity_values(fallback)
    metadata <- .pa_filter_metadata(source)
    if (anyNA(metadata$msLevel) || any(!metadata$msLevel %in% 2:4))
        stop("Purity parameter assignment requires requested MS2-MS4 events.")
    specs <- lapply(objects, as.list)
    matched <- .pa_parameter_matches(metadata, specs)
    higher <- metadata$msLevel > 2L
    for (i in seq_along(specs)) {
        if (any(matched[[i]] & higher) &&
                (specs[[i]]$iwNorm || !is.null(specs[[i]]$iwNormFun)))
            stop("Parameter set ", names(specs)[i],
                " requests isolation-efficiency weighting for MS3/MS4; ",
                "weighting is available only for MS2.")
    }
    settings <- lapply(specs, `[`, .PA_PURITY_FIELDS)
    higher_settings <- lapply(settings, .pa_purity_higher_values)
    higher_fallback <- fallback[.PA_PURITY_FIELDS]
    higher_fallback[c("ilim", "isotopes", "im")] <- list(0, FALSE, NULL)
    values <- c(settings, higher_settings,
        list(fallback[.PA_PURITY_FIELDS],
            .pa_purity_higher_values(higher_fallback)))
    matched <- c(lapply(matched, function(rows) rows & !higher),
        lapply(matched, function(rows) rows & higher))
    fallback_index <- ifelse(higher, length(values), length(values) - 1L)
    .pa_parameter_assign(source, metadata, objects, matched, values,
        fallback_index)
}
