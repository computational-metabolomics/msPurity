.PX_PARAMETER_FIELDS <- c("purityType", "offsets", "iwNorm", "iwNormFun",
    "ilim", "isotopes", "im")

.px_parameter_values <- function(object) {
    stats::setNames(lapply(.PX_PARAMETER_FIELDS, function(name) {
        methods::slot(object, name)
    }), .PX_PARAMETER_FIELDS)
}

.px_validate_parameters <- function(values) {
    modes <- c("purityFWHMmedian", "purityFWmedian", "purityFWHMmean",
        "purityFWmean")
    if (!is.character(values$purityType) || length(values$purityType) != 1L ||
            anyNA(values$purityType) || !values$purityType %in% modes)
        stop("purityType must name one documented purity prediction mode.")
    offsets <- values$offsets
    if (!is.numeric(offsets) || length(offsets) != 2L ||
            any(!is.finite(offsets)) || any(offsets < 0) || sum(offsets) <= 0)
        stop("offsets must contain nonnegative upper and lower offsets, ",
            "with positive total width.")
    for (name in c("iwNorm", "isotopes")) {
        value <- values[[name]]
        if (!is.logical(value) || length(value) != 1L || anyNA(value))
            stop(name, " must be TRUE or FALSE.")
    }
    if (!is.null(values$iwNormFun) && !is.function(values$iwNormFun))
        stop("iwNormFun must be a function or NULL.")
    ilim <- values$ilim
    if (!is.numeric(ilim) || length(ilim) != 1L || !is.finite(ilim) ||
            ilim < 0 || ilim > 1)
        stop("ilim must be between zero and one.")
    im <- values$im
    if (!is.null(im) && (!is.matrix(im) || !is.numeric(im) ||
                nrow(im) < 1L || ncol(im) != 8L || any(!is.finite(im))))
        stop("im must be NULL or a finite numeric matrix with eight columns.")
    invisible(values)
}

#' Parameters for chromatographic-feature purity prediction
#'
#' Supply one complete parameter object to \code{purityX(param = )}.
#' Explicit scientific arguments supplied alongside it must agree with its
#' settings. With \code{param = NULL}, \code{purityX()} uses its scalar
#' arguments. Sample/feature selection, workers and output choices remain
#' arguments of \code{purityX()}, rather than settings in this object.
#'
#' Prediction uses MS1 measurements. The assumed isolation window describes
#' the surrounding signal that would be isolated with the feature; it does not
#' require an acquired fragmentation spectrum. Use separate calls to compare
#' alternative settings.
#'
#' @param purityType Summary of per-scan purity: \code{"purityFWHMmedian"},
#'   \code{"purityFWmedian"}, \code{"purityFWHMmean"} or
#'   \code{"purityFWmean"}. FWHM uses the chromatographic peak's full width
#'   at half maximum; FW uses its full detected boundaries.
#' @param offsets Upper and lower isolation offsets in Da, in that order.
#'   Thus \code{c(0.7, 0.3)} describes a window from target minus 0.3 to target
#'   plus 0.7. This preserves the order used by \code{purityX()}, which differs
#'   from the lower/upper order of \code{PurityAParam}.
#' @param iwNorm Apply isolation-efficiency weighting.
#' @param iwNormFun Weighting function receiving offsets from the isolation
#'   window midpoint, or \code{NULL} for the default Gaussian.
#' @param ilim Fraction of the target intensity used as the low-intensity
#'   cutoff, between zero and one.
#' @param isotopes Exclude isotope peaks using the selected isotope model.
#' @param im Isotope-model matrix accepted by \code{\link{pcalc}}, or
#'   \code{NULL} for the default model.
#' @return A \code{PurityXParam} extending \code{ProtGenerics::Param}.
#'   \code{as.list()} retains all seven settings, including NULLs and functions.
#' @name PurityXParam
#' @aliases PurityXParam-class
#' @importClassesFrom ProtGenerics Param
#' @exportClass PurityXParam
#' @seealso \code{\link{purityX}}, \code{\link{processingParameters}}
#' @examples
#' prediction_param <- PurityXParam(offsets = c(0.5, 0.5), ilim = 0.05)
#' as.list(prediction_param)
#' prediction_param
NULL

methods::setClass("PurityXParam", contains = "Param",
    slots = c(purityType = "character", offsets = "numeric",
        iwNorm = "logical", iwNormFun = "ANY", ilim = "numeric",
        isotopes = "logical", im = "ANY"),
    prototype = list(purityType = "purityFWHMmedian", offsets = c(0.5, 0.5),
        iwNorm = FALSE, iwNormFun = NULL, ilim = 0.05, isotopes = FALSE,
        im = NULL),
    validity = function(object) {
        tryCatch({
            .px_validate_parameters(.px_parameter_values(object))
            TRUE
        }, error = function(error) conditionMessage(error))
    }
)

#' @rdname PurityXParam
#' @export
PurityXParam <- function(purityType = "purityFWHMmedian", offsets = c(0.5, 0.5),
    iwNorm = FALSE, iwNormFun = NULL, ilim = 0.05, isotopes = FALSE,
    im = NULL) {
    methods::new("PurityXParam", purityType = purityType, offsets = offsets,
        iwNorm = iwNorm, iwNormFun = iwNormFun, ilim = ilim,
        isotopes = isotopes, im = im)
}

#' @rdname PurityXParam
#' @param x A purity-prediction parameter object.
#' @param ... Ignored.
#' @export
methods::setMethod("as.list", "PurityXParam", function(x, ...) {
    .px_parameter_values(x)
})

.px_resolve_parameters <- function(param, scalars, explicit = character()) {
    if (is.null(param)) {
        .px_validate_parameters(scalars)
        return(scalars[.PX_PARAMETER_FIELDS])
    }
    if (!methods::is(param, "PurityXParam"))
        stop("param must be one PurityXParam object.")
    methods::validObject(param)
    settings <- as.list(param)
    for (name in intersect(explicit, .PX_PARAMETER_FIELDS)) {
        if ((name == "im" && !identical(dim(scalars$im), dim(settings$im))) ||
                !isTRUE(all.equal(scalars[[name]], settings[[name]],
                        check.attributes = FALSE, tolerance = 0)))
            stop("Explicit ", name, " conflicts with param$", name, ".")
    }
    settings
}
