.pa_feature_values <- function(values) {
    if (!is.numeric(values$ppm) || length(values$ppm) != 1L ||
            !is.finite(values$ppm) || values$ppm < 0)
        stop("ppm must be one finite, nonnegative number.")
    plim <- values$plim
    if (!is.numeric(plim) || length(plim) != 1L || is.nan(plim) ||
            (!is.na(plim) && (!is.finite(plim) || plim < 0 || plim > 1)))
        stop("plim must be between zero and one, or NA.")
    for (name in c("intense", "convert2RawRT", "useGroup")) {
        value <- values[[name]]
        if (!is.logical(value) || length(value) != 1L || anyNA(value))
            stop(name, " must be TRUE or FALSE.")
    }
    reference <- values$parentReference
    if (!is.character(reference) || length(reference) != 1L ||
            anyNA(reference) || !reference %in% c("auto", "acquisitionNum",
            "scanIndex", "spectrumId", "sourceSpectrumId"))
        stop("parentReference must be 'auto', 'acquisitionNum', 'scanIndex', ",
            "'spectrumId' or 'sourceSpectrumId'.")
    invisible(values)
}

#' Complete parameters for feature association
#'
#' Supply one complete object to \code{frag4feature(param = )}. Explicit
#' scalar association arguments supplied alongside the object must agree with
#' its settings. The deprecated \code{use_group} argument must also agree.
#' Named lists and acquisition-metadata selectors are not supported.
#'
#' The call-wide \code{msLevel} argument selects feature-link output levels.
#' MS3/MS4 inherit the feature membership of their validated MS2 ancestors;
#' association retains the source spectra needed to resolve that ancestry.
#' Purity is assessed separately by \code{\link{purityA}}.
#'
#' @param ppm Mass tolerance in ppm for matching an MS2 precursor to an
#'   XCMS feature.
#' @param plim Minimum MS2 precursor purity, between zero and one, or
#'   \code{NA} to skip exclusion. A positive threshold retains MS2 links with
#'   purity strictly greater than the threshold. Descendant links are retained
#'   independently of this exclusion.
#' @param intense Use the most intense precursor for individual-peak
#'   association instead of the closest precursor to the isolation-window
#'   centre.
#' @param convert2RawRT Convert adjusted XCMS feature retention-time bounds
#'   to raw retention times.
#' @param useGroup Use the full grouped-feature retention-time bounds instead
#'   of individual chromatographic peaks.
#' @param parentReference Parent-reference interpretation used by
#'   \code{\link{frag4feature}}.
#' @return A \code{Frag4FeatureParam} extending \code{ProtGenerics::Param}.
#'   \code{as.list()} returns all six settings.
#' @name Frag4FeatureParam
#' @aliases Frag4FeatureParam-class
#' @importClassesFrom ProtGenerics Param
#' @exportClass Frag4FeatureParam
#' @seealso \code{\link{frag4feature}}, \code{\link{processingParameters}}
#' @examples
#' feature_param <- Frag4FeatureParam(ppm = 5, plim = 0.7)
#' as.list(feature_param)
#' feature_param
NULL

methods::setClass("Frag4FeatureParam",
    contains = "Param",
    slots = c(ppm = "numeric", plim = "numeric", intense = "logical",
        convert2RawRT = "logical", useGroup = "logical",
        parentReference = "character"),
    prototype = list(ppm = 5, plim = NA_real_, intense = TRUE,
        convert2RawRT = TRUE, useGroup = FALSE, parentReference = "auto"),
    validity = function(object) {
        tryCatch({
            .pa_feature_values(as.list(object))
            TRUE
        }, error = function(error) conditionMessage(error))
    }
)

#' @rdname Frag4FeatureParam
#' @export
Frag4FeatureParam <- function(ppm = 5, plim = NA_real_, intense = TRUE,
    convert2RawRT = TRUE, useGroup = FALSE, parentReference = "auto") {
    if (identical(plim, NA)) plim <- NA_real_
    methods::new("Frag4FeatureParam", ppm = ppm, plim = plim,
        intense = intense, convert2RawRT = convert2RawRT,
        useGroup = useGroup, parentReference = parentReference)
}

.pa_feature_resolve <- function(param, supplied) {
    if (!methods::is(param, "Frag4FeatureParam"))
        stop("param must be one Frag4FeatureParam object.")
    methods::validObject(param)
    settings <- as.list(param)
    for (name in names(supplied)) {
        field <- if (name == "use_group") "useGroup" else name
        value <- supplied[[name]]
        if (field == "plim" && identical(value, NA)) value <- NA_real_
        if (!isTRUE(all.equal(value, settings[[field]],
                    check.attributes = FALSE, tolerance = 0)))
            stop("Explicit ", name, " conflicts with param$", field, ".")
    }
    settings
}
