# Complete parameters for acquisition-aware fragmentation-spectrum averaging.

.PA_AVERAGE_COMBINE <- c("activationType", "collisionEnergy")

.pa_average_param_values <- function(object) {
    list(select = object@select, combine = object@combine,
        minfrac = object@minfrac, minnum = object@minnum, ppm = object@ppm,
        tolerance = object@tolerance, snr = object@snr, ra = object@ra,
        av = object@av, sumi = object@sumi)
}

.pa_average_values <- function(values) {
    .pa_filter_select(values$select)
    combine <- values$combine
    if (!is.character(combine) || anyNA(combine) || anyDuplicated(combine) ||
            any(!combine %in% .PA_AVERAGE_COMBINE))
        stop("combine must contain unique values from ",
            paste(.PA_AVERAGE_COMBINE, collapse = " and "), ".")
    number <- function(name, lower, upper = Inf, inclusive = TRUE) {
        value <- values[[name]]
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
                value > upper || (if (inclusive) value < lower else
                    value <= lower))
            stop(name, " must be one finite number ",
                if (inclusive) ">= " else "> ", lower,
                if (is.finite(upper)) " and <= ", if (is.finite(upper)) upper,
                ".")
    }
    number("minfrac", 0, 1)
    number("minnum", 0)
    number("ppm", 0)
    number("tolerance", 0)
    if (values$ppm == 0 && values$tolerance == 0)
        stop("ppm or tolerance must be greater than 0.")
    number("snr", 0)
    number("ra", 0)
    if (!is.character(values$av) || length(values$av) != 1L ||
            !values$av %in% c("median", "mean"))
        stop("av must be \"median\" or \"mean\".")
    if (!is.logical(values$sumi) || length(values$sumi) != 1L ||
            is.na(values$sumi))
        stop("sumi must be TRUE or FALSE.")
    invisible(values)
}

#' Complete parameters for fragmentation-spectrum averaging
#'
#' \code{AverageFragSpectraParam} holds one complete set of settings for
#' \code{\link{averageAllFragSpectra}}, \code{\link{averageIntraFragSpectra}}
#' and \code{\link{averageInterFragSpectra}}. With a parameter object, averaging
#' keeps acquisition conditions and MS3/MS4 precursor paths apart, so that one
#' feature can have several averaged spectra.
#'
#' \code{select} chooses the feature-linked spectra this set averages, using
#' the acquisition fields and matching rules of
#' \code{\link{FilterFragSpectraParam}}: all named fields must match, any listed
#' value of a field matches, and \code{NA} matches an unknown value. Spectra
#' that no set selects are not averaged, and sets may select the same spectra.
#' An empty \code{select} selects every feature-linked spectrum.
#'
#' Spectra are averaged together only when they belong to the same feature and
#' MS level and share mass analyzer, scan-filter analyzer, activation type,
#' collision energy, energy unit and polarity. Activation types are compared
#' as for \code{select}, so names and accessions of one method match.
#' \code{combine} allows the activation type and/or collision energy of the
#' averaged stage to differ. It never relaxes the conditions of earlier stages:
#' MS3 and MS4 spectra are averaged together only when every earlier stage of
#' their precursor paths was acquired under the same conditions. Unknown values
#' match only unknown values.
#'
#' Two m/z values agree when they differ by at most
#' \code{tolerance + ppm * mz / 1e6}, where \code{mz} is their mean; the
#' absolute and relative parts add, as for \code{MsCoreUtils::closest()}. Each
#' averaged peak combines fragment peaks that all agree with one another. For
#' MS3 and MS4 spectra, the precursor m/z of each stage after MS2 must also
#' agree, for every pair of spectra in an average, so that a chain of small
#' differences cannot join two distant precursors. Repeated selections of one
#' precursor that differ by more than this tolerance form separate averages.
#' The MS2 precursor is not compared again: the association with the feature,
#' made by \code{\link{frag4feature}}, already decided it. Spectra are assigned
#' in order of file, scan number and feature link, each to the first
#' compatible average.
#'
#' Ion-trap m/z values vary by a similar absolute amount across the m/z range,
#' so a separate parameter set for ion-trap spectra, selected with
#' \code{scanFilterAnalyzer = "ITMS"} or \code{massAnalyzer}, can use
#' \code{tolerance} (for example 0.2 to 0.5) with \code{ppm = 0}. Parameter
#' sets selected by MS level or collision energy can likewise use their own
#' tolerances. A wide tolerance can join two peaks of the same spectrum into
#' one averaged peak; both then count towards \code{minnum} and
#' \code{minfrac}, and \code{sumi} adds both intensities.
#'
#' The other settings have the meaning of the corresponding arguments of
#' \code{\link{averageAllFragSpectra}}. A single object passed as
#' \code{param} is recorded under the parameter-set name \code{"default"}; a
#' named list uses its names.
#'
#' @param select named list of acquisition selectors, as for
#'   \code{\link{FilterFragSpectraParam}}.
#' @param combine character; any of \code{"activationType"} and
#'   \code{"collisionEnergy"}.
#' @param minfrac numeric; minimum ratio of the number of clustered peaks to
#'   the number of averaged spectra that have eligible peaks.
#' @param minnum numeric; minimum number of clustered peaks.
#' @param ppm numeric; relative m/z tolerance in ppm for fragment peaks and
#'   MS3/MS4 precursors, added to \code{tolerance}. \code{ppm = 0} requires
#'   \code{tolerance > 0}.
#' @param tolerance numeric; absolute m/z tolerance for fragment peaks and
#'   MS3/MS4 precursors, added to the \code{ppm} tolerance.
#' @param snr numeric; minimum signal-to-noise ratio of an averaged peak.
#' @param ra numeric; minimum relative abundance of an averaged peak.
#' @param av character; \code{"median"} or \code{"mean"}.
#' @param sumi logical; sum, rather than average, the peak intensities.
#' @return An \code{AverageFragSpectraParam} object. \code{as.list()} returns
#'   its complete settings as a named list.
#' @name AverageFragSpectraParam
#' @aliases AverageFragSpectraParam-class
#' @importClassesFrom ProtGenerics Param
#' @exportClass AverageFragSpectraParam
#' @seealso \code{\link{averageAllFragSpectra}},
#'   \code{\link{FilterFragSpectraParam}}
#' @examples
#' param <- AverageFragSpectraParam(select = list(msLevel = 3L),
#'     combine = "collisionEnergy", sumi = FALSE)
#' as.list(param)
#' alternatives <- list(
#'     hcd = AverageFragSpectraParam(select = list(activationType = "HCD")),
#'     cid = AverageFragSpectraParam(select = list(activationType = "CID")))
#' ## FTMS (for example Orbitrap) spectra with a relative and ITMS (ion-trap)
#' ## spectra with an absolute tolerance
#' by_analyzer <- list(
#'     ftms = AverageFragSpectraParam(
#'         select = list(scanFilterAnalyzer = "FTMS")),
#'     itms = AverageFragSpectraParam(
#'         select = list(scanFilterAnalyzer = "ITMS"), ppm = 0,
#'         tolerance = 0.3))
NULL

methods::setClass("AverageFragSpectraParam",
    contains = "Param",
    slots = c(select = "list", combine = "character", minfrac = "numeric",
        minnum = "numeric", ppm = "numeric", tolerance = "numeric",
        snr = "numeric", ra = "numeric", av = "character", sumi = "logical"),
    prototype = list(select = list(), combine = character(), minfrac = 0.5,
        minnum = 1, ppm = 5, tolerance = 0, snr = 0, ra = 0, av = "median",
        sumi = TRUE),
    validity = function(object) {
        tryCatch({
            .pa_average_values(.pa_average_param_values(object))
            TRUE
        }, error = function(error) conditionMessage(error))
    }
)

#' @rdname AverageFragSpectraParam
#' @export
AverageFragSpectraParam <- function(select = list(), combine = character(),
    minfrac = 0.5, minnum = 1, ppm = 5, tolerance = 0, snr = 0, ra = 0,
    av = "median", sumi = TRUE) {
    methods::new("AverageFragSpectraParam", select = select,
        combine = combine, minfrac = minfrac, minnum = minnum, ppm = ppm,
        tolerance = tolerance, snr = snr, ra = ra, av = av, sumi = sumi)
}

#' @rdname AverageFragSpectraParam
#' @param x An averaging parameter object.
#' @param ... Ignored.
#' @export
methods::setMethod("as.list", "AverageFragSpectraParam", function(x, ...) {
    .pa_average_param_values(x)
})

# One parameter object or a named list of them, each valid.
.pa_param_inputs <- function(param, class) {
    if (methods::is(param, class)) param <- list(default = param)
    param <- .pa_filter_named(param, "param", empty = FALSE)
    for (name in names(param)) {
        if (!methods::is(param[[name]], class))
            stop("Parameter set ", name, " requires ", class, ".")
        methods::validObject(param[[name]])
    }
    param
}

# Explicit scalar arguments must agree with every parameter set, as for
# frag4feature(); rmp and cores are call-wide controls.
.pa_average_resolve <- function(param, supplied, rmp) {
    param <- .pa_param_inputs(param, "AverageFragSpectraParam")
    settings <- lapply(param, as.list)
    for (name in names(supplied)) {
        for (set in names(settings)) {
            if (!isTRUE(all.equal(supplied[[name]], settings[[set]][[name]],
                        check.attributes = FALSE, tolerance = 0)))
                stop("Explicit ", name, " conflicts with param",
                    if (length(settings) > 1L) "$",
                    if (length(settings) > 1L) set, "$", name, ".")
        }
    }
    if (!is.logical(rmp) || length(rmp) != 1L || is.na(rmp))
        stop("rmp must be TRUE or FALSE.")
    settings
}
