#' Tolerances for scoring one pair of fragmentation spectra
#'
#' Each side has its own closed interval centred on its m/z. Its half-width
#' is Da + m/z * ppm / 1000000. Intersecting intervals, including a single
#' shared endpoint, match. Product endpoints are rounded to ten decimal places
#' before comparison, as in the original peak aligner. Precursor endpoints
#' are compared at native precision at every fragmentation stage.
#'
#' @param q_ppmPrec,l_ppmPrec Query and reference precursor tolerances in ppm.
#' @param q_daPrec,l_daPrec Query and reference precursor tolerances in Da.
#' @param q_ppmProd,l_ppmProd Query and reference product tolerances in ppm.
#' @param q_daProd,l_daProd Query and reference product tolerances in Da.
#' @param collisionEnergyTolerance Maximum absolute difference between collision
#'   energies on the same scale, inclusive. Zero requires equal energies.
#' @return A \code{SpectralMatchingParam} object. All nine settings must be
#'   finite, nonnegative numeric scalars. \code{as.list()} returns their values.
#' @name SpectralMatchingParam
#' @aliases SpectralMatchingParam-class
#' @importClassesFrom ProtGenerics Param
#' @exportClass SpectralMatchingParam
#' @seealso \code{\link{scoreSpectrumPair}}
#' @examples
#' as.list(SpectralMatchingParam(q_daProd = 0.01))
NULL

.sm_defaults <- c(q_ppmPrec = 5, l_ppmPrec = 5, q_daPrec = 0, l_daPrec = 0,
    q_ppmProd = 10, l_ppmProd = 10, q_daProd = 0, l_daProd = 0,
    collisionEnergyTolerance = 0)

methods::setClass("SpectralMatchingParam", contains = "Param",
    slots = stats::setNames(rep("numeric", length(.sm_defaults)), names(.sm_defaults)),
    prototype = as.list(.sm_defaults),
    validity = function(object) {
        for (name in names(.sm_defaults)) {
            x <- methods::slot(object, name)
            if (length(x) != 1L || !is.null(dim(x)) || !is.finite(x) || x < 0)
                return(paste(name, "must be one finite nonnegative number."))
        }
        TRUE
    })

#' @rdname SpectralMatchingParam
#' @export
SpectralMatchingParam <- function(q_ppmPrec = 5, l_ppmPrec = 5,
    q_daPrec = 0, l_daPrec = 0, q_ppmProd = 10, l_ppmProd = 10,
    q_daProd = 0, l_daProd = 0, collisionEnergyTolerance = 0) {
    methods::new("SpectralMatchingParam", q_ppmPrec = q_ppmPrec,
        l_ppmPrec = l_ppmPrec, q_daPrec = q_daPrec, l_daPrec = l_daPrec,
        q_ppmProd = q_ppmProd, l_ppmProd = l_ppmProd, q_daProd = q_daProd,
        l_daProd = l_daProd, collisionEnergyTolerance = collisionEnergyTolerance)
}

#' @rdname SpectralMatchingParam
#' @param x A matching parameter object.
#' @param ... Ignored.
#' @export
methods::setMethod("as.list", "SpectralMatchingParam", function(x, ...) {
    stats::setNames(lapply(names(.sm_defaults), function(n) methods::slot(x, n)),
        names(.sm_defaults))
})

.sm_index <- function(x, i, side) {
    if (!methods::is(x, "Spectra"))
        stop(side, " must be a Spectra object.")
    if (!is.numeric(i) || length(i) != 1L || !is.null(dim(i)) ||
            !is.finite(i) || i != trunc(i) || i < 1 || i > length(x))
        stop(side, "Index must select one spectrum in ", side, ".")
}

.sm_scores <- c("dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent")

.sm_pair_peaks <- function(x) {
    p <- .msp_peaks(x)[[1L]]
    keep <- is.finite(p[, "mz"]) & p[, "mz"] > 0 &
        is.finite(p[, "intensity"]) & p[, "intensity"] >= 0
    data.frame(mz = p[keep, "mz"], i = p[keep, "intensity"], row.names = NULL)
}

#' Check compatibility and score a pair of MS2, MS3 or MS4 spectra
#'
#' Both entries must have the same supported MS level and overlapping,
#' positive immediate precursor masses. At each stage, known polarity,
#' activation, analyzer and collision energy are compared. The instrument
#' analyzer and scan-filter analyzer are independent fields. Energy values
#' are compared only with a shared known unit (e.g. eV or normalized energy).
#' Uniform controlled-vocabulary identifiers resolve different source labels.
#' An unrecognized identifier paired with free text leaves that field partial.
#' Different known units conflict. Unknown fields, including absent units or
#' missing parents in an external library, leave the comparison partial.
#' Disjoint known values conflict. Overlapping mixed contributor values leave
#' the comparison partial because individual contributors may differ.
#'
#' Measured parents are resolved inside each supplied object using same-run
#' identifiers. Ambiguous or unavailable parents leave the stage unknown.
#' Consensus entries use the precursor path and ancestor conditions saved by
#' msPurity. Each earlier stage back to MS2 is compared by MS level, so
#' different source files can have compatible paths. Keep parent entries in
#' the supplied objects and select the desired children with the index arguments.
#'
#' @param query,reference Measured or consensus \code{Spectra} objects.
#' @param param Tolerances from \code{\link{SpectralMatchingParam}}.
#' @param queryIndex,referenceIndex One-based positions in the supplied objects.
#' @return A list with \code{eligible}, \code{reason}, a \code{checks} data
#'   frame (stage, field, outcome), and \code{scores} named dpc, rdpc, cdpc,
#'   mcount, allcount and mpercent. A rejected pair reports the first failed
#'   requirement and has missing scores. Eligible pairs report compatible or
#'   compatible_partial. If either selected spectrum has no usable signal
#'   after queued processing, eligibility is retained, the reason is
#'   no_usable_peaks and scores are missing.
#' @seealso \code{\link{spectralMatching}}
#' @examples
#' d <- S4Vectors::DataFrame(msLevel = 2L, precursorMz = 200)
#' d$mz <- list(c(50, 80))
#' d$intensity <- list(c(100, 20))
#' sp <- Spectra::Spectra(d, source = Spectra::MsBackendMemory())
#' scoreSpectrumPair(sp, sp)$scores
#' @export
scoreSpectrumPair <- function(query, reference, param = SpectralMatchingParam(),
    queryIndex = 1L, referenceIndex = 1L) {
    if (!methods::is(param, "SpectralMatchingParam"))
        stop("param must be a SpectralMatchingParam object.")
    methods::validObject(param)
    .sm_index(query, queryIndex, "query")
    .sm_index(reference, referenceIndex, "reference")
    q <- query[queryIndex]
    l <- reference[referenceIndex]
    level <- Spectra::msLevel(q)
    levels <- c(level, Spectra::msLevel(l))
    reason <- if (anyNA(levels)) "ms_level_missing" else if (
        any(!levels %in% 2:4)) "ms_level_unsupported" else if (
        levels[1L] != levels[2L]) "ms_level_mismatch" else "compatible"
    checks <- data.frame(stage = level, field = "msLevel",
        outcome = if (reason == "compatible") "match" else reason)
    scores <- stats::setNames(rep(NA_real_, length(.sm_scores)), .sm_scores)
    result <- function(eligible) list(eligible = eligible, reason = reason,
        checks = checks, scores = scores)
    if (reason != "compatible")
        return(result(FALSE))

    mass_check <- function(a, b) {
        if (!length(a) || !length(b) || any(!is.finite(c(a, b))) ||
                any(c(a, b) <= 0)) return("missing")
        if (.sm_mass_overlap(a, b, param@q_ppmPrec, param@l_ppmPrec,
                param@q_daPrec, param@l_daPrec)) "match" else "mismatch"
    }
    mass <- mass_check(Spectra::precursorMz(q), Spectra::precursorMz(l))
    checks <- rbind(checks, data.frame(stage = level, field = "precursorMz",
        outcome = mass))
    if (mass != "match") {
        reason <- if (mass == "missing") "missing_precursor" else "precursor_mismatch"
        return(result(FALSE))
    }

    qs <- .sm_stages(query, queryIndex, level)
    ls <- .sm_stages(reference, referenceIndex, level)
    for (stage in seq.int(level, 2L)) {
        a <- qs[[stage]]
        b <- ls[[stage]]
        conditions <- .sm_condition_checks(a, b, param@collisionEnergyTolerance)
        if (stage < level)
            conditions <- c(ancestorStage = if (is.null(a) || is.null(b)) "missing" else "match",
                precursorMz = mass_check(a$precursorMz, b$precursorMz), conditions)
        checks <- rbind(checks, data.frame(stage = stage,
            field = names(conditions), outcome = unname(conditions)))
    }
    conflict <- which(checks$outcome == "mismatch")
    if (length(conflict)) {
        first <- checks[conflict[1L], ]
        field <- if (first$field == "precursorMz") "precursor" else first$field
        reason <- paste0(if (first$stage < level) "ancestor_", field, "_mismatch")
        return(result(FALSE))
    }
    if (any(checks$outcome != "match"))
        reason <- "compatible_partial"
    qp <- .sm_pair_peaks(q)
    lp <- .sm_pair_peaks(l)
    if (!any(qp$i > 0) || !any(lp$i > 0)) {
        reason <- "no_usable_peaks"
        return(result(TRUE))
    }
    scores <- .sm_score_peaks(qp, lp, param@q_ppmProd, param@l_ppmProd,
        raW = 0.5, mzW = 2, q_daProd = param@q_daProd, l_daProd = param@l_daProd)
    result(TRUE)
}
