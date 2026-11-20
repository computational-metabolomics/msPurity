# SQLite output is the default but frozen (bug fixes only). Each write raises
# one `deprecatedWarning` via `.Deprecated()` naming the Parquet alternative.

#' Warn that an SQLite write is deprecated.
#'
#' @param old the call that wrote SQLite, e.g. `"createDatabase()"`.
#'
#' @param alternative what to use instead; by default `format = "parquet"`.
#'
#' @noRd
.msp_deprecate_sqlite <- function(old, alternative = NULL) {
    alternative <- alternative %||% paste0(
        "Use format = \"parquet\" to write a Parquet dataset ",
        "(experimental). See ?createDatabase.")
    .Deprecated(msg = paste0(
        "Writing SQLite output with ", old, " is deprecated and will ",
        "receive bug fixes only. ", alternative))
}


#' Deprecated interfaces in msPurity
#'
#' These interfaces are retained for compatibility during deprecation.
#'
#' @section DIMS retirement:
#' The DIMS workflow is deprecated and will receive no further bug fixes or
#' feature development. There is no direct replacement in msPurity. Calls
#' remain functional during deprecation; the interfaces will become defunct
#' and be removed in subsequent release cycles.
#'
#' @details
#' The deprecated interfaces are \code{\link{purityD}},
#' \code{\link{Getfiles}}, \code{\link{averageSpectraSingle}},
#' \code{\link{groupPeaksEx}}, \code{\link{subtractMZ}},
#' \code{\link{dimsPredictPuritySingle}} and
#' \code{\link{get_additional_mzml_meta}}, and the \code{purityD} methods
#' for \code{averageSpectra}, \code{filterp}, \code{subtract},
#' \code{dimsPredictPurity}, \code{groupPeaks}, \code{writeOut},
#' \code{averagedPeaks}, \code{getP}, \code{updatePeaks},
#' \code{updateObject} and \code{show}.
#'
#' Shared functions used by the LC-MS workflows remain supported.
#' @name msPurity-deprecated
#' @keywords internal
NULL

.msp_deprecate_dims <- function(old) {
    .Deprecated(old = old, package = "msPurity", msg = paste0(
        "'", old, "' is deprecated as part of the DIMS workflow retirement. ",
        "DIMS functionality will receive no further bug fixes or feature ",
        "development. There is no direct replacement in msPurity. ",
        "See help('msPurity-deprecated', package = 'msPurity')."))
}
