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
#' @section DIMS retirement:
#' The DIMS workflow is deprecated and will receive no further bug fixes or
#' feature development. There is no direct replacement in msPurity. Calls
#' remain functional during deprecation; the interfaces will become defunct
#' and be removed in subsequent release cycles.
#'
#' @description
#' This page documents the DIMS functions and the \code{purityD} class and
#' methods retained for compatibility. Shared functions used by the LC-MS
#' workflows remain supported. The \code{purityD} method for
#' \code{updatePeaks()} is also deprecated.
#'
#' @param fileList A data frame from \code{Getfiles()} with file paths and
#'   sample class information, or an \code{MsExperiment} whose sample data
#'   contains \code{filepth}, \code{name}, \code{sampleType} and \code{class}.
#'   Missing columns are filled from the spectra file names, with
#'   \code{sampleType = "sample"}.
#' @param sampleOnly Logical; apply filtering, grouping or purity prediction
#'   only to sample peak lists, excluding blanks.
#' @param mzML Logical; use mzML input when TRUE and CSV input when FALSE.
#' @return Return values depend on the function; see the individual sections.
#' @name msPurity-deprecated
NULL

.msp_deprecate_dims <- function(old) {
    .Deprecated(old = old, package = "msPurity", msg = paste0(
        "'", old, "' is deprecated as part of the DIMS workflow retirement. ",
        "DIMS functionality will receive no further bug fixes or feature ",
        "development. There is no direct replacement in msPurity. ",
        "See help('msPurity-deprecated', package = 'msPurity')."))
}
