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
