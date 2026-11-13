# Error conditions and the codec gate for the Parquet backend.
#
# Errors are the typed mzStack conditions of MsBackendParquet: class
# `mzstack_<class>`, inheriting `mzstack_error`, so a handler catches the same
# failure whether msPurity or MsBackendParquet raised it. Conditions msPurity
# raises also carry `msPurity_parquet_<class>` and `msPurity_parquet_error`.

.MZS_ERROR_CLASSES <- c("format", "archive", "reference", "study", "stale",
                        "unsupported", "semantic", "capability", "resource")

`%||%` <- function(x, y) if (is.null(x)) y else x

# Columns named in arrow/dplyr expressions.
utils::globalVariables("spectrum_id_")

#' Add msPurity's classes to an mzStack condition.
#'
#' @noRd
.mzs_classed <- function(cond, class) {
    cl <- class(cond)
    class(cond) <- c(paste0("msPurity_parquet_", class), cl[1L],
                     "msPurity_parquet_error", cl[-1L])
    cond
}

#' Raise a Parquet error.
#'
#' @param class one of `.MZS_ERROR_CLASSES`.
#'
#' @param ... pasted together to form the message.
#'
#' @param data named list of fields added to the condition, e.g. the source
#'     key a `reference` error names.
#'
#' @noRd
.mzs_abort <- function(class, ..., data = list()) {
    class <- match.arg(class, .MZS_ERROR_CLASSES)
    stop(.mzs_classed(
        MsBackendParquet::mzstackCondition(class, ..., data = data), class))
}

#' Stop unless arrow was built with the zstd and snappy codecs, which the
#' results tables and spectra are written with.
#'
#' @noRd
.mzs_require <- function(what = "format = \"parquet\"") {
    codecs <- c("zstd", "snappy")
    bad <- codecs[!vapply(codecs, arrow::codec_is_available, logical(1))]
    if (length(bad))
        stop(what, " requires an arrow build with the ",
             paste(bad, collapse = " and "), " codec(s).", call. = FALSE)
    invisible(TRUE)
}

#' The installed version of a package, or `NA`.
#'
#' @noRd
.mzs_pkg_version <- function(pkg) {
    tryCatch(as.character(utils::packageVersion(pkg)),
             error = function(e) NA_character_)
}
