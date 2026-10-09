# Error conditions and package gates for the Parquet backend.
#
# Errors are R conditions of class `msPurity_parquet_<class>`, each also
# inheriting `msPurity_parquet_error`, so callers can handle one kind of
# failure without matching message text.
#
# arrow and jsonlite are suggested, not imported; every entry point checks
# for them first.

.MZS_ERROR_CLASSES <- c("format", "archive", "reference", "study", "stale",
                        "unsupported", "semantic", "capability", "resource")

`%||%` <- function(x, y) if (is.null(x)) y else x

# Columns named in arrow/dplyr expressions.
utils::globalVariables("spectrum_id_")

#' Raise a Parquet error.
#'
#' @param class one of `.MZS_ERROR_CLASSES`.
#'
#' @param ... pasted together to form the message.
#'
#' @param data named list of fields added to the condition, e.g. the source
#'     key a `Reference` error names.
#'
#' @noRd
.mzs_abort <- function(class, ..., data = list()) {
    class <- match.arg(class, .MZS_ERROR_CLASSES)
    cond <- structure(
        c(list(message = paste0(...), call = NULL), data),
        class = c(paste0("msPurity_parquet_", class), "msPurity_parquet_error",
                  "error", "condition"))
    stop(cond)
}

#' Stop unless arrow (with zstd and snappy) and jsonlite are installed.
#'
#' @noRd
.mzs_require <- function(what = "format = \"parquet\"") {
    need <- c("arrow", "jsonlite")
    miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss))
        stop(what, " requires the suggested package(s) ",
             paste(miss, collapse = ", "), ". Install them with ",
             "install.packages(c(", paste0("\"", miss, "\"", collapse = ", "),
             ")).", call. = FALSE)
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
