# The mzStack capability set of msPurity, and the converter tables.
#
# Anything outside the capability set is refused rather than approximated.

.MZS_ENTITIES_WRITTEN <- c("feature", "chromatographic_peak",
                           "chromatographic_peak_feature", "abundance",
                           "assay", "sample", "spectrum_feature",
                           "merge_member", "compound", "compound_xref",
                           "compound_synonym", "evidence", "evidence_score",
                           "software", "coverage", "conversion",
                           "tool_provenance",
                           "source_identifier", "loss_ledger",
                           "x_mspurity_scan", "x_mspurity_scan_peak")

# Nothing is dropped as unmappable, so lossy_conversion is not declared.
.MZS_LOSSY <- FALSE

.MZS_SOURCE_FORMATS <- list(
    list(format = "msPurity createDatabase SQLite",
         versions = ">= 1.16.0",
         route = "mspurity:sqlite-file",
         coverage_manifest = list(
             fileinfo = "mapped", s_peak_meta = "mapped_partial",
             s_peaks = "mapped_partial", c_peak_X_s_peak_meta = "mapped",
             c_peak_group_X_s_peak_meta = "mapped", source = "mapped",
             c_peaks = "mapped", c_peak_groups = "mapped_partial",
             c_peak_X_c_peak_group = "mapped",
             metab_compound = "mapped", sm_matches = "mapped",
             l_s_peak_meta = "mapped_partial", xcms_match = "dropped",
             parameters = "not_present_in_source")),
    list(format = "msp2db SQLite library",
         versions = ">= 0.0.5",
         route = "library:msp2db",
         coverage_manifest = list(
             library_spectra_meta = "mapped", library_spectra = "mapped",
             library_spectra_source = "mapped",
             metab_compound = "mapped_partial",
             library_spectra_annotation = "refused")),
    list(format = "MSP (dialects: massbank, mona, mspurity)",
         versions = "any; the dialect is declared by the caller",
         route = "library:msp",
         coverage_manifest = list(records = "mapped", peaks = "mapped",
                                  unknown_keys = "refused")),
    list(format = "mzSpecLib",
         versions = "none",
         route = "library:mzspeclib",
         coverage_manifest = list(library = "refused")),
    list(format = "msPurity purityA",
         versions = ">= 1.37.4",
         route = "mspurity:live-object",
         coverage_manifest = list(
             puritydf = "mapped", grped_df = "mapped_partial",
             grped_ms2 = "mapped_partial", av_spectra = "mapped",
             parameters = "mapped")))

#' The mzStack capability set of msPurity
#'
#' @description
#'
#' Reports what msPurity's mzStack backend can do, as a capability set
#' conforming to `capabilities.schema.json`: the conformance classes it
#' claims, the run kinds it reads and writes, the results entities it writes,
#' the source formats it converts and how each error class is surfaced.
#'
#' The backend is experimental. A capability set is a declaration, not a
#' conformance claim.
#'
#' @return A named list; `jsonlite::toJSON(x, auto_unbox = TRUE)` gives the
#'     JSON document.
#'
#' @examples
#' str(mzstackCapabilities(), max.level = 1)
#' @export
mzstackCapabilities <- function() {
    errors <- stats::setNames(
        as.list(paste0("mzstack_", .MZS_ERROR_CLASSES)),
        c("Format", "Archive", "Reference", "Study", "Stale", "Unsupported",
          "Semantic", "Capability", "Resource"))
    results <- list(
        results_major_versions = .mzs_array(.MZS_RESULTS_MAJOR),
        entities_read = .mzs_array(.MZS_ENTITIES_WRITTEN),
        entities_written = .mzs_array(.MZS_ENTITIES_WRITTEN),
        reference_resolution = "by_id",
        source_freshness_enforced = TRUE,
        truncation = list(written = TRUE, honoured = TRUE),
        self_contained = FALSE)
    if (length(.MZS_SOURCE_FORMATS))
        results$converter <- list(
            source_formats = .MZS_SOURCE_FORMATS,
            run_identity = list(
                tuple = .mzs_array(c("source_checksum", "ordinal")),
                `function` = paste(
                    "run_id = 'run-' + first 12 hex digits of",
                    "SHA-256(<source_checksum, or the normalised path when",
                    "the file is unreadable> + newline + <1-based",
                    "ordinal>)")),
            lossy_conversion = .MZS_LOSSY)
    list(
        mzstack_version = .MZS_SPEC_VERSION,
        implementation = list(
            name = "msPurity",
            version = .mzs_pkg_version("msPurity"),
            language = "R",
            language_version = paste(R.version$major, R.version$minor,
                                     sep = "."),
            homepage = "https://github.com/computational-metabolomics/msPurity"),
        classes = .mzs_array(c("Reader", "Writer", "Results Reader",
                               "Results Writer")),
        reader = list(
            mzstack_major_versions = .mzs_array(
                .mzs_semver_major(.MZS_SPEC_VERSION)),
            run_kinds = .mzs_array("native"),
            signal_layouts = .mzs_array("list"),
            sample_metadata = FALSE,
            usi = FALSE,
            filters = list(),
            portable_form = "Apache Arrow (arrow::Table)",
            validate = TRUE,
            id_width = "int32"),
        writer = list(
            run_kinds = .mzs_array("native"),
            sample_metadata = FALSE),
        results = results,
        errors = errors)
}

# ---------------------------------------------------------------------------
# Converter tables carried by every conversion.
#
#   conversion         who converted what, from where, with which coverage
#   tool_provenance    upstream processing steps, with parameter completeness
#   source_identifier  the source's own identifier for every minted key
#   loss_ledger        one row per distinct loss, by disposition
#
# Anything not mapped in full has a ledger row saying why.
# ---------------------------------------------------------------------------

#' The `conversion` row.
#'
#' @param coverage named character, construct -> coverage state.
#'
#' @noRd
.mzs_conversion_row <- function(route, source_format, source_format_version,
                                coverage, source_uri = NA_character_,
                                checksum = NA_character_) {
    bad <- setdiff(coverage, .MZS_COVERAGE_STATES)
    if (length(bad))
        stop("Unknown coverage state(s): ", paste(bad, collapse = ", "),
             call. = FALSE)
    data.frame(
        conversion_id_ = 1L,
        converter_name = "msPurity",
        converter_version = .mzs_pkg_version("msPurity"),
        route = route,
        mzstack_version = .MZS_SPEC_VERSION,
        results_version = .MZS_RESULTS_VERSION,
        source_format = source_format,
        source_format_version = source_format_version,
        source_uri = source_uri,
        source_checksum = checksum,
        source_checksum_algorithm = if (is.na(checksum)) NA_character_
                                    else "sha256",
        converted_at = .mzs_now(),
        coverage_manifest = as.character(.mzs_to_json(as.list(coverage),
                                                      pretty = FALSE)),
        loss_ledger = "loss_ledger",
        stringsAsFactors = FALSE)
}

#' `tool_provenance` rows.
#'
#' @param steps list of `list(tool_name, tool_version, step_name, parameters,
#'     completeness)`.
#'
#' @noRd
.mzs_tool_provenance <- function(steps) {
    steps <- Filter(Negate(is.null), steps)
    if (!length(steps))
        return(NULL)
    data.frame(
        tool_provenance_id_ = seq_along(steps),
        tool_name = vapply(steps, function(s) s$tool_name, ""),
        tool_version = vapply(steps, function(s)
            as.character(s$tool_version %||% NA_character_), ""),
        step_name = vapply(steps, function(s) s$step_name, ""),
        parameters = vapply(steps, function(s)
            if (is.null(s$parameters)) NA_character_
            else as.character(.mzs_to_json(.mzs_json_safe(s$parameters),
                                           pretty = FALSE)), ""),
        parameters_completeness = vapply(steps, function(s)
            match.arg(s$completeness, .MZS_COMPLETENESS), ""),
        stringsAsFactors = FALSE)
}

#' `source_identifier` rows for one set of minted keys.
#'
#' @param value the source's own lexical form, unnormalised.
#'
#' @param type e.g. `"positional"`, `"name"`, `"path"`, `"native_id"`.
#'
#' @noRd
.mzs_source_ids <- function(target_table, target_key, source_table,
                            source_column, value, type) {
    n <- length(target_key)
    if (!n || is.null(value))
        return(NULL)
    data.frame(target_table = rep(target_table, n),
               target_key = as.numeric(target_key),
               source_table = rep(source_table, n),
               source_column = rep(source_column, n),
               source_value = as.character(value),
               source_value_type = rep(type, n),
               stringsAsFactors = FALSE)
}

#' Bind `source_identifier` rows and number them.
#'
#' @noRd
.mzs_bind_source_ids <- function(...) {
    df <- do.call(.mzs_rbind_fill, list(...))
    if (is.null(df))
        return(NULL)
    df <- df[!is.na(df$source_value), , drop = FALSE]
    cbind(source_identifier_id_ = seq_len(nrow(df)), df)
}

#' `loss_ledger` rows.
#'
#' @param losses list of `list(construct, disposition, affected_rows,
#'     target_ref, reason)`.
#'
#' @noRd
.mzs_loss_ledger <- function(losses) {
    losses <- Filter(Negate(is.null), losses)
    if (!length(losses))
        return(NULL)
    df <- data.frame(
        loss_ledger_id_ = seq_along(losses),
        source_construct = vapply(losses, function(l) l$construct, ""),
        disposition = vapply(losses, function(l) l$disposition, ""),
        affected_rows = vapply(losses, function(l)
            as.numeric(l$affected_rows %||% NA_real_), 0),
        target_ref = vapply(losses, function(l)
            as.character(l$target_ref %||% NA_character_), ""),
        reason = vapply(losses, function(l)
            as.character(l$reason %||% NA_character_), ""),
        stringsAsFactors = FALSE)
    bad <- setdiff(df$disposition, .MZS_DISPOSITIONS)
    if (length(bad))
        stop("Unknown disposition(s): ", paste(bad, collapse = ", "),
             call. = FALSE)
    df
}

#' Ledger rows for the specified columns a conversion leaves wholly null.
#'
#' @param reasons named list, `"<table>.<column>"` -> `list(disposition,
#'     reason)` for columns whose emptiness has a specific cause; the rest
#'     are concepts the source format does not have.
#'
#' @noRd
.mzs_null_column_losses <- function(tables, reasons = list()) {
    out <- list()
    for (name in names(tables)) {
        df <- tables[[name]]
        if (is.null(df) || !nrow(df))
            next
        cols <- .mzs_table_schema(name)$columns
        for (c in setdiff(cols$name, "activity")) {
            v <- df[[c]]
            if (!is.null(v) && !all(is.na(v) & !is.nan(v)))
                next
            ref <- paste0(name, ".", c)
            r <- reasons[[ref]] %||% reasons[[paste0("*.", c)]] %||%
                list(disposition = "no_such_concept_in_source",
                     reason = "msPurity does not compute this quantity.")
            out[[length(out) + 1L]] <- list(
                construct = ref, disposition = r$disposition,
                affected_rows = nrow(df), target_ref = ref,
                reason = r$reason)
        }
    }
    out
}
