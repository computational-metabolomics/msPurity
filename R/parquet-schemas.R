# The results tables as data: columns with Arrow types and nullability,
# total sort order, Hive partitioning and surrogate key. The writer casts to
# these schemas and the validator checks against them.
#
# Type codes, a trailing `!` marking NOT NULL:
#
#   i64 int64   i32 int32   f64 float64   f32 float32   bool boolean
#   str utf8    dict dictionary<int32, utf8>
#
# Every table gets an `activity` column and every spectrum reference all
# four reference columns. `spectrum_feature` and `coverage` sort orders end
# with tie-breakers so they stay total when a spectrum id is null. Tables
# with no mzStack counterpart are named `x_mspurity_<name>`.

#' The columns of a spectrum reference under one role prefix.
#'
#' @noRd
.mzs_ref_columns <- function(prefix, usi = FALSE) {
    cols <- c("dict!", "dict", "i64", "str", if (usi) "str")
    names(cols) <- paste0(prefix, c("_source", "_run_id", "_spectrum_id_",
                                    "_native_id", if (usi) "_usi"))
    cols
}

#' @param key the table's own integer surrogate key, or `NA` for a table
#'     whose rows are identified by the keys of others.
#'
#' @noRd
.mzs_table <- function(columns, sorted_by, partitioning = character(),
                       key = names(columns)[1L], extension = FALSE) {
    columns <- c(columns, activity = "dict!")
    data <- data.frame(name = names(columns),
                       type = sub("!$", "", unname(columns)),
                       nullable = !endsWith(unname(columns), "!"),
                       stringsAsFactors = FALSE)
    list(columns = data, sorted_by = sorted_by, partitioning = partitioning,
         key = key, extension = extension)
}

.MZS_TABLES <- list(
    feature = .mzs_table(c(
        feature_id_ = "i64!",
        exp_mass_to_charge = "f64!",
        mass_to_charge_min = "f64",
        mass_to_charge_max = "f64",
        retention_time_in_seconds = "f64",
        retention_time_in_seconds_start = "f64",
        retention_time_in_seconds_end = "f64",
        charge = "i32",
        adduct_ion = "dict",
        isotopomer = "dict",
        theoretical_neutral_mass = "f64",
        n_chromatographic_peaks = "i32",
        feature_name = "str"),
        sorted_by = c("exp_mass_to_charge", "feature_id_")),

    chromatographic_peak = .mzs_table(c(
        chromatographic_peak_id_ = "i64!",
        run_id = "dict!",
        peak_index_in_run = "i32!",
        exp_mass_to_charge = "f64!",
        mass_to_charge_min = "f64",
        mass_to_charge_max = "f64",
        retention_time_in_seconds = "f64!",
        retention_time_in_seconds_start = "f64",
        retention_time_in_seconds_end = "f64",
        retention_time_in_seconds_start_adjusted = "f64",
        retention_time_in_seconds_end_adjusted = "f64",
        signal_to_noise = "f32"),
        sorted_by = c("run_id", "retention_time_in_seconds",
                      "chromatographic_peak_id_"),
        partitioning = "run_id"),

    chromatographic_peak_feature = .mzs_table(c(
        chromatographic_peak_id_ = "i64!",
        feature_id_ = "i64!",
        peak_index_in_run = "i32!",
        is_representative = "bool!"),
        sorted_by = c("feature_id_", "chromatographic_peak_id_"),
        key = NA_character_),

    abundance = .mzs_table(c(
        feature_id_ = "i64!",
        chromatographic_peak_id_ = "i64",
        assay_id_ = "i32!",
        quantity_kind = "dict!",
        value = "f64",
        value_source = "dict"),
        sorted_by = c("feature_id_", "assay_id_", "quantity_kind"),
        key = NA_character_),

    assay = .mzs_table(c(
        assay_id_ = "i32!",
        assay_name = "str!",
        sample_id_ = "i32",
        source = "dict",
        run_id = "dict"),
        sorted_by = "assay_id_"),

    sample = .mzs_table(c(
        sample_id_ = "i32!",
        sample_name = "str",
        sample_class = "dict"),
        sorted_by = "sample_id_"),

    spectrum_feature = .mzs_table(c(
        association_id_ = "i64!",
        .mzs_ref_columns("ms2", usi = TRUE),
        feature_id_ = "i64!",
        chromatographic_peak_id_ = "i64",
        link_mode = "dict!",
        run_scoped = "bool!",
        precursor_mz_used = "f64",
        precursor_mz_error = "f32",
        precursor_retention_time = "f64",
        precursor_purity = "f32"),
        sorted_by = c("feature_id_", "ms2_run_id", "ms2_spectrum_id_",
                      "ms2_native_id", "association_id_"),
        partitioning = "ms2_run_id"),

    merge_member = .mzs_table(c(
        .mzs_ref_columns("merged"),
        .mzs_ref_columns("member"),
        member_rank = "i32!",
        members_complete = "bool!"),
        sorted_by = c("merged_spectrum_id_", "member_rank"),
        key = NA_character_),

    compound = .mzs_table(c(
        compound_id_ = "i64!",
        inchikey = "str",
        inchikey_block1 = "str",
        inchi = "str",
        smiles = "str",
        chemical_formula = "str",
        chemical_name = "str",
        theoretical_neutral_mass = "f64",
        average_molecular_weight = "f64",
        compound_class = "dict",
        uri = "str"),
        sorted_by = c("inchikey", "compound_id_")),

    compound_xref = .mzs_table(c(
        compound_id_ = "i64!",
        database = "dict!",
        identifier = "str!",
        database_version = "dict"),
        sorted_by = c("compound_id_", "database", "identifier"),
        key = NA_character_),

    compound_synonym = .mzs_table(c(
        compound_id_ = "i64!",
        synonym = "str!"),
        sorted_by = c("compound_id_", "synonym"),
        key = NA_character_),

    evidence = .mzs_table(c(
        evidence_id_ = "i64!",
        evidence_input_id = "i64!",
        .mzs_ref_columns("query"),
        .mzs_ref_columns("reference"),
        feature_id_ = "i64",
        compound_id_ = "i64",
        compound_match_level = "dict!",
        identification_method = "dict!",
        ms_level = "i32",
        adduct_ion = "dict",
        exp_mass_to_charge = "f64",
        theoretical_mass_to_charge = "f64",
        charge = "i32",
        rank = "i32",
        database_identifier = "str",
        chemical_name = "str",
        smiles = "str",
        inchi = "str",
        chemical_formula = "str",
        uri = "str",
        software_id_ = "i32!"),
        sorted_by = c("query_spectrum_id_", "rank", "evidence_id_"),
        partitioning = "query_run_id"),

    evidence_score = .mzs_table(c(
        evidence_id_ = "i64!",
        score_term = "dict!",
        score_value = "f64",
        score_kind = "dict!",
        rescaling_scope = "dict"),
        sorted_by = c("evidence_id_", "score_term"),
        key = NA_character_),

    software = .mzs_table(c(
        software_id_ = "i32!",
        name = "str!",
        version = "str",
        parameters = "str"),
        sorted_by = "software_id_"),

    coverage = .mzs_table(c(
        .mzs_ref_columns("query"),
        candidates_considered = "i64",
        hits_retained = "i32!",
        truncated = "bool!",
        boundary_score = "f64"),
        sorted_by = c("query_source", "query_run_id", "query_spectrum_id_",
                      "query_native_id", "activity"),
        key = NA_character_),

    conversion = .mzs_table(c(
        conversion_id_ = "i32!",
        converter_name = "str!",
        converter_version = "str!",
        route = "dict!",
        mzstack_version = "str!",
        results_version = "str!",
        source_format = "dict!",
        source_format_version = "str",
        source_uri = "str",
        source_checksum = "str",
        source_checksum_algorithm = "dict",
        converted_at = "str!",
        coverage_manifest = "str!",
        loss_ledger = "str!"),
        sorted_by = "conversion_id_"),

    tool_provenance = .mzs_table(c(
        tool_provenance_id_ = "i64!",
        tool_name = "str!",
        tool_version = "str",
        step_name = "str!",
        parameters = "str",
        parameters_completeness = "dict!"),
        sorted_by = "tool_provenance_id_"),

    source_identifier = .mzs_table(c(
        source_identifier_id_ = "i64!",
        target_table = "dict!",
        target_key = "i64!",
        source_table = "dict",
        source_column = "dict",
        source_value = "str!",
        source_value_type = "dict!"),
        sorted_by = c("target_table", "target_key", "source_identifier_id_")),

    loss_ledger = .mzs_table(c(
        loss_ledger_id_ = "i64!",
        source_construct = "str!",
        disposition = "dict!",
        affected_rows = "i64",
        target_ref = "str",
        reason = "str"),
        sorted_by = "loss_ledger_id_"),

    ## msPurity extensions.

    ## Every assessed MS2 scan and its precursor purity, linked to a feature
    ## or not. Raw scans are referenced, not copied.
    x_mspurity_scan = .mzs_table(c(
        scan_annotation_id_ = "i64!",
        .mzs_ref_columns("ms2", usi = TRUE),
        precursor_mz = "f64",
        precursor_retention_time = "f64",
        precursor_intensity = "f64",
        retention_time = "f64",
        scan_number = "i32",
        precursor_scan_number = "i32",
        precursor_nearest = "i32",
        a_mz = "f64",
        a_purity = "f64",
        a_peak_count = "f64",
        i_mz = "f64",
        i_purity = "f64",
        i_peak_count = "f64",
        in_purity = "f64",
        in_peak_count = "f64",
        purity_pass_flag = "bool"),
        sorted_by = c("ms2_run_id", "ms2_spectrum_id_", "ms2_native_id",
                      "scan_annotation_id_"),
        partitioning = "ms2_run_id", extension = TRUE),

    ## filterFragSpectra() peak flags; a peak is identified by its m/z and
    ## intensity within the scan.
    x_mspurity_scan_peak = .mzs_table(c(
        scan_annotation_id_ = "i64!",
        peak_rank = "i32!",
        mz = "f64!",
        intensity = "f64!",
        snr = "f64",
        relative_intensity = "f64",
        purity_pass_flag = "bool",
        intensity_pass_flag = "bool",
        ra_pass_flag = "bool",
        snr_pass_flag = "bool"),
        sorted_by = c("scan_annotation_id_", "peak_rank"),
        key = NA_character_, extension = TRUE),

    ## Feature-level annotations (MetFrag, SIRIUS CSI:FingerID, ProbMetab,
    ## MS1 lookup), one row per tool, feature, compound and score.
    x_mspurity_feature_annotation = .mzs_table(c(
        feature_annotation_id_ = "i64!",
        feature_id_ = "i64!",
        compound_id_ = "i64",
        tool = "dict!",
        score_term = "dict!",
        score_value = "f64",
        score_kind = "dict!",
        rescaling_scope = "dict",
        tool_rank = "i32"),
        sorted_by = c("feature_id_", "tool", "feature_annotation_id_"),
        extension = TRUE),

    ## combineAnnotations() weighted scores per feature and compound; tied
    ## scores share a rank.
    x_mspurity_combined_annotation = .mzs_table(c(
        combined_annotation_id_ = "i64!",
        feature_id_ = "i64!",
        compound_id_ = "i64!",
        weighted_score = "f64",
        rank = "i32",
        sm_score = "f64",
        metfrag_score = "f64",
        sirius_csifingerid_score = "f64",
        probmetab_score = "f64",
        ms1_lookup_score = "f64",
        biosim_score = "f64"),
        sorted_by = c("feature_id_", "rank", "combined_annotation_id_"),
        extension = TRUE)
)

#' @noRd
.mzs_table_schema <- function(name) {
    s <- .MZS_TABLES[[name]]
    if (is.null(s))
        stop("No results table '", name, "' is defined.", call. = FALSE)
    s
}

#' The reference prefixes of a table, e.g. `c("query", "reference")`.
#'
#' @noRd
.mzs_ref_prefixes <- function(name) {
    cols <- .mzs_table_schema(name)$columns$name
    p <- sub("_source$", "", grep("_source$", cols, value = TRUE))
    p[paste0(p, "_spectrum_id_") %in% cols]
}

# Allowed values of discriminator columns.
.MZS_LINK_MODES <- c("chromatographic_peak", "feature_width",
                     "isolation_window", "manual")
.MZS_MATCH_LEVELS <- c("exact_structure", "skeleton", "formula", "mass_only")
.MZS_SCORE_KINDS <- c("raw", "rescaled", "constant", "count")
.MZS_DISPOSITIONS <- c("not_present_in_source", "no_such_concept_in_source",
                       "no_such_concept_in_target", "preserved_opaque",
                       "dropped_by_configuration", "dropped_unmappable",
                       "refused")
.MZS_COVERAGE_STATES <- c("mapped", "mapped_partial",
                          "not_present_in_source", "dropped", "refused")
.MZS_COMPLETENESS <- c("complete", "partial", "absent")

#' @noRd
.mzs_arrow_type <- function(code) {
    switch(code,
           i64 = arrow::int64(), i32 = arrow::int32(),
           f64 = arrow::float64(), f32 = arrow::float32(),
           bool = arrow::boolean(), str = arrow::utf8(),
           dict = arrow::dictionary(index_type = arrow::int32(),
                                    value_type = arrow::utf8()),
           stop("Unknown column type code '", code, "'.", call. = FALSE))
}

#' Coerce an R vector to the R type of a type code.
#'
#' @noRd
.mzs_as_r_type <- function(x, code) {
    if (is.factor(x))
        x <- as.character(x)
    switch(code,
           i64 = , f64 = , f32 = as.double(x),
           i32 = as.integer(x),
           bool = as.logical(x),
           str = , dict = as.character(x),
           stop("Unknown column type code '", code, "'.", call. = FALSE))
}

#' Column names breaking the naming rules: a trailing underscore is only for
#' integer `_id_` keys, and only run ids, native ids and names inherited
#' from mzTab-M may end in `_id`.
#'
#' @noRd
.mzs_bad_names <- function(names, types) {
    integer_key <- types %in% c("i32", "i64")
    trailing <- endsWith(names, "_") & !(endsWith(names, "_id_") &
                                         integer_key)
    ends_id <- grepl("_id$", names) &
        !(names == "run_id" | grepl("_(run|native)_id$", names) |
          names %in% c("evidence_input_id"))
    names[trailing | ends_id]
}
