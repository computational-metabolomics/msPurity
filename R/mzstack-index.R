# The results index: per results table and per run with peak-annotation
# columns, the entity, data kind, sort order and CV bindings of columns.
# It lives under `results/`, not `index/`, because it is authored, not
# derived.
#
# PSI-MS and UO accessions are checked against inst/mzstack/cv/terms.tsv.
# Quantities with no PSI-MS term use the local MSPURITY namespace
# (inst/mzstack/cv/MSPURITY.obo).

.MZS_CV <- list(
    MS = list(id = "MS", full_name = "PSI-MS",
              uri = "http://purl.obolibrary.org/obo/ms/psi-ms.obo",
              version = "4.1.249"),
    UO = list(id = "UO", full_name = "Units of Measurement Ontology",
              uri = "http://purl.obolibrary.org/obo/uo.obo",
              version = "releases/2023-05-25"),
    MSPURITY = list(
        id = "MSPURITY", full_name = "msPurity local terms",
        uri = paste0("https://raw.githubusercontent.com/computational-",
                     "metabolomics/msPurity/master/inst/mzstack/cv/",
                     "MSPURITY.obo"),
        version = "1"))

# Column -> c(accession, unit), wherever the column appears.
.MZS_BINDINGS <- list(
    exp_mass_to_charge = c("MS:1000744", NA),
    precursor_mz_used = c("MS:1000744", NA),
    precursor_mz = c("MS:1000744", NA),
    retention_time_in_seconds = c("MS:1000894", "UO:0000010"),
    retention_time_in_seconds_start = c("MS:1000894", "UO:0000010"),
    retention_time_in_seconds_end = c("MS:1000894", "UO:0000010"),
    retention_time_in_seconds_start_adjusted = c("MS:1000894", "UO:0000010"),
    retention_time_in_seconds_end_adjusted = c("MS:1000894", "UO:0000010"),
    precursor_retention_time = c("MS:1000894", "UO:0000010"),
    retention_time = c("MS:1000894", "UO:0000010"),
    charge = c("MS:1000041", NA),
    ms_level = c("MS:1000511", NA),
    precursor_purity = c("MSPURITY:0000001", NA),
    in_purity = c("MSPURITY:0000001", NA),
    a_purity = c("MSPURITY:0000002", NA),
    a_mz = c("MSPURITY:0000003", NA),
    a_peak_count = c("MSPURITY:0000004", NA),
    i_purity = c("MSPURITY:0000005", NA),
    i_mz = c("MSPURITY:0000006", NA),
    i_peak_count = c("MSPURITY:0000007", NA),
    in_peak_count = c("MSPURITY:0000008", NA),
    precursor_nearest = c("MSPURITY:0000009", NA),
    sn = c("MSPURITY:0000016", NA),
    contributor_count = c("MSPURITY:0000017", NA),
    contributor_total = c("MSPURITY:0000018", NA),
    intensity_rsd = c("MSPURITY:0000019", NA),
    cluster = c("MSPURITY:0000020", NA),
    x_mspurity_ra = c("MSPURITY:0000021", "UO:0000187"),
    x_mspurity_in_purity = c("MSPURITY:0000022", NA),
    x_mspurity_snr_pass_flag = c("MSPURITY:0000023", NA),
    x_mspurity_minnum_pass_flag = c("MSPURITY:0000024", NA),
    x_mspurity_minfrac_pass_flag = c("MSPURITY:0000025", NA),
    x_mspurity_ra_pass_flag = c("MSPURITY:0000026", NA),
    weighted_score = c("MSPURITY:0000027", NA),
    metfrag_score = c("MSPURITY:0000028", NA),
    sirius_csifingerid_score = c("MSPURITY:0000029", NA),
    probmetab_score = c("MSPURITY:0000030", NA),
    ms1_lookup_score = c("MSPURITY:0000031", NA),
    biosim_score = c("MSPURITY:0000032", NA),
    sm_score = c("MSPURITY:0000033", NA))

# Score terms written as `evidence_score.score_term` values.
.MZS_SCORE_TERMS <- c(dpc = "MSPURITY:0000010", rdpc = "MSPURITY:0000011",
                      cdpc = "MSPURITY:0000012", mcount = "MSPURITY:0000013",
                      allcount = "MSPURITY:0000014",
                      mpercent = "MSPURITY:0000015")

# Identification method of every evidence row spectralMatching() writes.
.MZS_LIBRARY_SEARCH <- "MS:1001031"

#' @noRd
.mzs_index_bindings <- function(columns) {
    out <- list()
    for (col in columns) {
        b <- .MZS_BINDINGS[[col]]
        if (is.null(b))
            next
        e <- list(path = col, accession = b[[1L]])
        if (!is.na(b[[2L]]))
            e$unit <- b[[2L]]
        out[[length(out) + 1L]] <- e
    }
    out
}

#' @noRd
.mzs_index_table <- function(name, columns) {
    list(name = name, entity_type = name, data_kind = "table",
         sorted_by = .mzs_array(.mzs_table_schema(name)$sorted_by),
         column_mapping = .mzs_index_bindings(columns))
}

#' The index entry declaring a run's peak-annotation columns.
#'
#' @noRd
.mzs_index_peaks <- function(run_id, columns) {
    list(name = run_id, entity_type = "spectrum",
         data_kind = "peak_annotation",
         column_mapping = .mzs_index_bindings(columns))
}

#' Merge entries into an index and recompute its ontology list.
#'
#' @noRd
.mzs_index_merge <- function(old, entries) {
    files <- old$files %||% list()
    have <- vapply(files, function(f) paste(f$data_kind, f$name),
                   character(1))
    for (e in entries) {
        i <- match(paste(e$data_kind, e$name), have)
        if (is.na(i)) {
            files[[length(files) + 1L]] <- e
            have <- c(have, paste(e$data_kind, e$name))
        } else {
            files[[i]] <- e
        }
    }
    used <- unique(unlist(lapply(files, function(f)
        lapply(f$column_mapping, function(b)
            sub(":.*$", "", c(b$accession, b$unit))))))
    cv <- unname(.MZS_CV[intersect(names(.MZS_CV), used)])
    if (!length(cv))
        cv <- list(.MZS_CV$MS)
    list(files = files, cv_list = cv)
}

#' Read a dataset's results index, or `NULL`.
#'
#' @noRd
.mzs_index_read <- function(path, m) {
    rel <- m$results$column_mapping_path
    if (is.null(rel))
        return(NULL)
    fl <- file.path(path, rel)
    if (!file.exists(fl))
        return(NULL)
    jsonlite::fromJSON(fl, simplifyVector = FALSE)
}

#' Accessions known to the vendored term lists.
#'
#' @noRd
.mzs_known_terms <- function() {
    fl <- system.file("mzstack", "cv", "terms.tsv", package = "msPurity")
    utils::read.delim(fl, stringsAsFactors = FALSE)$accession
}
