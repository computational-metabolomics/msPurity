# combineAnnotations() on an mzStack results dataset.
#
# The SQLite scoring is reused unchanged: the tables it reads are rebuilt in a
# temporary database, so both formats rank compounds identically.
# Annotations are per feature, not per spectrum, so they go in extension
# tables: x_mspurity_feature_annotation (one row per tool, feature, compound
# and score) and x_mspurity_combined_annotation (weighted score and rank).

# Per-tool score columns, terms and kinds.
.MZS_CA_TOOLS <- list(
    spectral_matching = c(col = "sm_score", term = "MSPURITY:0000033",
                          kind = "raw"),
    metfrag = c(col = "metfrag_score", term = "MSPURITY:0000028",
                kind = "raw"),
    sirius_csifingerid = c(col = "sirius_score", term = "MSPURITY:0000029",
                           kind = "rescaled"),
    probmetab = c(col = "probmetab_score", term = "MSPURITY:0000030",
                  kind = "raw"),
    ms1_lookup = c(col = "ms1_lookup_score", term = "MSPURITY:0000031",
                   kind = "constant"),
    biosim = c(col = "biosim_max_score", term = "MSPURITY:0000032",
               kind = "raw"))

#' Rebuild the tables combineAnnotations() reads in a temporary SQLite
#' database.
#'
#' @noRd
.mzs_ca_database <- function(path, m, db) {
    can <- .mzs_canonical_mzstack(path)
    tab <- function(n) .mzs_read_table(path, m, n)
    ev <- tab("evidence")
    es <- tab("evidence_score")
    cmp <- if (!is.null(m$results$tables$compound)) tab("compound")
           else .mzs_empty_table("compound")
    si <- tab("source_identifier")
    grp_of <- function(fid) as.numeric(si$source_value[match(
        paste("feature", "grpid", fid),
        paste(si$target_table, si$source_column, si$target_key))])
    ev <- ev[!is.na(ev$feature_id_), , drop = FALSE]
    ref <- paste(ev$reference_source, ev$reference_run_id,
                 ev$reference_spectrum_id_, ev$reference_native_id)
    lpid <- match(ref, unique(ref))
    dpc <- es$score_value[es$score_term == .MZS_SCORE_TERMS[["dpc"]]][match(
        ev$evidence_id_, es$evidence_id_[es$score_term ==
                                         .MZS_SCORE_TERMS[["dpc"]]])]
    ik <- cmp$inchikey[match(ev$compound_id_, cmp$compound_id_)]
    xcms_match <- data.frame(grpid = grp_of(ev$feature_id_), lpid = lpid,
                             mid = ev$evidence_id_, dpc = dpc,
                             library_precursor_type = ev$adduct_ion,
                             stringsAsFactors = FALSE)
    l_s_peak_meta <- unique(data.frame(id = lpid, inchikey_id = ik,
                                       accession = ev$reference_native_id,
                                       stringsAsFactors = FALSE))
    keys <- unique(cmp$inchikey[!is.na(cmp$inchikey)])
    metab_compound <- data.frame(
        inchikey_id = keys,
        name = cmp$chemical_name[match(keys, cmp$inchikey)],
        pubchem_id = NA_character_, chemspider_id = NA_character_,
        other_names = NA_character_,
        exact_mass = cmp$theoretical_neutral_mass[match(keys, cmp$inchikey)],
        molecular_formula = cmp$chemical_formula[match(keys, cmp$inchikey)],
        molecular_weight = NA_real_, compound_class = NA_character_,
        smiles = cmp$smiles[match(keys, cmp$inchikey)],
        created_at = NA_character_, updated_at = NA_character_,
        stringsAsFactors = FALSE)
    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    on.exit(DBI::dbDisconnect(con))
    w <- function(name, df) DBI::dbWriteTable(con, name, df,
                                              row.names = FALSE)
    w("c_peak_groups", can$c_peak_groups)
    w("c_peaks", can$c_peaks)
    w("c_peak_X_c_peak_group", can$c_peak_X_c_peak_group)
    if ("cid" %in% names(can$scan_link))
        w("c_peak_X_s_peak_meta", can$scan_link)
    w("s_peak_meta", can$scans)
    w("xcms_match", xcms_match)
    w("l_s_peak_meta", l_s_peak_meta)
    w("metab_compound", metab_compound)
    list(compound = cmp, lpid = lpid, evidence = ev)
}

#' combineAnnotations(format = "mzstack").
#'
#' @noRd
.combineAnnotations_mzstack <- function(a) {
    .mzs_require("combineAnnotations(format = \"mzstack\")")
    started <- .mzs_now()
    if (isTRUE(a$ms1_lookup_checkAdducts))
        .mzs_abort("unsupported", "ms1_lookup_checkAdducts = TRUE needs ",
                   "CAMERA annotation, which the mzStack route does not ",
                   "map.", data = list(constructs = "ms1_lookup_checkAdducts"))
    if (!.mzs_is_dataset(a$sm_resultPth))
        .mzs_abort("format", "'", a$sm_resultPth, "' is not an mzStack ",
                   "dataset.")
    path <- normalizePath(a$sm_resultPth)
    m <- .mzs_read_manifest(path)
    if (!identical(m$role, "results") || is.null(m$results$tables$evidence))
        .mzs_abort("capability", "'", path, "' holds no spectral-matching ",
                   "evidence to combine; run spectralMatching(format = ",
                   "\"mzstack\", updateDb = TRUE) first.")

    tmp <- tempfile("mzs-combine-")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    db <- file.path(tmp, "combine.sqlite")
    built <- .mzs_ca_database(path, m, db)
    core <- a[names(formals(.combineAnnotations_sqlite))]
    core$sm_resultPth <- db
    core$outPth <- NA
    summary <- do.call(.combineAnnotations_sqlite, core)

    con <- DBI::dbConnect(RSQLite::SQLite(), db)
    ca <- suppressWarnings(DBI::dbGetQuery(con,
                                           "SELECT * FROM combined_annotations"))
    mc <- suppressWarnings(DBI::dbGetQuery(con, "SELECT * FROM metab_compound"))
    DBI::dbDisconnect(con)

    target <- path
    forked <- NULL
    if (!is.na(a$outPth)) {
        forked <- .mzs_fork(path, a$outPth, overwrite = FALSE)
        target <- normalizePath(a$outPth)
    }
    tm <- .mzs_read_manifest(target)
    si <- .mzs_read_table(target, tm, "source_identifier")
    fid_of <- function(grpid) {
        s <- si[si$target_table == "feature" & si$source_column == "grpid", ]
        s$target_key[match(as.numeric(grpid), as.numeric(s$source_value))]
    }

    ## New compound revision: existing ids kept and enriched, plus new
    ## compounds from the tools.
    old <- built$compound
    mc$inchikey <- .mzs_inchikey(mc$inchikey %||% mc$inchikey_id)
    k <- match(old$inchikey, mc$inchikey)
    fill <- function(v, w) ifelse(is.na(v), w, v)
    if (nrow(old)) {
        old$chemical_name <- fill(old$chemical_name, mc$name[k])
        old$chemical_formula <- fill(old$chemical_formula,
                                     mc$molecular_formula[k])
        old$theoretical_neutral_mass <- fill(
            old$theoretical_neutral_mass,
            suppressWarnings(as.numeric(mc$exact_mass[k])))
        old$smiles <- fill(old$smiles, mc$smiles_canonical[k] %||% NA)
        old$inchi <- fill(old$inchi, mc$inchi[k] %||% NA)
    }
    add <- mc[!is.na(mc$inchikey) & !mc$inchikey %in% old$inchikey, ,
              drop = FALSE]
    add <- add[!duplicated(add$inchikey), , drop = FALSE]
    new <- if (nrow(add)) data.frame(
        compound_id_ = max(c(0, old$compound_id_)) + seq_len(nrow(add)),
        inchikey = add$inchikey, inchikey_block1 = substr(add$inchikey, 1, 14),
        inchi = add$inchi %||% NA_character_,
        smiles = add$smiles_canonical %||% NA_character_,
        chemical_formula = add$molecular_formula %||% NA_character_,
        chemical_name = add$name %||% NA_character_,
        theoretical_neutral_mass = suppressWarnings(as.numeric(
            add$exact_mass %||% NA)),
        average_molecular_weight = NA_real_, compound_class = NA_character_,
        uri = NA_character_, stringsAsFactors = FALSE)
    keep <- intersect(.mzs_table_schema("compound")$columns$name, names(old))
    compound <- .mzs_rbind_fill(old[, setdiff(keep, "activity"), drop = FALSE],
                                new)
    compound <- compound[order(compound$compound_id_), , drop = FALSE]
    cid_of <- function(ik) compound$compound_id_[match(.mzs_inchikey(ik),
                                                       compound$inchikey)]

    ## Cross-references from the compound database.
    xdb <- c(pubchem_cids = "PubChem Compound", kegg_cids = "KEGG Compound",
             hmdb_ids = "HMDB")
    have_x <- if (!is.null(tm$results$tables$compound_xref))
        .mzs_read_table(target, tm, "compound_xref")
    xref <- do.call(.mzs_rbind_fill, lapply(names(xdb), function(c) {
        if (is.null(mc[[c]]))
            return(NULL)
        v <- .mzs_na_empty(mc[[c]])
        ok <- !is.na(v) & !is.na(mc$inchikey)
        if (!any(ok))
            return(NULL)
        ids <- strsplit(v[ok], "[,;][[:space:]]*")
        data.frame(compound_id_ = rep(cid_of(mc$inchikey[ok]), lengths(ids)),
                   database = xdb[[c]], identifier = unlist(ids),
                   stringsAsFactors = FALSE)
    }))
    if (!is.null(xref)) {
        xref <- unique(xref[!is.na(xref$compound_id_) & nzchar(xref$identifier),
                            , drop = FALSE])
        if (!is.null(have_x))
            xref <- xref[!paste(xref$compound_id_, xref$database,
                                xref$identifier) %in%
                         paste(have_x$compound_id_, have_x$database,
                               have_x$identifier), , drop = FALSE]
    }

    ## Combined scores and per-tool contributions.
    num <- function(c) if (c %in% names(ca)) suppressWarnings(
        as.numeric(ca[[c]])) else rep(NA_real_, nrow(ca))
    comb <- data.frame(
        combined_annotation_id_ = seq_len(nrow(ca)),
        feature_id_ = fid_of(ca$grpid), compound_id_ = cid_of(ca$inchikey),
        weighted_score = num("wscore"), rank = as.integer(num("rank")),
        sm_score = num("sm_score"), metfrag_score = num("metfrag_score"),
        sirius_csifingerid_score = num("sirius_score"),
        probmetab_score = num("probmetab_score"),
        ms1_lookup_score = num("ms1_lookup_score"),
        biosim_score = num("biosim_max_score"),
        x_mspurity_adduct_overall = .mzs_na_empty(ca$adduct_overall),
        stringsAsFactors = FALSE)
    if (anyNA(comb$feature_id_) || anyNA(comb$compound_id_))
        .mzs_abort("semantic", "combineAnnotations() returned features or ",
                   "compounds absent from the dataset.")
    fa <- do.call(rbind, lapply(names(.MZS_CA_TOOLS), function(t) {
        spec <- .MZS_CA_TOOLS[[t]]
        v <- num(spec[["col"]])
        ok <- !is.na(v) & v != 0
        if (!any(ok))
            return(NULL)
        data.frame(feature_id_ = comb$feature_id_[ok],
                   compound_id_ = comb$compound_id_[ok], tool = t,
                   score_term = spec[["term"]], score_value = v[ok],
                   score_kind = spec[["kind"]],
                   rescaling_scope = if (spec[["kind"]] == "rescaled")
                       "feature" else NA_character_,
                   stringsAsFactors = FALSE)
    }))
    if (!is.null(fa))
        fa$feature_annotation_id_ <- seq_len(nrow(fa))

    files <- c(metfrag = a$metfrag_resultPth,
               sirius_csifingerid = a$sirius_csi_resultPth,
               probmetab = a$probmetab_resultPth,
               ms1_lookup = a$ms1_lookup_resultPth,
               compound_database = if (identical(a$compoundDbType, "sqlite"))
                   a$compoundDbPth else NA)
    files <- files[!is.na(files) & file.exists(files)]
    append <- list(x_mspurity_combined_annotation = comb)
    if (!is.null(fa))
        append$x_mspurity_feature_annotation <- fa
    if (!is.null(xref) && nrow(xref))
        append$compound_xref <- xref
    if (!is.null(tm$results$tables$loss_ledger))
        append$loss_ledger <- .mzs_loss_ledger(list(
            list(construct = "combineAnnotations: per-tool hits",
                 disposition = "dropped_by_configuration",
                 target_ref = "x_mspurity_feature_annotation",
                 reason = paste("As combineAnnotations() combines them, only",
                                "each tool's best hit per feature and",
                                "compound is kept; the tool result files are",
                                "recorded as inputs with their digests.")),
            list(construct = "combineAnnotations: kegg, hmdb, pubchem tables",
                 disposition = "no_such_concept_in_target",
                 target_ref = "compound_xref",
                 reason = paste("Compound database lookups; the identifiers",
                                "found are in compound_xref."))))
    params <- a[c("ms1_lookup_dbSource", "ms1_lookup_checkAdducts",
                  "ms1_lookup_keepAdducts", "weights", "compoundDbType",
                  "summaryOutput")]
    .mzs_commit_update(
        target, append = append, supersede = list(compound = compound),
        activity = list(
            action = "annotate", started = started,
            fn = "combineAnnotations",
            inputs = lapply(files, .mzs_file_input),
            parameters = c(params, list(
                forked_from = forked,
                tool_files = as.list(stats::setNames(basename(files),
                                                     names(files)))))))
    summary
}
