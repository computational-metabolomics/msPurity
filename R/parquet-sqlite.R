# Reading msPurity's SQLite output for conversion to Parquet.
#
# SQLite column types depend on the first row (flags may be TEXT '0'/'1', ids
# REAL or TEXT), so every mapped column is cast explicitly. Unknown tables and
# columns are refused, except peak-picker/xcms-dependent columns, which pass
# through.

# Tables the converter maps.
.MZS_SQLITE_REQUIRED <- c("fileinfo", "c_peaks", "c_peak_groups",
                          "c_peak_X_c_peak_group", "s_peak_meta", "s_peaks")
.MZS_SQLITE_OPTIONAL <- c("source", "metab_compound", "sm_matches",
                          "l_s_peak_meta", "xcms_match",
                          "c_peak_X_s_peak_meta",
                          "c_peak_group_X_s_peak_meta")

# Refused tables: CAMERA annotation and combineAnnotations() outputs (use
# combineAnnotations(format = "parquet") instead).
.MZS_SQLITE_REFUSED <- c("adduct_rules", "neutral_masses",
                         "adduct_annotations", "isotope_annotations",
                         "combined_annotations", "metfrag_results",
                         "sirius_csifingerid_results", "probmetab_results",
                         "ms1_lookup_results", "kegg", "hmdb", "pubchem",
                         "eics")

# Allowed columns of fixed-schema tables.
.MZS_SQLITE_COLUMNS <- list(
    fileinfo = c("fileid", "filename", "filepth", "nm_save", "class"),
    s_peak_meta = c(
        "pid", "seqNum", "acquisitionNum", "precursorIntensity",
        "precursorMZ", "precursorRT", "precursorScanNum", "retentionTime",
        "precursorNearest", "aMz", "aPurity", "apkNm", "iMz", "iPurity",
        "ipkNm", "inPkNm", "inPurity", "purity_pass_flag", "name",
        "collision_energy", "ms_level", "accession", "resolution",
        "polarity", "fragmentation_type", "precursor_type",
        "instrument_type", "instrument", "copyright", "column",
        "mass_accuracy", "mass_error", "origin", "splash",
        "retention_index", "retention_time", "inchikey_id", "sourceid",
        "precursor_mz", "spectrum_type", "grpid", "fileid"),
    s_peaks = c(
        "sid", "fileid", "mz", "i", "snr", "ra", "type", "scan",
        "purity_pass_flag", "intensity_pass_flag", "ra_pass_flag",
        "snr_pass_flag", "pass_flag", "cl", "rsd", "count", "total",
        "inPurity", "frac", "minnum_pass_flag", "minfrac_pass_flag",
        "grpid", "pid"),
    c_peak_X_c_peak_group = c("cXg_id", "idi", "bestpeak", "grpid", "cid"),
    c_peak_X_s_peak_meta = c("cXp_id", "pid", "cid"),
    c_peak_group_X_s_peak_meta = c("gXp_id", "pid", "grpid"),
    source = c("id", "name", "parsing_software"),
    metab_compound = c(
        "inchikey_id", "name", "pubchem_id", "chemspider_id", "other_names",
        "exact_mass", "molecular_formula", "molecular_weight",
        "compound_class", "smiles", "created_at", "updated_at"),
    sm_matches = c(
        "mid", "lpid", "qpid", "dpc", "rdpc", "cdpc", "mcount", "allcount",
        "mpercent", "library_rt", "query_rt", "rtdiff",
        "library_precursor_mz", "query_precursor_mz",
        "library_precursor_ion_purity", "query_precursor_ion_purity",
        "library_accession", "library_precursor_type",
        "library_entry_name", "inchikey", "library_source_name",
        "library_compound_name"))

#' Read the tables of an msPurity database, typed.
#'
#' @return list of data.frames by table name, plus `lexical` (source text
#'     form of identifiers) and `tables` (tables present).
#'
#' @noRd
.mzs_read_sqlite <- function(db) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
    on.exit(DBI::dbDisconnect(con))
    have <- DBI::dbListTables(con)
    miss <- setdiff(.MZS_SQLITE_REQUIRED, have)
    if (length(miss))
        .mzs_abort("format", "'", db, "' is not a database written by ",
                   "msPurity's createDatabase(): it has no table(s) ",
                   paste(miss, collapse = ", "), ".")
    if (!any(c("c_peak_X_s_peak_meta", "c_peak_group_X_s_peak_meta") %in%
             have))
        .mzs_abort("format", "'", db, "' links no MS2 spectrum to a ",
                   "feature: it has neither c_peak_X_s_peak_meta nor ",
                   "c_peak_group_X_s_peak_meta.")
    refused <- intersect(have, .MZS_SQLITE_REFUSED)
    unknown <- setdiff(have, c(.MZS_SQLITE_REQUIRED, .MZS_SQLITE_OPTIONAL,
                               .MZS_SQLITE_REFUSED))
    if (length(refused) || length(unknown))
        .mzs_abort("unsupported", "'", db, "' holds table(s) this converter ",
                   "does not map: ", paste(c(refused, unknown),
                                          collapse = ", "),
                   ". Refusing to convert the rest silently (R-081).",
                   data = list(constructs = c(refused, unknown)))
    rd <- function(t) if (t %in% have) suppressWarnings(DBI::dbGetQuery(
        con, sprintf('SELECT * FROM "%s" ORDER BY rowid', t)))
    tabs <- intersect(c(.MZS_SQLITE_REQUIRED, .MZS_SQLITE_OPTIONAL), have)
    src <- lapply(stats::setNames(tabs, tabs), rd)
    for (t in intersect(names(.MZS_SQLITE_COLUMNS), tabs)) {
        bad <- setdiff(names(src[[t]]), .MZS_SQLITE_COLUMNS[[t]])
        if (length(bad))
            .mzs_abort("unsupported", "Table ", t, " of '", db, "' has ",
                       "column(s) this converter does not map: ",
                       paste(bad, collapse = ", "), " (R-081).",
                       data = list(constructs = paste0(t, ".", bad)))
    }

    lex <- function(t, c) if (t %in% have && c %in% names(src[[t]]))
        DBI::dbGetQuery(con, sprintf(
            'SELECT CAST("%s" AS TEXT) AS v FROM "%s" ORDER BY rowid',
            c, t))$v
    src$lexical <- list(
        fileid = lex("fileinfo", "fileid"),
        grpid = lex("c_peak_groups", "grpid"),
        cid = lex("c_peaks", "cid"),
        pid = lex("s_peak_meta", "pid"),
        mid = lex("sm_matches", "mid"),
        inchikey_id = lex("metab_compound", "inchikey_id"))
    src$tables <- have

    int <- function(x) as.integer(round(suppressWarnings(as.numeric(x))))
    num <- function(x) suppressWarnings(as.numeric(x))
    flag <- function(x) {
        v <- as.character(x)
        out <- rep(NA, length(v))
        out[v %in% c("1", "1.0", "TRUE", "True", "true")] <- TRUE
        out[v %in% c("0", "0.0", "FALSE", "False", "false")] <- FALSE
        out
    }
    cast <- function(df, ints = character(), nums = character(),
                     flags = character(), chrs = character()) {
        for (c in intersect(ints, names(df))) df[[c]] <- int(df[[c]])
        for (c in intersect(nums, names(df))) df[[c]] <- num(df[[c]])
        for (c in intersect(flags, names(df))) df[[c]] <- flag(df[[c]])
        for (c in intersect(chrs, names(df)))
            df[[c]] <- .mzs_na_empty(df[[c]])
        df
    }
    src$fileinfo <- cast(src$fileinfo, ints = "fileid",
                         chrs = c("filename", "filepth", "nm_save", "class"))
    src$c_peaks <- cast(src$c_peaks, ints = c("cid", "fileid"),
                        nums = setdiff(names(src$c_peaks),
                                       c("cid", "fileid")))
    src$c_peak_groups <- cast(src$c_peak_groups, ints = "grpid",
                              nums = c("mz", "mzmin", "mzmax", "rt",
                                       "rtmin", "rtmax", "npeaks"),
                              chrs = c("grp_name", "peakidx"))
    src$c_peak_X_c_peak_group <- cast(src$c_peak_X_c_peak_group,
                                      ints = c("cXg_id", "idi", "bestpeak",
                                               "grpid", "cid"))
    if (!is.null(src$c_peak_X_s_peak_meta))
        src$c_peak_X_s_peak_meta <- cast(src$c_peak_X_s_peak_meta,
                                         ints = c("cXp_id", "pid", "cid"))
    if (!is.null(src$c_peak_group_X_s_peak_meta))
        src$c_peak_group_X_s_peak_meta <- cast(
            src$c_peak_group_X_s_peak_meta, ints = c("gXp_id", "pid",
                                                     "grpid"))
    src$s_peak_meta <- cast(
        src$s_peak_meta,
        ints = c("pid", "seqNum", "acquisitionNum", "precursorScanNum",
                 "precursorNearest", "grpid", "fileid", "sourceid",
                 "ms_level"),
        nums = c("precursorIntensity", "precursorMZ", "precursorRT",
                 "retentionTime", "aMz", "aPurity", "apkNm", "iMz",
                 "iPurity", "ipkNm", "inPkNm", "inPurity", "precursor_mz",
                 "retention_time", "collision_energy", "resolution",
                 "mass_accuracy", "mass_error", "retention_index"),
        flags = "purity_pass_flag",
        chrs = c("name", "accession", "polarity", "fragmentation_type",
                 "precursor_type", "instrument_type", "instrument",
                 "copyright", "column", "origin", "splash", "inchikey_id",
                 "spectrum_type"))
    src$s_peaks <- cast(
        src$s_peaks,
        ints = c("sid", "fileid", "scan", "cl", "count", "total", "grpid",
                 "pid"),
        nums = c("mz", "i", "snr", "ra", "rsd", "inPurity", "frac"),
        flags = grep("_flag$", names(src$s_peaks), value = TRUE),
        chrs = "type")
    if (!is.null(src$sm_matches))
        src$sm_matches <- cast(
            src$sm_matches, ints = c("mid", "lpid", "qpid"),
            nums = c("dpc", "rdpc", "cdpc", "mcount", "allcount",
                     "mpercent", "library_rt", "query_rt", "rtdiff",
                     "library_precursor_mz", "query_precursor_mz",
                     "library_precursor_ion_purity",
                     "query_precursor_ion_purity"),
            chrs = c("library_accession", "library_precursor_type",
                     "library_entry_name", "inchikey",
                     "library_source_name", "library_compound_name"))
    src
}

#' Empty strings as NA.
#'
#' @noRd
.mzs_na_empty <- function(x) {
    if (is.null(x))
        return(NULL)
    x <- as.character(x)
    x[!is.na(x) & !nzchar(trimws(x))] <- NA_character_
    x
}

#' Convert an msPurity SQLite database to a Parquet results dataset
#'
#' @description
#'
#' Converts a database written by [createDatabase()] to a Parquet results
#' dataset, as `createDatabase(format = "parquet")` would write it. The
#' database is not modified.
#'
#' MS2 scans are referenced, not copied. With `study`, each scan resolves to
#' its spectrum in that dataset; otherwise each source file is an external
#' source and scans are identified by mzML native id.
#'
#' Processing parameters are not stored in the database, which is recorded
#' in `tool_provenance` and the loss ledger. Unmapped tables are an error.
#'
#' Requires the suggested packages arrow and jsonlite.
#'
#' @param db `character(1)`, path of the SQLite database.
#'
#' @param path `character(1)`, the results dataset to create.
#'
#' @param study `character(1)` or `NULL`, path of the Parquet dataset the
#'     source files were converted into.
#'
#' @param studyKey `character(1)`, source key for the study.
#'
#' @param fileMap named `character`, `fileinfo.filepth` -> run id. Required
#'     when source file names collide or recorded paths no longer exist.
#'
#' @param overwrite `logical(1)`, whether to replace an existing dataset.
#'
#' @return The path of the dataset, invisibly.
#'
#' @seealso [createDatabase()], [validateParquet()]
#'
#' @examples
#' if (requireNamespace("arrow", quietly = TRUE) &&
#'     requireNamespace("jsonlite", quietly = TRUE)) {
#'     db <- system.file("extdata", "tests", "db",
#'                       "createDatabase_example.sqlite", package = "msPurity")
#'     out <- file.path(tempdir(), "converted.parquet")
#'     convertSqliteToParquet(db, out, overwrite = TRUE)
#'     validateParquet(out)
#' }
#' @export
convertSqliteToParquet <- function(db, path, study = NULL, studyKey = "study",
                                   fileMap = NULL, overwrite = FALSE) {
    .mzs_require("convertSqliteToParquet()")
    started <- .mzs_now()
    if (!is.character(db) || length(db) != 1L || !file.exists(db))
        stop("'db' must be the path of an msPurity SQLite database.",
             call. = FALSE)
    if (!is.null(study))
        study <- normalizePath(study, mustWork = TRUE)
    .mzs_check_destination(path, overwrite, study)
    src <- .mzs_read_sqlite(db)
    input <- .mzs_file_input(db)
    .mzs_mspurity_convert(
        src, path, study, studyKey, fileMap, overwrite,
        route = "sqlite-file", live = NULL,
        fn = "convertSqliteToParquet", started = started,
        inputs = list(input), parameters = list(),
        source_uri = db, checksum = input$sha256)
    invisible(normalizePath(path))
}

# ---------------------------------------------------------------------------
# Round trip SQLite -> Parquet, compared semantically.
#
# Both sides are normalised to one canonical table set: fixed column types,
# fixed row order, per-sample columns by position, msPurity ids recovered
# from source_identifier. Dropped positional columns are omitted; dropped
# recomputable ones are recomputed on the Parquet side.
# ---------------------------------------------------------------------------

.MZS_CANON_SCAN <- c("pid", "fileid", "seqNum", "acquisitionNum",
                     "precursorIntensity", "precursorMZ", "precursorRT",
                     "precursorScanNum", "retentionTime", "precursorNearest",
                     "aMz", "aPurity", "apkNm", "iMz", "iPurity", "ipkNm",
                     "inPkNm", "inPurity")
.MZS_CANON_AV_PEAK <- c("pid", "mz", "i", "snr", "ra", "rsd", "count",
                        "total", "inPurity", "cl", "snr_pass_flag",
                        "minnum_pass_flag", "minfrac_pass_flag",
                        "ra_pass_flag", "frac", "pass_flag")
.MZS_CANON_SCAN_PEAK <- c("pid", "mz", "i", "snr", "ra", "purity_pass_flag",
                          "intensity_pass_flag", "ra_pass_flag",
                          "snr_pass_flag", "pass_flag")

#' Order rows and reset row names.
#'
#' @noRd
.mzs_canon_order <- function(df, by) {
    if (nrow(df)) {
        o <- do.call(order, c(unname(as.list(df[by])),
                              list(method = "radix", na.last = TRUE)))
        df <- df[o, , drop = FALSE]
    }
    rownames(df) <- NULL
    df
}

#' Keep `cols`, coercing integer-valued and flag columns to fixed types.
#'
#' @noRd
.mzs_canon_cols <- function(df, cols, ints = character(),
                            flags = character()) {
    cols <- cols[cols %in% names(df)]
    df <- df[, cols, drop = FALSE]
    for (c in names(df)) {
        if (c %in% flags) df[[c]] <- as.logical(df[[c]])
        else if (c %in% ints) df[[c]] <- as.integer(round(as.numeric(df[[c]])))
        else if (!is.character(df[[c]])) df[[c]] <- as.numeric(df[[c]])
    }
    df
}

#' The canonical form of an msPurity database.
#'
#' @noRd
.mzs_canonical_sqlite <- function(db) {
    src <- .mzs_read_sqlite(db)
    fi <- .mzs_canon_cols(src$fileinfo, c("fileid", "filename", "filepth",
                                          "nm_save", "class"), "fileid")
    cp <- src$c_peaks
    names(cp) <- sub("^_+", "", names(cp))
    cp <- .mzs_canon_cols(cp, names(cp), c("cid", "fileid"))
    grp <- src$c_peak_groups
    sc <- .mzs_sample_columns(grp, nrow(fi))
    g <- .mzs_canon_cols(grp, c("grpid", "mz", "mzmin", "mzmax", "rt",
                                "rtmin", "rtmax", "npeaks", "grp_name",
                                "ms_level"), c("grpid", "npeaks", "ms_level"))
    for (j in seq_along(sc$samples))
        g[[paste0("sample_", j)]] <- suppressWarnings(as.numeric(
            grp[[sc$samples[j]]]))
    cxg <- .mzs_canon_cols(src$c_peak_X_c_peak_group,
                           c("grpid", "cid", "idi", "bestpeak"),
                           c("grpid", "cid", "idi", "bestpeak"))
    link <- .mzs_scan_links(src)
    link <- .mzs_canon_cols(link, c("pid", if (!all(is.na(link$cid))) "cid"
                                    else "grpid"), c("pid", "cid", "grpid"))
    scans <- .mzs_scans(src)
    sm <- .mzs_canon_cols(scans, c(.MZS_CANON_SCAN, "purity_pass_flag"),
                          c("pid", "fileid", "seqNum", "acquisitionNum",
                            "precursorScanNum", "precursorNearest"),
                          "purity_pass_flag")
    if ("purity_pass_flag" %in% names(sm) && all(is.na(sm$purity_pass_flag)))
        sm$purity_pass_flag <- NULL
    m <- src$s_peak_meta
    files <- data.frame(fileid = fi$fileid, run_id = as.character(fi$fileid),
                        data_origin = fi$filepth, source = NA_character_,
                        stringsAsFactors = FALSE)
    av <- .mzs_derived_index(src, files, "x")
    avm <- data.frame(pid = av$pid, method = av$method, grpid = av$grpid,
                      fileid = ifelse(av$method == "intra", av$fileid, NA),
                      precursor_mz = av$precursor_mz,
                      retention_time = av$retention_time,
                      inPurity = av$inPurity,
                      polarity = tolower(as.character(av$polarity)),
                      stringsAsFactors = FALSE)
    avm <- .mzs_canon_cols(avm, names(avm), c("pid", "grpid", "fileid"))
    pk <- src$s_peaks
    avp <- pk[pk$pid %in% av$pid, , drop = FALSE]
    avp <- .mzs_canon_cols(avp, .MZS_CANON_AV_PEAK,
                           c("pid", "count", "total", "cl"),
                           grep("_flag$", .MZS_CANON_AV_PEAK, value = TRUE))
    sp <- pk[pk$type %in% "scan", , drop = FALSE]
    scp <- NULL
    if (any(c("pass_flag", "snr") %in% names(sp)) &&
        any(!is.na(sp$snr %||% NA))) {
        sp$pid <- .mzs_scan_peak_pid(src, sp)
        scp <- .mzs_canon_cols(sp, .MZS_CANON_SCAN_PEAK, "pid",
                               grep("_flag$", .MZS_CANON_SCAN_PEAK,
                                    value = TRUE))
        scp$rank <- as.integer(stats::ave(seq_len(nrow(scp)), scp$pid,
                                          FUN = seq_along))
    }
    list(
        fileinfo = .mzs_canon_order(fi, "fileid"),
        c_peaks = .mzs_canon_order(cp[, sort(names(cp))], "cid"),
        c_peak_groups = .mzs_canon_order(g, "grpid"),
        c_peak_X_c_peak_group = .mzs_canon_order(cxg, c("grpid", "cid")),
        scan_link = .mzs_canon_order(link, names(link)),
        scans = .mzs_canon_order(sm, "pid"),
        averaged = .mzs_canon_order(avm, "pid"),
        averaged_peaks = .mzs_canon_order(avp, c("pid", "mz", "cl")),
        scan_peaks = if (!is.null(scp))
            .mzs_canon_order(scp, c("pid", "rank")))
}

#' The canonical form of a Parquet results dataset written from msPurity.
#'
#' @noRd
.mzs_canonical_parquet <- function(path) {
    m <- .mzs_read_manifest(path)
    tab <- function(n) .mzs_read_table(path, m, n)
    si <- tab("source_identifier")
    src_id <- function(table, column, keys) {
        s <- si[si$target_table == table & si$source_column == column, ]
        s$source_value[match(keys, s$target_key)]
    }
    num <- function(x) as.numeric(x)
    assay <- tab("assay")
    sample <- tab("sample")
    fi <- data.frame(
        fileid = num(src_id("assay", "fileid", assay$assay_id_)),
        filename = src_id("sample", "filename", assay$sample_id_),
        filepth = src_id("assay", "filepth", assay$assay_id_),
        nm_save = src_id("assay", "nm_save", assay$assay_id_),
        class = sample$sample_class[match(assay$sample_id_,
                                          sample$sample_id_)],
        stringsAsFactors = FALSE)
    fi <- .mzs_canon_cols(fi, names(fi), "fileid")
    fileid_of_run <- function(r) fi$fileid[match(r, assay$run_id)]

    pk <- tab("chromatographic_peak")
    cp <- data.frame(cid = num(src_id("chromatographic_peak", "cid",
                                      pk$chromatographic_peak_id_)),
                     fileid = fileid_of_run(pk$run_id),
                     mz = pk$exp_mass_to_charge, mzmin = pk$mass_to_charge_min,
                     mzmax = pk$mass_to_charge_max,
                     rt = pk$retention_time_in_seconds,
                     rtmin = pk$retention_time_in_seconds_start,
                     rtmax = pk$retention_time_in_seconds_end)
    for (c in grep("^x_mspurity_", names(pk), value = TRUE))
        cp[[sub("^x_mspurity_", "", c)]] <- pk[[c]]
    cp <- .mzs_canon_cols(cp, names(cp), c("cid", "fileid"))

    ft <- tab("feature")
    grpid_of <- function(fid) num(src_id("feature", "grpid", fid))
    g <- data.frame(grpid = grpid_of(ft$feature_id_),
                    mz = ft$exp_mass_to_charge, mzmin = ft$mass_to_charge_min,
                    mzmax = ft$mass_to_charge_max,
                    rt = ft$retention_time_in_seconds,
                    rtmin = ft$retention_time_in_seconds_start,
                    rtmax = ft$retention_time_in_seconds_end,
                    npeaks = ft$n_chromatographic_peaks,
                    grp_name = ft$feature_name, stringsAsFactors = FALSE)
    if ("x_mspurity_ms_level" %in% names(ft))
        g$ms_level <- ft$x_mspurity_ms_level
    g <- .mzs_canon_cols(g, names(g), c("grpid", "npeaks", "ms_level"))
    ab <- tab("abundance")
    for (j in sort(unique(ab$assay_id_))) {
        a <- ab[ab$assay_id_ == j, ]
        g[[paste0("sample_", j)]] <- a$value[match(
            ft$feature_id_[match(g$grpid, grpid_of(ft$feature_id_))],
            a$feature_id_)]
    }
    cid_of <- function(id) num(src_id("chromatographic_peak", "cid", id))
    cpf <- tab("chromatographic_peak_feature")
    cxg <- .mzs_canon_cols(data.frame(
        grpid = grpid_of(cpf$feature_id_),
        cid = cid_of(cpf$chromatographic_peak_id_),
        idi = cpf$x_mspurity_idi, bestpeak = as.numeric(cpf$is_representative)),
        c("grpid", "cid", "idi", "bestpeak"), c("grpid", "cid", "idi",
                                               "bestpeak"))

    sc <- tab("x_mspurity_scan")
    scan_pid <- num(src_id("x_mspurity_scan", "pid", sc$scan_annotation_id_))
    sm <- data.frame(pid = scan_pid, fileid = fileid_of_run(sc$ms2_run_id),
                     seqNum = sc$x_mspurity_seq_num)
    for (c in names(.MZS_SCAN_COLUMNS))
        if (.MZS_SCAN_COLUMNS[[c]] %in% names(sc))
            sm[[c]] <- sc[[.MZS_SCAN_COLUMNS[[c]]]]
    sm <- .mzs_canon_cols(sm, c(.MZS_CANON_SCAN, "purity_pass_flag"),
                          c("pid", "fileid", "seqNum", "acquisitionNum",
                            "precursorScanNum", "precursorNearest"),
                          "purity_pass_flag")
    if ("purity_pass_flag" %in% names(sm) && all(is.na(sm$purity_pass_flag)))
        sm$purity_pass_flag <- NULL

    sf <- tab("spectrum_feature")
    lk <- sf[sf$ms2_source != "self", ]
    lk_pid <- scan_pid[match(paste(lk$ms2_run_id, lk$x_mspurity_scan_number),
                             paste(sc$ms2_run_id, sc$scan_number))]
    link <- if (all(lk$link_mode == "chromatographic_peak"))
        data.frame(pid = lk_pid, cid = cid_of(lk$chromatographic_peak_id_))
    else data.frame(pid = lk_pid, grpid = grpid_of(lk$feature_id_))
    link <- .mzs_canon_cols(link, names(link), names(link))

    ## Averaged spectra.
    rr <- m$results$runs
    avm <- list()
    avp <- list()
    self <- sf[sf$ms2_source == "self", ]
    for (r in m$runs) {
        sp <- as.data.frame(arrow::read_parquet(file.path(
            .mzs_run_dir(path, r), "part-0.parquet")))
        if (!nrow(sp))
            next
        method <- rr[[r$run_id]]$method
        pid <- num(src_id("spectra", "pid", sp$spectrum_id_))
        grpid <- grpid_of(self$feature_id_[match(sp$spectrum_id_,
                                                 self$ms2_spectrum_id_)])
        fid <- if (method == "intra")
            fi$fileid[match(rr[[r$run_id]]$source_run_id, assay$run_id)]
        else NA
        pol <- c("-1" = "negative", "1" = "positive")[as.character(
            sp$scan_polarity)]
        avm[[length(avm) + 1L]] <- data.frame(
            pid = pid, method = method, grpid = grpid,
            fileid = rep(fid, nrow(sp)),
            precursor_mz = sp$selected_ion_mz,
            retention_time = ft$retention_time_in_seconds[match(
                self$feature_id_[match(sp$spectrum_id_,
                                       self$ms2_spectrum_id_)],
                ft$feature_id_)],
            inPurity = sp$x_mspurity_precursor_purity,
            polarity = unname(pol), stringsAsFactors = FALSE)
        n <- lengths(sp$mz)
        lst <- function(c) if (c %in% names(sp))
            unlist(lapply(seq_len(nrow(sp)), function(i)
                if (is.null(sp[[c]][[i]])) rep(NA, n[i]) else sp[[c]][[i]]))
            else rep(NA, sum(n))
        p <- data.frame(pid = rep(pid, n), mz = unlist(sp$mz),
                        i = unlist(sp$intensity), snr = lst("sn"),
                        ra = lst("x_mspurity_ra"), rsd = lst("intensity_rsd"),
                        count = lst("contributor_count"),
                        total = lst("contributor_total"),
                        inPurity = lst("x_mspurity_in_purity"),
                        cl = lst("cluster"),
                        snr_pass_flag = lst("x_mspurity_snr_pass_flag"),
                        minnum_pass_flag = lst("x_mspurity_minnum_pass_flag"),
                        minfrac_pass_flag = lst("x_mspurity_minfrac_pass_flag"),
                        ra_pass_flag = lst("x_mspurity_ra_pass_flag"))
        ## Not stored; recomputed.
        p$frac <- p$count / p$total
        p$pass_flag <- p$minfrac_pass_flag & p$snr_pass_flag &
            p$ra_pass_flag & p$minnum_pass_flag
        avp[[length(avp) + 1L]] <- p
    }
    avm <- do.call(rbind, avm)
    avm <- .mzs_canon_cols(avm, names(avm), c("pid", "grpid", "fileid"))
    avp <- do.call(rbind, avp)
    avp <- .mzs_canon_cols(avp, .MZS_CANON_AV_PEAK,
                           c("pid", "count", "total", "cl"),
                           grep("_flag$", .MZS_CANON_AV_PEAK, value = TRUE))

    scp <- NULL
    if (!is.null(m$results$tables$x_mspurity_scan_peak)) {
        s <- tab("x_mspurity_scan_peak")
        scp <- data.frame(
            pid = scan_pid[match(s$scan_annotation_id_,
                                 sc$scan_annotation_id_)],
            mz = s$mz, i = s$intensity, snr = s$snr,
            ra = s$relative_intensity,
            purity_pass_flag = s$purity_pass_flag,
            intensity_pass_flag = s$intensity_pass_flag,
            ra_pass_flag = s$ra_pass_flag, snr_pass_flag = s$snr_pass_flag)
        scp$pass_flag <- scp$purity_pass_flag & scp$intensity_pass_flag &
            scp$ra_pass_flag & scp$snr_pass_flag
        scp <- .mzs_canon_cols(scp, .MZS_CANON_SCAN_PEAK, "pid",
                               grep("_flag$", .MZS_CANON_SCAN_PEAK,
                                    value = TRUE))
        scp$rank <- s$peak_rank
    }
    list(
        fileinfo = .mzs_canon_order(fi, "fileid"),
        c_peaks = .mzs_canon_order(cp[, sort(names(cp))], "cid"),
        c_peak_groups = .mzs_canon_order(g, "grpid"),
        c_peak_X_c_peak_group = .mzs_canon_order(cxg, c("grpid", "cid")),
        scan_link = .mzs_canon_order(link, names(link)),
        scans = .mzs_canon_order(sm, "pid"),
        averaged = .mzs_canon_order(avm, "pid"),
        averaged_peaks = .mzs_canon_order(avp, c("pid", "mz", "cl")),
        scan_peaks = if (!is.null(scp))
            .mzs_canon_order(scp, c("pid", "rank")))
}
