# Read only the projection a validation query needs. DuckDB spills joins
# and aggregations to a private temporary directory when its budget fills.
.mzs_validation_store <- function(path, manifest, dir) {
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(dir, "validation"),
        config = list(threads = "1", memory_limit = "128MB"))
    id <- function(x) as.character(DBI::dbQuoteIdentifier(con, x))
    str <- function(x) as.character(DBI::dbQuoteString(con, x))
    query <- function(sql) DBI::dbGetQuery(con, sql)
    views <- new.env(parent = emptyenv())
    view <- function(name) {
        if (!exists(name, views, inherits = FALSE)) {
            e <- manifest$results$tables[[name]]
            files <- list.files(file.path(path, e$path), "\\.parquet$",
                                full.names = TRUE, recursive = TRUE)
            if (length(files)) {
                sql <- paste0("SELECT * FROM read_parquet([", paste(str(files), collapse = ","),
                    "], union_by_name=true, hive_types_autocast=false, hive_partitioning=",
                    if (length(e$partitioning)) "true" else "false", ")")
                DBI::dbExecute(con, paste("CREATE TEMP VIEW", id(name), "AS", sql))
            } else {
                DBI::dbWriteTable(con, name, .mzs_empty_table(name), temporary = TRUE)
            }
            assign(name, TRUE, views)
        }
        id(name)
    }
    rows <- function(name) query(paste("SELECT count(*) AS n FROM", view(name)))$n
    schema <- function(name) query(paste("SELECT * FROM", view(name), "LIMIT 0"))
    batches <- function(name, visit) {
        rs <- DBI::dbSendQuery(con, paste("SELECT * FROM", view(name)))
        on.exit(DBI::dbClearResult(rs))
        repeat {
            batch <- DBI::dbFetch(rs, n = .mzs_batch_size())
            if (!nrow(batch)) break
            visit(batch)
        }
    }
    empty <- function(name) {
        t <- view(name)
        cols <- names(schema(name))
        sql <- paste0("count(", id(cols), ") AS ", id(cols))
        n <- query(paste("SELECT", paste(sql, collapse = ","), "FROM", t))
        names(n)[vapply(n, function(v) v == 0, logical(1))]
    }
    list(view = view, query = query, id = id, str = str, schema = schema,
         rows = rows, batches = batches, empty = empty, con = con)
}

# Relational checks for the result families createDatabase writes. Joins
# include every part of each table, including references across batches.
.mzs_v_integrity <- function(path, m, read, have) {
    db <- attr(read, "store")
    t <- db$view
    f <- list()
    check <- function(table, condition, message, req = "R-099", level = "MUST") {
        if (nrow(db$query(paste("SELECT 1 FROM", condition, "LIMIT 1"))))
            f[[length(f) + 1L]] <<- .mzs_finding(level, req, table, message)
    }
    orphan <- function(table, field, target, key = field) {
        if (!have(table) || !field %in% names(db$schema(table))) return()
        target_test <- if (have(target)) paste0(" AND NOT EXISTS (SELECT 1 FROM ",
            t(target), " y WHERE x.", field, " = y.", key, ")") else ""
        check(table, paste0(t(table), " x WHERE x.", field, " IS NOT NULL", target_test),
              paste(field, "names a row absent from", target))
    }
    for (table in c("chromatographic_peak_feature", "abundance", "spectrum_feature",
                    "x_mspurity_msn_feature", "evidence", "x_mspurity_feature_annotation",
                    "x_mspurity_combined_annotation"))
        orphan(table, "feature_id_", "feature")
    for (table in c("chromatographic_peak_feature", "spectrum_feature", "x_mspurity_msn_feature"))
        orphan(table, "chromatographic_peak_id_", "chromatographic_peak")
    for (table in c("evidence", "compound_xref", "compound_synonym",
                    "x_mspurity_combined_annotation"))
        orphan(table, "compound_id_", "compound")
    orphan("evidence", "software_id_", "software")
    orphan("evidence_score", "evidence_id_", "evidence")
    orphan("abundance", "assay_id_", "assay")
    orphan("assay", "sample_id_", "sample")
    for (table in c("x_mspurity_scan_peak", "x_mspurity_msn_feature"))
        orphan(table, "scan_annotation_id_", "x_mspurity_scan")
    for (table in c("x_mspurity_msn_feature", "x_mspurity_scan"))
        orphan(table, "root_scan_annotation_id_", "x_mspurity_scan", "scan_annotation_id_")
    orphan("x_mspurity_scan", "parent_scan_annotation_id_", "x_mspurity_scan", "scan_annotation_id_")

    acts <- vapply(m$provenance, `[[`, "", "id")
    for (table in names(m$results$tables)) {
        if (!have(table)) next
        tbl <- t(table)
        check(table, paste(tbl, "WHERE activity NOT IN (",
            paste(db$str(acts), collapse = ","), ")"), "Unknown provenance activity.")
        key <- if (table %in% names(.MZS_TABLES)) .mzs_table_schema(table)$key else NA_character_
        if (!is.na(key))
            check(table, paste(tbl, "HAVING count(*) <> count(DISTINCT", db$id(key),
                ") OR min(", db$id(key), ") <> 1 OR max(", db$id(key), ") <> count(*)"),
                "Surrogate keys must be contiguous from one.")
        cols <- unlist(m$results$tables[[table]]$sorted_by)
        if (length(cols))
            check(table, paste(tbl, "GROUP BY", paste(db$id(cols), collapse = ","),
                "HAVING count(*) > 1"), "Rows tie in the declared sort order.", "R-038")
    }
    if (have("spectrum_feature")) {
        tbl <- t("spectrum_feature")
        check("spectrum_feature", paste(tbl, "WHERE (link_mode = 'chromatographic_peak')
            <> (chromatographic_peak_id_ IS NOT NULL)"), "Invalid chromatographic peak reference.")
        check("spectrum_feature", paste(tbl, "WHERE link_mode NOT IN (",
            paste(db$str(.MZS_LINK_MODES), collapse = ","), ")"), "Unknown link mode.", "R-044")
        check("spectrum_feature", paste(tbl, "WHERE link_mode = 'feature_width' AND run_scoped"),
              "Feature-width association is run-scoped.", "R-045")
    }
    if (have("x_mspurity_scan") && "ms_level" %in% names(db$schema("x_mspurity_scan"))) {
        tbl <- t("x_mspurity_scan")
        check("x_mspurity_scan", paste(tbl, "WHERE ms_level NOT IN (2,3,4)"), "Invalid scan level.")
        for (prefix in c("parent", "root")) {
            join <- paste0(tbl, " x JOIN ", tbl, " y ON x.", prefix,
                "_scan_annotation_id_ = y.scan_annotation_id_ WHERE ")
            wrong <- if (prefix == "parent") "y.ms_level <> x.ms_level - 1" else "y.ms_level <> 2"
            cols <- c("source", "run_id", "spectrum_id_", "native_id")
            ref <- paste0("x.", prefix, "_", cols, " IS DISTINCT FROM y.ms2_", cols)
            check("x_mspurity_scan", paste0(join, paste(c(wrong, ref), collapse = " OR ")),
                  paste("Inconsistent", prefix, "level or reference."))
        }
        states <- c("valid", "root", "missing_reference", "unresolved_reference",
                    "invalid_reference", "missing_stage", "missing_precursor")
        for (col in c("parent_status", "path_status"))
            check("x_mspurity_scan", paste(tbl, "WHERE", col, "NOT IN (",
                paste(db$str(states), collapse = ","), ")"), paste("Unknown", col), level = "SHOULD")
        if (have("x_mspurity_msn_feature"))
            check("x_mspurity_msn_feature", paste(t("x_mspurity_msn_feature"),
                "x JOIN", tbl, "y USING (scan_annotation_id_) WHERE x.ms_level <> y.ms_level OR
                x.root_scan_annotation_id_ <> y.root_scan_annotation_id_"),
                "An MSn feature link disagrees with its scan.")
    }
    if (have("source_identifier")) {
        tbl <- t("source_identifier")
        targets <- db$query(paste("SELECT DISTINCT target_table FROM", tbl))$target_table
        for (target in targets) {
            cond <- if (target == "spectra")
                paste("x.target_key < 1 OR x.target_key >", sum(.mzs_runs_frame(m)$n_spectra))
            else if (have(target) && !is.na(.mzs_table_schema(target)$key))
                paste0("NOT EXISTS (SELECT 1 FROM ", t(target), " y WHERE x.target_key = y.",
                       db$id(.mzs_table_schema(target)$key), ")") else "TRUE"
            check("source_identifier", paste(tbl, "x WHERE target_table =", db$str(target),
                "AND (", cond, ")"), paste("Invalid target key for", target))
        }
    }
    if (have("merge_member")) {
        mm <- t("merge_member")
        check("merge_member", paste(mm, "WHERE member_source = 'self' AND
            member_spectrum_id_ = merged_spectrum_id_"), "An average is its own member.")
        if (length(m$runs)) {
            files <- unlist(lapply(m$runs, function(r)
                list.files(.mzs_run_dir(path, r), "\\.parquet$", recursive = TRUE, full.names = TRUE)))
            meta <- paste0("read_parquet([", paste(db$str(files), collapse = ","),
                           "], union_by_name=true)")
            check("merge_member", paste(mm, "m JOIN", meta, "x ON m.merged_spectrum_id_ = x.spectrum_id_
                JOIN", meta, "y ON m.member_spectrum_id_ = y.spectrum_id_
                WHERE m.member_source = 'self' AND x.ms_level <> y.ms_level"), "Member MS levels differ.")
        }
        if (!have("x_mspurity_merge_peak") &&
            "x_mspurity_contributes" %in% names(db$schema("merge_member")))
            check("merge_member", paste(mm, "WHERE x_mspurity_contributes"),
                  "A contributing member has no peak edges.")
        if (have("x_mspurity_merge_peak")) {
            edge <- t("x_mspurity_merge_peak")
            join <- "m.merged_spectrum_id_ = e.merged_spectrum_id_ AND m.member_rank = e.member_rank"
            check("x_mspurity_merge_peak", paste(edge, "e LEFT JOIN", mm, "m ON", join,
                "WHERE m.member_rank IS NULL OR e.basis NOT IN ('source_reference','filter_input','derived','unknown')
                OR (e.member_peak_rank IS NULL) <> (e.basis = 'unknown')
                OR (m.member_source = 'self') <> (e.basis = 'derived')
                OR e.merged_peak_rank < 1 OR e.member_peak_rank < 1 OR e.member_rank < 1"),
                "Invalid peak rank, basis or membership reference.")
            if ("x_mspurity_contributes" %in% names(db$schema("merge_member")))
                check("merge_member", paste(mm, "m WHERE x_mspurity_contributes <> EXISTS
                    (SELECT 1 FROM", edge, "e WHERE", join, ")"), "Member contribution disagrees with its edges.")
        }
    }
    do.call(.mzs_findings, f)
}
