# Writing Parquet datasets: native runs, results tables, the results index
# and the manifest.
#
# A new dataset is built and validated in a staging directory beside its
# destination, then moved into place. The manifest is written last, so an
# interrupted write never looks like a dataset. Updates add new Parquet parts
# or a new table revision, never rewriting existing parts, then commit a new
# manifest atomically if no other writer has committed meanwhile.

.MZS_ROW_GROUP <- 65536L
.MZS_SPECTRA_ROW_GROUP <- 250L

# ---------------------------------------------------------------------------
# Preparing tables.
# ---------------------------------------------------------------------------

#' Cast one table to its schema, add its `activity`, and sort it.
#'
#' @param df data.frame. Columns outside the schema must be named
#'     `x_<tool>_<name>`; a missing nullable column is added as a typed null.
#'
#' @return the typed data.frame: schema columns in schema order, then
#'     passthrough columns, rows in the declared sort order.
#'
#' @noRd
.mzs_prepare_table <- function(name, df, activity) {
    s <- .mzs_table_schema(name)
    cols <- s$columns
    df <- as.data.frame(df, stringsAsFactors = FALSE, optional = TRUE)
    n <- nrow(df)
    df$activity <- rep(activity, n)

    extra <- setdiff(names(df), cols$name)
    ## Drop wholly-null passthrough columns.
    empty <- extra[vapply(extra, function(c) {
        v <- df[[c]]
        !is.list(v) && all(is.na(v) & !is.nan(v))
    }, logical(1))]
    if (n && length(empty)) {
        df <- df[, setdiff(names(df), empty), drop = FALSE]
        extra <- setdiff(extra, empty)
    }
    bad <- extra[!grepl("^x_[A-Za-z0-9]+_[A-Za-z0-9_]+$", extra)]
    if (length(bad))
        stop("Table '", name, "' has column(s) mzStack-4 does not define: ",
             paste(bad, collapse = ", "), ".", call. = FALSE)
    bad <- .mzs_bad_names(extra, rep("", length(extra)))
    if (length(bad))
        stop("Column name(s) of table '", name, "' break the naming rule ",
             "of mzStack-0 G-016: ", paste(bad, collapse = ", "), ".",
             call. = FALSE)

    for (i in seq_len(nrow(cols))) {
        nm <- cols$name[i]
        v <- if (nm %in% names(df)) df[[nm]] else rep(NA, n)
        v <- .mzs_as_r_type(v, cols$type[i])
        # NaN and Inf are values, not missing markers.
        miss <- is.na(v) & !is.nan(v)
        if (!cols$nullable[i] && any(miss))
            stop("Column '", nm, "' of table '", name, "' is NOT NULL but ",
                 "has ", sum(miss), " missing value(s).", call. = FALSE)
        df[[nm]] <- v
    }
    df <- df[, c(cols$name, extra), drop = FALSE]
    if (n) {
        o <- do.call(order, c(unname(as.list(df[s$sorted_by])),
                              list(method = "radix", na.last = TRUE)))
        df <- df[o, , drop = FALSE]
    }
    rownames(df) <- NULL
    if (n > 1L && anyDuplicated(df[s$sorted_by]))
        .mzs_abort("format", "The declared sort order of '", name, "' (",
                   paste(s$sorted_by, collapse = ", "), ") is not total over ",
                   "the rows written: two rows tie (R-038).")
    df
}

#' Check a table's own surrogate key runs `offset + 1 .. offset + n`.
#'
#' @noRd
.mzs_check_keys <- function(name, df, offset = 0) {
    key <- .mzs_table_schema(name)$key
    if (is.na(key) || !nrow(df))
        return(invisible())
    want <- offset + seq_len(nrow(df))
    if (!identical(sort(as.numeric(df[[key]])), as.numeric(want)))
        stop("Surrogate key '", key, "' of table '", name, "' must run ",
             offset + 1, "..", offset + nrow(df), ".", call. = FALSE)
    invisible()
}

#' The Arrow table for a prepared data.frame.
#'
#' @noRd
.mzs_to_arrow <- function(name, df) {
    cols <- .mzs_table_schema(name)$columns
    arrays <- lapply(names(df), function(nm) {
        i <- match(nm, cols$name)
        x <- df[[nm]]
        if (is.na(i))
            return(arrow::Array$create(x))
        code <- cols$type[i]
        if (code == "dict")
            return(arrow::Array$create(factor(x))$cast(
                .mzs_arrow_type("dict")))
        arrow::Array$create(x, type = .mzs_arrow_type(code))
    })
    names(arrays) <- names(df)
    do.call(arrow::arrow_table, arrays)
}

#' Write a prepared table's rows as new Parquet parts under `root`.
#'
#' Part names carry the generation and a random suffix, so an append never
#' touches an existing part.
#'
#' @return the paths of the parts written, relative to `dir`.
#'
#' @noRd
.mzs_write_parts <- function(dir, root, name, df, generation) {
    s <- .mzs_table_schema(name)
    part <- sprintf("part-g%06d-%s.parquet", as.integer(generation),
                    .mzs_rand())
    groups <- if (length(s$partitioning) && nrow(df)) {
        key <- s$partitioning
        v <- as.character(df[[key]])
        if (anyNA(v))
            stop("Partition column '", key, "' of table '", name, "' has ",
                 "missing values.", call. = FALSE)
        lapply(sort(unique(v), method = "radix"), function(u)
            list(dir = file.path(root, paste0(key, "=", u)),
                 rows = which(v == u)))
    } else {
        list(list(dir = root, rows = seq_len(nrow(df))))
    }
    written <- character()
    for (g in groups) {
        dir.create(file.path(dir, g$dir), recursive = TRUE,
                   showWarnings = FALSE)
        rel <- file.path(g$dir, part)
        sub <- df[g$rows, setdiff(names(df), s$partitioning), drop = FALSE]
        arrow::write_parquet(.mzs_to_arrow(name, sub), file.path(dir, rel),
                             chunk_size = .MZS_ROW_GROUP,
                             compression = "zstd")
        written <- c(written, rel)
    }
    written
}

#' The `results.tables` entry of a freshly written table.
#'
#' @noRd
.mzs_table_entry <- function(name, root, rows, parts, generation,
                             revision = 1L, supersedes = NULL,
                             fill_policy = NULL, selection = NULL) {
    s <- .mzs_table_schema(name)
    e <- list(entity = name, path = root, rows = as.integer(rows),
              parts = as.integer(parts), generation = as.integer(generation),
              revision = as.integer(revision),
              sorted_by = .mzs_array(s$sorted_by))
    if (!is.null(supersedes))
        e$supersedes <- as.integer(supersedes)
    if (length(s$partitioning))
        e$partitioning <- .mzs_array(s$partitioning)
    if (!is.null(fill_policy))
        e$fill_policy <- fill_policy
    if (!is.null(selection))
        e$selection <- selection
    e
}

# ---------------------------------------------------------------------------
# Native runs.
# ---------------------------------------------------------------------------

#' Write native runs of derived or library spectra.
#'
#' @param runs list of runs in manifest order, each `list(run_id, meta,
#'     peaks)`: `meta` a data.frame of storage-vocabulary scalar columns, one
#'     row per spectrum in write order; `peaks` a named list of list columns
#'     (`mz`, `intensity` and any peak-annotation columns), each a list with
#'     one element per spectrum and an Arrow element type in its `type`
#'     attribute.
#'
#' @param types named character vector: Arrow type code of each `meta`
#'     column.
#'
#' @return the run entries for the manifest.
#'
#' @noRd
.mzs_write_native_runs <- function(dir, runs, types, generation) {
    base <- 1L
    entries <- list()
    for (r in runs) {
        .mzs_check_key(r$run_id, "A run_id")
        n <- nrow(r$meta)
        rel <- file.path("spectra", paste0("run_id=", r$run_id))
        dir.create(file.path(dir, rel), recursive = TRUE, showWarnings = FALSE)
        arrays <- list(
            spectrum_id_ = arrow::Array$create(base + seq_len(n) - 1L,
                                               type = arrow::int32()),
            spectrum_index = arrow::Array$create(seq_len(n) - 1,
                                                 type = arrow::int64()))
        for (c in names(r$meta)) {
            code <- types[[c]]
            arrays[[c]] <- arrow::Array$create(
                .mzs_as_r_type(r$meta[[c]], code),
                type = .mzs_arrow_type(code))
        }
        for (c in names(r$peaks)) {
            el <- attr(r$peaks[[c]], "type") %||% "f64"
            v <- lapply(unclass(r$peaks[[c]]), function(e)
                if (is.null(e)) NULL else .mzs_as_r_type(e, el))
            arrays[[c]] <- arrow::Array$create(
                v, type = arrow::list_of(.mzs_arrow_type(el)))
        }
        arrow::write_parquet(do.call(arrow::arrow_table, arrays),
                             file.path(dir, rel, "part-0.parquet"),
                             chunk_size = .MZS_SPECTRA_ROW_GROUP,
                             compression = "snappy")
        entries[[length(entries) + 1L]] <- list(
            run_id = r$run_id, kind = "native", path = rel,
            n_spectra = as.integer(n), uid_base = as.integer(base),
            ingested_at = as.integer(generation),
            signal = list(layout = "list", profile = NULL, centroid = NULL),
            projections = .mzs_object())
        base <- base + n
    }
    entries
}

#' A list column for `.mzs_write_native_runs()`.
#'
#' @noRd
.mzs_peak_list <- function(x, type = "f64") structure(x, type = type)

# ---------------------------------------------------------------------------
# Committing.
# ---------------------------------------------------------------------------

#' Check a dataset's destination before anything is written.
#'
#' The destination may not overlap a source dataset.
#'
#' @noRd
.mzs_check_destination <- function(path, overwrite, sources = character()) {
    if (!is.character(path) || length(path) != 1L || !nzchar(path))
        stop("'path' must be a single directory path.", call. = FALSE)
    if (!dir.exists(dirname(path)))
        stop("The parent directory of '", path, "' does not exist.",
             call. = FALSE)
    full <- .mzs_full_path(path)
    for (s in sources[!is.na(sources)]) {
        sf <- .mzs_full_path(s)
        if (identical(full, sf) || startsWith(full, paste0(sf, "/")) ||
            startsWith(sf, paste0(full, "/")))
            .mzs_abort("semantic", "'", path, "' would overlap the source ",
                       "dataset '", s, "'; a results dataset never writes ",
                       "into a dataset it reads (R-069).")
    }
    if (file.exists(path) && !isTRUE(overwrite))
        stop("'", path, "' already exists. Use overwrite = TRUE to replace ",
             "it.", call. = FALSE)
    if (file.exists(path) && !dir.exists(path))
        stop("'", path, "' exists and is not a directory.", call. = FALSE)
    invisible(full)
}

#' A staging directory beside the destination, on the same filesystem so
#' that moving it into place is a rename.
#'
#' @noRd
.mzs_stage <- function(path) {
    stage <- file.path(dirname(path), paste0(".", basename(path),
                                             ".mzs-staging-", Sys.getpid(),
                                             "-", .mzs_rand()))
    dir.create(stage)
    stage
}

#' Move a completed staging directory into place.
#'
#' @noRd
.mzs_finalise <- function(stage, path, overwrite) {
    old <- NULL
    if (file.exists(path)) {
        if (!isTRUE(overwrite))
            stop("'", path, "' appeared while it was being written.",
                 call. = FALSE)
        old <- paste0(path, ".old-", .mzs_rand())
        if (!file.rename(path, old))
            stop("Could not move the existing '", path, "' aside.",
                 call. = FALSE)
    }
    if (!file.rename(stage, path)) {
        if (!is.null(old))
            file.rename(old, path)
        stop("Could not move the new dataset into '", path, "'.",
             call. = FALSE)
    }
    if (!is.null(old))
        unlink(old, recursive = TRUE)
    invisible(path)
}

#' Write a new dataset.
#'
#' @param path destination directory.
#'
#' @param plan list with elements:
#'   - `role`: `"results"` or `"library"`;
#'   - `runs`: native runs for `.mzs_write_native_runs()`, and `run_types`
#'     the type codes of their metadata columns;
#'   - `results_runs`: named list, run_id -> `list(method, scope, source,
#'     source_run_id)`;
#'   - `peak_columns`: named list, run_id -> its peak-annotation columns;
#'   - `tables`: named list of data.frames;
#'   - `omitted`: OPTIONAL tables deliberately not written;
#'   - `fill_policy`, `selection` (named list, table -> selection);
#'   - `sources`: source entries;
#'   - `library`: `list(version, digest)` for a library dataset;
#'   - `activity`: `list(action, started, fn, inputs, parameters)`.
#'
#' @return the committed manifest, invisibly.
#'
#' @noRd
.mzs_commit_new <- function(path, plan, overwrite = FALSE) {
    stage <- .mzs_stage(path)
    on.exit(unlink(stage, recursive = TRUE), add = TRUE)
    generation <- 1L
    act <- "act-0001"

    m <- .mzs_new_manifest()
    m$uid <- plan$uid %||% .mzs_uid()
    m$role <- plan$role %||% "results"
    if (!is.null(plan$library))
        m$library <- plan$library
    m$runs <- .mzs_write_native_runs(stage, plan$runs %||% list(),
                                     plan$run_types %||% character(),
                                     generation)

    entries <- list()
    index <- list()
    for (name in names(plan$tables)) {
        df <- plan$tables[[name]]
        .mzs_check_keys(name, df)
        df <- .mzs_prepare_table(name, df, act)
        root <- file.path("results", name, "rev-1")
        parts <- .mzs_write_parts(stage, root, name, df, generation)
        entries[[name]] <- .mzs_table_entry(
            name, root, nrow(df), length(parts), generation,
            fill_policy = if (name == "abundance") plan$fill_policy,
            selection = plan$selection[[name]])
        index[[length(index) + 1L]] <- .mzs_index_table(name, names(df))
    }
    for (name in plan$omitted %||% character())
        entries[[name]] <- list(entity = name, omitted = TRUE)
    for (r in names(plan$peak_columns %||% list()))
        if (length(plan$peak_columns[[r]]))
            index[[length(index) + 1L]] <- .mzs_index_peaks(
                r, plan$peak_columns[[r]])

    if (length(plan$sources))
        m$sources <- .mzs_merge_sources(list(), plan$sources)
    a <- plan$activity
    m$provenance <- list(.mzs_activity(
        id = act, generation = generation,
        started = a$started %||% .mzs_now(), action = a$action,
        agent = .mzs_agent(a$fn), inputs = a$inputs %||% list(),
        outputs = list(runs = vapply(m$runs, function(r) r$run_id,
                                     character(1)),
                       tables = names(plan$tables)),
        parameters = a$parameters))
    if (identical(m$role, "results")) {
        rr <- plan$results_runs %||% list()
        m$results <- list(
            version = .MZS_RESULTS_VERSION,
            column_mapping_path = .MZS_COLUMN_MAP_PATH,
            runs = if (length(rr)) rr else .mzs_object(),
            tables = if (length(entries)) entries else .mzs_object())
        .mzs_write_json(.mzs_index_merge(NULL, index),
                        file.path(stage, .MZS_COLUMN_MAP_PATH))
    }

    .mzs_write_manifest(stage, m)
    report <- validateParquet(stage, sourcePaths = plan$sourcePaths)
    must <- report[report$level == "MUST", , drop = FALSE]
    if (nrow(must))
        .mzs_abort("format", "Refusing to commit an invalid dataset:\n",
                   paste0("  [", must$requirement, "] ",
                          ifelse(is.na(must$table), "", paste0(must$table,
                                                               ": ")),
                          must$message, collapse = "\n"))
    if (isTRUE(getOption("msPurity.parquet.fail_before_commit")))
        stop("msPurity.parquet.fail_before_commit is set: stopping before ",
             "the dataset is moved into place.", call. = FALSE)
    .mzs_finalise(stage, path, overwrite)
    invisible(m)
}

#' Add to an existing results dataset: append rows to tables, supersede
#' tables with a new revision, record one new activity, and commit.
#'
#' Appending writes new parts beside the existing ones; superseding writes a
#' new revision to a new path. New parts are removed if the commit fails.
#' Surrogate keys continue from the table's declared row count.
#'
#' @param append named list of data.frames, keys numbered from 1; offset
#'     here.
#'
#' @param supersede named list of data.frames: each replaces its table.
#'
#' @param selection named list, table -> selection declaration. Appending
#'     under a different declaration is Unsupported.
#'
#' @return the committed manifest, invisibly.
#'
#' @noRd
.mzs_commit_update <- function(path, append = list(), supersede = list(),
                               selection = list(), sources = list(),
                               activity, sourcePaths = NULL) {
    m <- .mzs_read_manifest(path)
    if (!identical(m$role, "results"))
        .mzs_abort("unsupported", "'", path, "' is not a results dataset.")
    if (!identical(.mzs_semver_major(m$results$version), .MZS_RESULTS_MAJOR))
        .mzs_abort("unsupported", "'", path, "' holds results version ",
                   m$results$version %||% "<missing>", "; msPurity writes ",
                   .MZS_RESULTS_MAJOR, ".x (R-035).")
    expect <- m$generation
    generation <- expect + 1L
    act <- .mzs_next_activity_id(m)
    ## Add ledger rows for wholly-null columns of the tables written here.
    tb0 <- m$results$tables %||% list()
    if (!is.null(tb0$loss_ledger) && !isTRUE(tb0$loss_ledger$omitted)) {
        have <- .mzs_read_table(path, m, "loss_ledger")$target_ref
        new <- c(append, supersede)
        new <- new[setdiff(names(new), "loss_ledger")]
        fn <- activity$fn %||% "This step"
        losses <- Filter(function(l) !l$target_ref %in% have,
                         .mzs_null_column_losses(new, list()))
        losses <- lapply(losses, function(l) {
            l$reason <- paste0(fn, "() does not record this quantity.")
            l
        })
        given <- append$loss_ledger
        losses <- c(lapply(seq_len(NROW(given)), function(i)
            list(construct = given$source_construct[i],
                 disposition = given$disposition[i],
                 affected_rows = given$affected_rows[i],
                 target_ref = given$target_ref[i],
                 reason = given$reason[i])), losses)
        if (length(losses))
            append$loss_ledger <- .mzs_loss_ledger(losses)
    }
    written <- character()
    dirs <- character()
    on.exit({
        unlink(file.path(path, written))
        unlink(file.path(path, dirs), recursive = TRUE)
    }, add = TRUE)
    tables <- m$results$tables %||% .mzs_object()
    index <- list()

    for (name in names(append)) {
        df <- append[[name]]
        old <- tables[[name]]
        have <- !is.null(old) && !isTRUE(old$omitted)
        if (have && !is.null(selection[[name]]) &&
            !.mzs_json_equal(old$selection %||% list(policy = "complete"),
                             selection[[name]]))
            .mzs_abort("unsupported", "Table '", name, "' is declared with ",
                       "selection ", .mzs_to_json(old$selection, FALSE),
                       "; rows selected under another policy cannot be ",
                       "added to it.")
        offset <- if (have) as.numeric(old$rows) else 0
        key <- .mzs_table_schema(name)$key
        if (!is.na(key) && nrow(df))
            df[[key]] <- df[[key]] + offset
        .mzs_check_keys(name, df, offset)
        df <- .mzs_prepare_table(name, df, act)
        root <- if (have) old$path else file.path("results", name, "rev-1")
        if (!have)
            dirs <- c(dirs, root)
        parts <- .mzs_write_parts(path, root, name, df, generation)
        written <- c(written, parts)
        e <- if (have) old else .mzs_table_entry(name, root, 0L, 0L,
                                                 generation)
        e$rows <- as.integer(offset + nrow(df))
        e$parts <- as.integer((e$parts %||% 0L) + length(parts))
        e$generation <- generation
        if (!is.null(selection[[name]]))
            e$selection <- selection[[name]]
        tables[[name]] <- e
        index[[length(index) + 1L]] <- .mzs_index_table(name, names(df))
    }
    for (name in names(supersede)) {
        df <- supersede[[name]]
        old <- tables[[name]]
        rev <- if (is.null(old) || isTRUE(old$omitted)) 1L
               else as.integer(old$revision) + 1L
        .mzs_check_keys(name, df)
        df <- .mzs_prepare_table(name, df, act)
        root <- file.path("results", name, paste0("rev-", rev))
        if (file.exists(file.path(path, root)))
            stop("Revision directory '", root, "' already exists.",
                 call. = FALSE)
        dirs <- c(dirs, root)
        parts <- .mzs_write_parts(path, root, name, df, generation)
        tables[[name]] <- .mzs_table_entry(
            name, root, nrow(df), length(parts), generation, revision = rev,
            supersedes = if (rev > 1L) rev - 1L,
            fill_policy = if (name == "abundance") tables[[name]]$fill_policy,
            selection = selection[[name]] %||% tables[[name]]$selection)
        index[[length(index) + 1L]] <- .mzs_index_table(name, names(df))
    }

    m$results$tables <- tables
    if (length(sources))
        m$sources <- .mzs_merge_sources(m$sources, sources)
    m$provenance <- c(m$provenance %||% list(), list(.mzs_activity(
        id = act, generation = generation,
        started = activity$started %||% .mzs_now(), action = activity$action,
        agent = .mzs_agent(activity$fn), inputs = activity$inputs %||% list(),
        outputs = list(runs = character(),
                       tables = c(names(append), names(supersede))),
        parameters = activity$parameters)))
    cm <- .mzs_index_merge(.mzs_index_read(path, m), index)
    rel <- sprintf("results/column_map.g%06d.json", generation)
    .mzs_write_json(cm, file.path(path, rel))
    written <- c(written, rel)
    m$results$column_mapping_path <- rel
    m <- .mzs_commit_manifest(path, m, expect)
    written <- character()
    dirs <- character()
    invisible(m)
}
