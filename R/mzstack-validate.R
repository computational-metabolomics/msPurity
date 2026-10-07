# Validation of mzStack datasets.
#
# - a minimal JSON-Schema checker for the draft 2020-12 keywords the vendored
#   schemas use, avoiding a further dependency;
# - `validateMzstack()`, which reports every failure, at MUST or SHOULD level.

.mzs_schema_cache <- new.env(parent = emptyenv())

#' A vendored schema, parsed once per session.
#'
#' @noRd
.mzs_schema <- function(name) {
    if (is.null(.mzs_schema_cache[[name]])) {
        fl <- system.file("mzstack", "schemas", name, package = "msPurity")
        if (!nzchar(fl))
            stop("Schema '", name, "' is not installed with msPurity.",
                 call. = FALSE)
        .mzs_schema_cache[[name]] <- jsonlite::fromJSON(
            fl, simplifyVector = FALSE)
    }
    .mzs_schema_cache[[name]]
}

#' Parse JSON text, or turn an R value into what parsing its JSON would give.
#'
#' @noRd
.mzs_as_json_value <- function(x) {
    txt <- if (inherits(x, "json") || (is.character(x) && length(x) == 1L &&
                                       grepl("^\\s*[[{]", x)))
        x else .mzs_to_json(x, pretty = FALSE)
    jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

#' Check a JSON value against a vendored schema.
#'
#' @param x an R value (converted through its JSON form) or JSON text.
#'
#' @return character vector of failures, each prefixed by the JSON pointer of
#'     the offending value; empty when valid.
#'
#' @noRd
.mzs_schema_check <- function(x, schema) {
    root <- if (is.character(schema)) .mzs_schema(schema) else schema
    .mzs_jsv(.mzs_as_json_value(x), root, root, "")
}

#' @noRd
.mzs_json_type <- function(v) {
    if (is.null(v)) return("null")
    if (is.list(v)) return(if (is.null(names(v))) "array" else "object")
    if (is.logical(v) && length(v) == 1L) return("boolean")
    if (is.character(v) && length(v) == 1L) return("string")
    if (is.numeric(v) && length(v) == 1L)
        return(if (is.finite(v) && v == round(v)) "integer" else "number")
    "unknown"
}

#' @noRd
.mzs_type_ok <- function(v, t) {
    vt <- .mzs_json_type(v)
    vt == t || (t == "number" && vt == "integer")
}

#' @noRd
.mzs_resolve_ref <- function(root, ref) {
    if (!startsWith(ref, "#"))
        stop("Only local $ref is supported, got '", ref, "'.", call. = FALSE)
    parts <- strsplit(sub("^#/?", "", ref), "/", fixed = TRUE)[[1L]]
    node <- root
    for (p in parts[nzchar(parts)]) {
        p <- gsub("~1", "/", gsub("~0", "~", p, fixed = TRUE), fixed = TRUE)
        node <- node[[p]]
    }
    node
}

#' @noRd
.mzs_json_equal <- function(a, b) {
    identical(jsonlite::toJSON(a, auto_unbox = TRUE, null = "null",
                               digits = NA),
              jsonlite::toJSON(b, auto_unbox = TRUE, null = "null",
                               digits = NA))
}

#' Validate `v` against schema node `s`.
#'
#' @noRd
.mzs_jsv <- function(v, s, root, ptr) {
    if (isTRUE(s) || (is.list(s) && !length(s)))
        return(character())
    if (isFALSE(s))
        return(paste0(ptr, ": no value is allowed here"))
    err <- character()
    add <- function(...) err <<- c(err, paste0(ptr, ": ", ...))
    sub_ok <- function(node) !length(.mzs_jsv(v, node, root, ptr))

    if (!is.null(s[["$ref"]]))
        err <- c(err, .mzs_jsv(v, .mzs_resolve_ref(root, s[["$ref"]]), root,
                               ptr))
    if (!is.null(s$type)) {
        ts <- unlist(s$type)
        if (!any(vapply(ts, function(t) .mzs_type_ok(v, t), logical(1)))) {
            add("expected ", paste(ts, collapse = " or "), ", got ",
                .mzs_json_type(v))
            return(err)
        }
    }
    if (!is.null(s$const) || "const" %in% names(s))
        if (!.mzs_json_equal(v, s$const))
            add("must equal ", .mzs_to_json(s$const, pretty = FALSE))
    if (!is.null(s$enum) &&
        !any(vapply(s$enum, function(e) .mzs_json_equal(v, e), logical(1))))
        add("must be one of ", .mzs_to_json(s$enum, pretty = FALSE))

    vt <- .mzs_json_type(v)
    if (vt == "string") {
        if (!is.null(s$minLength) && nchar(v) < s$minLength)
            add("shorter than ", s$minLength)
        if (!is.null(s$pattern) && !grepl(s$pattern, v, perl = TRUE))
            add("'", v, "' does not match ", s$pattern)
    }
    if (vt %in% c("integer", "number") && !is.null(s$minimum) &&
        v < s$minimum)
        add(v, " is below the minimum ", s$minimum)
    if (vt == "array") {
        if (!is.null(s$minItems) && length(v) < s$minItems)
            add("needs at least ", s$minItems, " item(s)")
        if (isTRUE(s$uniqueItems)) {
            enc <- vapply(v, function(e) as.character(jsonlite::toJSON(
                e, auto_unbox = TRUE, null = "null", digits = NA)),
                character(1))
            if (anyDuplicated(enc))
                add("items must be unique")
        }
        if (!is.null(s$items))
            for (i in seq_along(v))
                err <- c(err, .mzs_jsv(v[[i]], s$items, root,
                                       paste0(ptr, "/", i - 1L)))
        if (!is.null(s$contains) &&
            !any(vapply(v, function(e)
                !length(.mzs_jsv(e, s$contains, root, ptr)), logical(1))))
            add("no item matches the required 'contains' schema")
    }
    if (vt == "object") {
        for (r in unlist(s$required))
            if (!r %in% names(v))
                add("missing required key '", r, "'")
        props <- s$properties %||% list()
        for (k in names(v)) {
            kp <- paste0(ptr, "/", k)
            if (!is.null(s$propertyNames))
                err <- c(err, .mzs_jsv(k, s$propertyNames, root, kp))
            if (k %in% names(props)) {
                err <- c(err, .mzs_jsv(v[[k]], props[[k]], root, kp))
            } else if (!is.null(s$additionalProperties)) {
                err <- c(err, .mzs_jsv(v[[k]], s$additionalProperties, root,
                                       kp))
            }
        }
    }
    for (node in s$allOf %||% list())
        err <- c(err, .mzs_jsv(v, node, root, ptr))
    if (!is.null(s$anyOf) && !any(vapply(s$anyOf, sub_ok, logical(1))))
        add("matches none of the 'anyOf' alternatives")
    if (!is.null(s$oneOf)) {
        n <- sum(vapply(s$oneOf, sub_ok, logical(1)))
        if (n != 1L)
            add("matches ", n, " of the 'oneOf' alternatives; exactly one ",
                "is required")
    }
    if (!is.null(s$not) && sub_ok(s$not))
        add("matches a schema it must not match")
    if (!is.null(s[["if"]])) {
        if (sub_ok(s[["if"]])) {
            if (!is.null(s[["then"]]))
                err <- c(err, .mzs_jsv(v, s[["then"]], root, ptr))
        } else if (!is.null(s[["else"]])) {
            err <- c(err, .mzs_jsv(v, s[["else"]], root, ptr))
        }
    }
    unique(err)
}

# ---------------------------------------------------------------------------
# validateMzstack()
# ---------------------------------------------------------------------------

#' A validation finding.
#'
#' @noRd
.mzs_finding <- function(level, requirement, table, message) {
    data.frame(level = level, requirement = requirement,
               table = table %||% NA_character_, message = message,
               stringsAsFactors = FALSE)
}

#' @noRd
.mzs_findings <- function(...) {
    parts <- Filter(function(d) !is.null(d) && nrow(d), list(...))
    if (!length(parts))
        return(.mzs_finding(character(), character(), character(),
                            character()))
    do.call(rbind, parts)
}

#' Validate an mzStack dataset written by msPurity
#'
#' @description
#'
#' Checks an mzStack dataset against the validation rules of the mzStack
#' specification: the manifest, runs and spectra files of every dataset, and
#' the tables, references, keys, sort orders, coverage, loss ledger, naming
#' and CV bindings of a results or library dataset. Every failure is
#' reported, not only the first.
#'
#' Referenced sources need not be resolvable; an unresolvable source is
#' reported at SHOULD level.
#'
#' The append-only rule for results tables cannot be checked: a rewritten
#' Parquet part is indistinguishable from a new one.
#'
#' Requires the suggested packages arrow and jsonlite.
#'
#' @param path `character(1)`, the dataset directory.
#'
#' @param sourcePaths named `character`, source key -> location, overriding
#'     the paths the manifest records.
#'
#' @return A `data.frame` with columns `level` (`"MUST"`, `"SHOULD"` or
#'     `"INFO"`), `requirement` (the identifier in the mzStack
#'     specification), `table` and `message`; zero rows for a dataset with no
#'     findings.
#'
#' @examples
#' if (requireNamespace("arrow", quietly = TRUE) &&
#'     requireNamespace("jsonlite", quietly = TRUE)) {
#'     db <- system.file("extdata", "tests", "db",
#'                       "createDatabase_example.sqlite", package = "msPurity")
#'     out <- file.path(tempdir(), "validate-example.mzstack")
#'     convertSqliteToMzstack(db, out, overwrite = TRUE)
#'     validateMzstack(out)
#' }
#' @export
validateMzstack <- function(path, sourcePaths = NULL) {
    .mzs_require("validateMzstack()")
    if (!.mzs_is_dataset(path))
        return(.mzs_finding("MUST", "D-001", NA_character_,
                            paste0("'", path, "' has no ", .MZS_MANIFEST,
                                   ".")))
    path <- normalizePath(path)
    raw <- tryCatch(jsonlite::fromJSON(.mzs_manifest_path(path),
                                       simplifyVector = FALSE),
                    error = function(e) e)
    if (inherits(raw, "error"))
        return(.mzs_finding("MUST", "D-005", NA_character_,
                            paste("The manifest is not JSON:",
                                  conditionMessage(raw))))
    m <- tryCatch(.mzs_read_manifest(path), mzstack_format = function(e) e)
    if (inherits(m, "error"))
        return(.mzs_finding("MUST", "G-009", NA_character_,
                            conditionMessage(m)))
    out <- .mzs_findings(
        .mzs_check_manifest_schema(raw),
        .mzs_check_runs(path, m))
    if (identical(m$role, "results") || identical(m$role, "library") ||
        !is.null(m$results))
        out <- .mzs_findings(out, .mzs_validate_results(path, m, sourcePaths))
    rownames(out) <- NULL
    out
}

#' @noRd
.mzs_check_manifest_schema <- function(raw) {
    err <- .mzs_schema_check(raw, "manifest.schema.json")
    if (!length(err))
        return(NULL)
    .mzs_finding("MUST", "D-005", NA_character_,
                 paste("manifest.schema.json:", err))
}

#' Checks of the run entries and their files.
#'
#' @noRd
.mzs_check_runs <- function(path, m) {
    f <- list()
    add <- function(level, req, msg, table = NA_character_)
        f[[length(f) + 1L]] <<- .mzs_finding(level, req, table, msg)
    rf <- .mzs_runs_frame(m)
    if (!nrow(rf))
        return(NULL)
    kinds <- unique(rf$kind[!is.na(rf$kind)])
    if (length(kinds) > 1L)
        add("MUST", "G-007", paste("The dataset mixes run kinds:",
                                   paste(kinds, collapse = ", ")))
    if (anyNA(rf$run_id) || anyNA(rf$kind))
        add("MUST", "D-013", "Every run entry must carry run_id and kind.")
    bad <- rf$run_id[!is.na(rf$run_id) &
                     !grepl("^[A-Za-z0-9._-]{1,128}$", rf$run_id)]
    if (length(bad))
        add("MUST", "D-014", paste("Invalid run_id(s):",
                                   paste(bad, collapse = ", ")))
    if (anyDuplicated(rf$run_id))
        add("MUST", "D-014", paste("Duplicate run_id(s):",
                                   paste(unique(rf$run_id[duplicated(
                                       rf$run_id)]), collapse = ", ")))
    for (r in m$runs) {
        lay <- r$signal$layout %||% NA_character_
        want <- if (identical(r$kind, "native")) "list" else "point"
        if (!identical(lay, want))
            add("MUST", "D-015", paste0("Run '", r$run_id, "' of kind ",
                                        r$kind, " declares layout '", lay,
                                        "'; expected '", want, "'."))
        for (pt in names(r$projections %||% list()))
            if (!identical(as.integer(r$projections[[pt]]$generation),
                           as.integer(r$ingested_at)))
                add("MUST", "D-057", paste0("Projection '", pt, "' of run '",
                                            r$run_id, "' is stale."))
    }
    expect <- cumsum(c(1L, utils::head(rf$n_spectra, -1L)))
    if (!isTRUE(all(rf$uid_base == expect)))
        add("MUST", "D-019", paste("uid_base values are not contiguous",
                                   "blocks in manifest order:",
                                   paste(rf$uid_base, collapse = ", ")))
    for (i in seq_len(nrow(rf))) {
        r <- m$runs[[i]]
        if (!identical(r$kind, "native"))
            next
        res <- tryCatch(.mzs_check_native_run(path, r),
                        error = function(e) .mzs_finding(
                            "MUST", "D-048", NA_character_,
                            paste0("Run '", r$run_id, "' could not be read: ",
                                   conditionMessage(e))))
        f[[length(f) + 1L]] <- res
    }
    do.call(.mzs_findings, f)
}

#' Checks of one native run's files: row count, identities, list lengths.
#'
#' @noRd
.mzs_check_native_run <- function(path, r) {
    f <- list()
    add <- function(level, req, msg)
        f[[length(f) + 1L]] <<- .mzs_finding(level, req, NA_character_, msg)
    dir <- .mzs_run_dir(path, r)
    ds <- arrow::open_dataset(dir, format = "parquet", unify_schemas = TRUE)
    cols <- names(ds)
    miss <- setdiff(c("spectrum_id_", "spectrum_index", "mz", "intensity"),
                    cols)
    if (length(miss)) {
        add("MUST", "D-050", paste0("Run '", r$run_id, "' lacks column(s) ",
                                    paste(miss, collapse = ", "), "."))
        return(do.call(.mzs_findings, f))
    }
    ## List lengths are computed in Arrow: peaks never become R objects.
    tab <- arrow::as_arrow_table(ds)
    if (tab$num_rows != as.integer(r$n_spectra))
        add("MUST", "D-063", paste0("Run '", r$run_id, "' holds ",
                                    tab$num_rows, " spectra; the manifest ",
                                    "declares ", r$n_spectra, "."))
    if (!isTRUE(all(as.vector(tab$spectrum_id_) ==
                    as.integer(r$uid_base) + as.vector(tab$spectrum_index))))
        add("MUST", "D-018", paste0("Run '", r$run_id, "': spectrum_id_ is ",
                                    "not uid_base + spectrum_index."))
    len <- function(c) as.vector(arrow::call_function("list_value_length",
                                                      tab[[c]]))
    lens <- len("mz")
    li <- len("intensity")
    lens[is.na(lens)] <- 0L
    li[is.na(li)] <- 0L
    if (!identical(lens, li))
        add("MUST", "D-051", paste0("Run '", r$run_id, "': mz and intensity ",
                                    "differ in length in some row."))
    lists <- names(tab)[vapply(names(tab), function(c)
        inherits(tab[[c]]$type, "ListType"), logical(1))]
    for (c in setdiff(lists, c("mz", "intensity"))) {
        l <- len(c)
        if (any(!is.na(l) & l != lens))
            add("MUST", "R-059", paste0("Run '", r$run_id, "': peak ",
                                        "annotation '", c, "' differs in ",
                                        "length from mz."))
    }
    sch <- arrow::open_dataset(dir, format = "parquet")$schema
    child <- vapply(c("mz", "intensity"), function(c) {
        t <- sch$GetFieldByName(c)$type
        if (inherits(t, "ListType")) t$value_field$name else NA_character_
    }, character(1))
    if (any(child != "item", na.rm = TRUE))
        add("INFO", "D-049", paste0("Run '", r$run_id, "': list child field ",
                                    "is '", child[1L], "', not 'item'. ",
                                    "arrow's R writer cannot set it."))
    do.call(.mzs_findings, f)
}

#' Results-layer checks.
#'
#' @noRd
.mzs_validate_results <- function(path, m, sourcePaths) {
    f <- list()
    add <- function(level, req, msg, table = NA_character_)
        f[[length(f) + 1L]] <<- .mzs_finding(level, req, table, msg)
    if (identical(m$role, "results")) {
        if (is.null(m$results))
            add("MUST", "D-010", "role 'results' without a results block.")
        maj <- .mzs_semver_major(m$results$version)
        if (is.na(maj) || maj != .MZS_RESULTS_MAJOR)
            add("MUST", "R-034", paste0("results.version ",
                                        m$results$version %||% "<missing>",
                                        " is not ", .MZS_RESULTS_MAJOR,
                                        ".x."))
    }
    keys <- vapply(m$sources %||% list(), function(s)
        as.character(s$key %||% NA_character_), character(1))
    if (anyDuplicated(keys))
        add("MUST", "R-024", paste("Duplicate source key(s):",
                                   paste(unique(keys[duplicated(keys)]),
                                         collapse = ", ")))
    if ("self" %in% keys)
        add("MUST", "R-020", "A source is declared with the reserved key 'self'.")
    for (k in setdiff(unique(keys), c("self", NA))) {
        st <- tryCatch({
            .mzs_open_source(path, m, k, sourcePaths)
            NULL
        }, mzstack_capability = function(e) NULL,
        mzstack_reference = function(e)
            .mzs_finding("SHOULD", "R-100", NA_character_,
                         conditionMessage(e)),
        mzstack_stale = function(e)
            .mzs_finding("MUST", "R-029", NA_character_,
                         conditionMessage(e)))
        f[[length(f) + 1L]] <- st
    }
    cmp <- m$results$column_mapping_path
    if (!is.null(cmp)) {
        if (startsWith(cmp, "index/"))
            add("MUST", "R-062", "The results index is stored under index/.")
        fl <- file.path(path, cmp)
        if (!file.exists(fl)) {
            add("MUST", "R-061", paste0("No results index at '", cmp, "'."))
        } else {
            err <- .mzs_schema_check(
                jsonlite::fromJSON(fl, simplifyVector = FALSE),
                "results_index.schema.json")
            if (length(err))
                add("MUST", "R-061", paste("results_index.schema.json:", err))
        }
    }
    f[[length(f) + 1L]] <- .mzs_validate_tables(path, m)
    do.call(.mzs_findings, f)
}

#' Table-level checks.
#'
#' An error inside one check group is reported as a finding and the others
#' still run.
#'
#' @noRd
.mzs_validate_tables <- function(path, m) {
    tb <- m$results$tables %||% list()
    if (!length(tb))
        return(NULL)
    tabs <- list()
    read <- function(n) {
        if (is.null(tabs[[n]]))
            tabs[[n]] <<- .mzs_read_table(path, m, n)
        tabs[[n]]
    }
    have <- function(n) !is.null(tb[[n]]) && !isTRUE(tb[[n]]$omitted)
    checks <- list(
        structure = .mzs_v_structure, references = .mzs_v_references,
        integrity = .mzs_v_integrity, derived = .mzs_v_derived,
        sort_orders = .mzs_v_sort_orders, truncation = .mzs_v_truncation,
        ledger = .mzs_v_ledger, naming = .mzs_v_naming,
        cv_binding = .mzs_v_cv_binding)
    out <- lapply(names(checks), function(nm)
        tryCatch(checks[[nm]](path, m, read, have),
                 error = function(e) .mzs_finding(
                     "MUST", "R-099", NA_character_,
                     paste0("The ", nm, " checks could not run: ",
                            conditionMessage(e)))))
    do.call(.mzs_findings, out)
}

#' Every declared table present, readable, of its declared row count, not
#' under index/, and declared in the results index.
#'
#' @noRd
.mzs_v_structure <- function(path, m, read, have) {
    f <- list()
    add <- function(req, msg, table = NA_character_, level = "MUST")
        f[[length(f) + 1L]] <<- .mzs_finding(level, req, table, msg)
    idx <- .mzs_index_read(path, m)
    declared <- vapply(idx$files %||% list(), function(x)
        paste(x$data_kind, x$name), "")
    for (n in names(m$results$tables)) {
        e <- m$results$tables[[n]]
        if (isTRUE(e$omitted))
            next
        if (startsWith(e$path %||% "", "index/"))
            add("B-038", "A results table is stored under index/.", n)
        rows <- tryCatch(nrow(read(n)), error = function(e) NA)
        if (is.na(rows)) {
            add("R-037", "The table could not be read.", n)
        } else if (!identical(as.integer(rows), as.integer(e$rows))) {
            add("R-037", paste0("Holds ", rows, " rows; the manifest ",
                                "declares ", e$rows, "."), n)
        }
        if (!paste("table", n) %in% declared)
            add("R-061", "The results index does not declare this table.", n)
    }
    if (identical(m$role, "results")) {
        for (r in m$runs) {
            sch <- arrow::open_dataset(.mzs_run_dir(path, r),
                                       format = "parquet")$schema
            lists <- setdiff(names(sch)[vapply(names(sch), function(c)
                inherits(sch$GetFieldByName(c)$type, "ListType"),
                logical(1))], c("mz", "intensity"))
            if (length(lists) && !paste("peak_annotation", r$run_id) %in%
                declared)
                add("R-058", paste0("Run '", r$run_id, "' carries ",
                                    "peak-annotation column(s) the results ",
                                    "index does not declare."))
            if (is.null(m$results$runs[[r$run_id]]))
                add("R-010", paste0("Run '", r$run_id, "' has no ",
                                    "results.runs entry declaring its ",
                                    "provenance scope."))
        }
    }
    do.call(.mzs_findings, f)
}

#' Spectrum references: declared sources, complete components, self
#' references inside a declared run.
#'
#' @noRd
.mzs_v_references <- function(path, m, read, have) {
    f <- list()
    add <- function(req, msg, table)
        f[[length(f) + 1L]] <<- .mzs_finding("MUST", req, table, msg)
    keys <- c("self", vapply(m$sources %||% list(), function(s)
        as.character(s$key), ""))
    rf <- .mzs_runs_frame(m)
    for (n in names(m$results$tables)) {
        if (!have(n) || !n %in% names(.MZS_TABLES))
            next
        df <- read(n)
        for (p in .mzs_ref_prefixes(n)) {
            src <- df[[paste0(p, "_source")]]
            run <- df[[paste0(p, "_run_id")]]
            id <- df[[paste0(p, "_spectrum_id_")]]
            bad <- unique(src[!src %in% keys])
            if (length(bad))
                add("R-018", paste0(p, "_source names undeclared source(s): ",
                                    paste(bad, collapse = ", ")), n)
            miss <- !is.na(id) & (is.na(src) | is.na(run))
            if (any(miss))
                add("R-021", paste0(sum(miss), " ", p, "_spectrum_id_ ",
                                    "value(s) lack their source or run."), n)
            self <- which(src %in% "self" & !is.na(id))
            if (length(self)) {
                k <- match(run[self], rf$run_id)
                out <- is.na(k) | id[self] < rf$uid_base[k] |
                    id[self] >= rf$uid_base[k] + rf$n_spectra[k]
                if (any(out))
                    add("R-022", paste0(sum(out), " self reference(s) of ",
                                        p, " fall outside their run's ",
                                        "spectrum ids, e.g. ",
                                        id[self][out][1L], "."), n)
            }
        }
    }
    do.call(.mzs_findings, f)
}

#' Foreign keys, reported as orphan rows.
#'
#' @noRd
.mzs_v_integrity <- function(path, m, read, have) {
    f <- list()
    orphan <- function(table, col, keys, target) {
        if (!have(table))
            return()
        v <- read(table)[[col]]
        bad <- unique(v[!is.na(v) & !v %in% keys])
        if (length(bad))
            f[[length(f) + 1L]] <<- .mzs_finding(
                "MUST", "R-099", table,
                paste0(col, " names ", length(bad), " row(s) absent from ",
                       target, ", e.g. ", bad[1L], "."))
    }
    key <- function(t, k) if (have(t)) read(t)[[k]] else numeric()
    fid <- key("feature", "feature_id_")
    pid <- key("chromatographic_peak", "chromatographic_peak_id_")
    orphan("chromatographic_peak_feature", "feature_id_", fid, "feature")
    orphan("chromatographic_peak_feature", "chromatographic_peak_id_", pid,
           "chromatographic_peak")
    orphan("abundance", "feature_id_", fid, "feature")
    orphan("abundance", "assay_id_", key("assay", "assay_id_"), "assay")
    orphan("assay", "sample_id_", key("sample", "sample_id_"), "sample")
    orphan("spectrum_feature", "feature_id_", fid, "feature")
    orphan("spectrum_feature", "chromatographic_peak_id_", pid,
           "chromatographic_peak")
    if (have("spectrum_feature")) {
        sf <- read("spectrum_feature")
        cp <- sf$link_mode == "chromatographic_peak"
        if (any(cp != !is.na(sf$chromatographic_peak_id_)))
            f[[length(f) + 1L]] <- .mzs_finding(
                "MUST", "R-099", "spectrum_feature",
                paste("chromatographic_peak_id_ must be set exactly where",
                      "link_mode is chromatographic_peak."))
        bad <- setdiff(unique(sf$link_mode), .MZS_LINK_MODES)
        if (length(bad))
            f[[length(f) + 1L]] <- .mzs_finding(
                "MUST", "R-044", "spectrum_feature",
                paste("Unknown link mode(s):", paste(bad, collapse = ", ")))
        fw <- sf$link_mode == "feature_width" & sf$run_scoped
        if (any(fw))
            f[[length(f) + 1L]] <- .mzs_finding(
                "MUST", "R-045", "spectrum_feature",
                "feature_width associations must not be run-scoped.")
    }
    orphan("evidence", "feature_id_", fid, "feature")
    orphan("evidence", "compound_id_", key("compound", "compound_id_"),
           "compound")
    orphan("evidence", "software_id_", key("software", "software_id_"),
           "software")
    orphan("evidence_score", "evidence_id_", key("evidence", "evidence_id_"),
           "evidence")
    orphan("compound_xref", "compound_id_", key("compound", "compound_id_"),
           "compound")
    orphan("compound_synonym", "compound_id_", key("compound", "compound_id_"),
           "compound")
    orphan("x_mspurity_feature_annotation", "feature_id_", fid, "feature")
    orphan("x_mspurity_combined_annotation", "feature_id_", fid, "feature")
    orphan("x_mspurity_combined_annotation", "compound_id_",
           key("compound", "compound_id_"), "compound")
    orphan("x_mspurity_scan_peak", "scan_annotation_id_",
           key("x_mspurity_scan", "scan_annotation_id_"), "x_mspurity_scan")
    acts <- vapply(m$provenance %||% list(), function(a) as.character(a$id),
                   "")
    for (n in names(m$results$tables))
        orphan(n, "activity", acts, "provenance")
    if (have("source_identifier")) {
        si <- read("source_identifier")
        for (t in unique(si$target_table)) {
            k <- si$target_key[si$target_table == t]
            ok <- if (t == "spectra")
                k >= 1 & k <= sum(.mzs_runs_frame(m)$n_spectra)
            else if (have(t) && !is.na(.mzs_table_schema(t)$key))
                k %in% read(t)[[.mzs_table_schema(t)$key]]
            else rep(FALSE, length(k))
            if (any(!ok))
                f[[length(f) + 1L]] <- .mzs_finding(
                    "MUST", "R-099", "source_identifier",
                    paste0(sum(!ok), " target_key(s) absent from ", t, "."))
        }
    }
    do.call(.mzs_findings, f)
}

#' Derived spectra: every one has a merge_member row, members_complete is
#' set, and aggregation methods are in separate runs.
#'
#' @noRd
.mzs_v_derived <- function(path, m, read, have) {
    if (!identical(m$role, "results") || !length(m$runs))
        return(NULL)
    f <- list()
    rr <- m$results$runs %||% list()
    methods <- vapply(m$runs, function(r)
        as.character(rr[[r$run_id]]$method %||% NA), "")
    scope <- vapply(m$runs, function(r)
        as.character(rr[[r$run_id]]$scope %||% NA), "")
    ds <- methods[scope %in% "dataset"]
    if (anyDuplicated(ds[!is.na(ds)]))
        f[[length(f) + 1L]] <- .mzs_finding(
            "MUST", "R-011", NA_character_,
            "Two dataset-scoped runs hold the same aggregation method.")
    n <- sum(.mzs_runs_frame(m)$n_spectra)
    if (!have("merge_member")) {
        if (n)
            f[[length(f) + 1L]] <- .mzs_finding(
                "MUST", "R-047", NA_character_,
                "Derived spectra are present but merge_member is not.")
        return(do.call(.mzs_findings, f))
    }
    mm <- read("merge_member")
    missing <- setdiff(seq_len(n), mm$merged_spectrum_id_)
    if (length(missing))
        f[[length(f) + 1L]] <- .mzs_finding(
            "MUST", "R-047", "merge_member",
            paste0(length(missing), " derived spectra have no members, e.g. ",
                   missing[1L], "."))
    if (anyNA(mm$members_complete))
        f[[length(f) + 1L]] <- .mzs_finding(
            "MUST", "R-099", "merge_member", "members_complete is null.")
    do.call(.mzs_findings, f)
}

#' Declared sort orders hold within every part, and no sort column holds
#' NaN.
#'
#' @noRd
.mzs_v_sort_orders <- function(path, m, read, have) {
    f <- list()
    for (n in names(m$results$tables)) {
        e <- m$results$tables[[n]]
        if (!have(n))
            next
        by <- unlist(e$sorted_by)
        parts <- list.files(file.path(path, e$path), pattern = "\\.parquet$",
                            recursive = TRUE, full.names = TRUE)
        for (p in parts) {
            df <- as.data.frame(arrow::read_parquet(p))
            cols <- intersect(by, names(df))
            if (!length(cols) || nrow(df) < 2L)
                next
            for (c in cols)
                if (is.factor(df[[c]])) df[[c]] <- as.character(df[[c]])
            if (any(vapply(df[cols], function(v) is.double(v) &&
                           any(is.nan(v)), logical(1))))
                f[[length(f) + 1L]] <- .mzs_finding(
                    "MUST", "R-099", n,
                    paste0("A sort column holds NaN in ", basename(p), "."))
            o <- do.call(order, c(unname(as.list(df[cols])),
                                  list(method = "radix", na.last = TRUE)))
            if (!identical(o, seq_len(nrow(df))))
                f[[length(f) + 1L]] <- .mzs_finding(
                    "MUST", "R-038", n,
                    paste0("Part ", basename(p), " is not in its declared ",
                           "order (", paste(by, collapse = ", "), ")."))
            else if (anyDuplicated(df[cols]) && identical(cols, by))
                f[[length(f) + 1L]] <- .mzs_finding(
                    "MUST", "R-038", n,
                    paste0("The declared order admits ties in ",
                           basename(p), "."))
        }
    }
    do.call(.mzs_findings, f)
}

#' Truncation and coverage: a cut records its boundary, coverage counts
#' equal each query's evidence rows, and evidence has coverage or a scope
#' predicate.
#'
#' @noRd
.mzs_v_truncation <- function(path, m, read, have) {
    if (!have("evidence"))
        return(NULL)
    f <- list()
    e <- m$results$tables$evidence
    if (is.null(e$selection))
        f[[length(f) + 1L]] <- .mzs_finding(
            "MUST", "R-075", "evidence", "No selection is declared.")
    if (!have("coverage") && is.null(e$scope_predicate))
        return(do.call(.mzs_findings, c(f, list(.mzs_finding(
            "MUST", "R-074", "evidence",
            "Neither coverage nor a scope predicate accompanies evidence.")))))
    if (!have("coverage"))
        return(do.call(.mzs_findings, f))
    cv <- read("coverage")
    ev <- read("evidence")
    if (any(cv$truncated & is.na(cv$boundary_score)))
        f[[length(f) + 1L]] <- .mzs_finding(
            "MUST", "R-076", "coverage",
            "A truncated query records no boundary score.")
    key <- function(d) paste(d$query_source, d$query_run_id,
                             d$query_spectrum_id_, d$query_native_id,
                             d$activity)
    n <- table(factor(key(ev), levels = key(cv)))
    bad <- as.integer(n) != cv$hits_retained
    if (any(bad))
        f[[length(f) + 1L]] <- .mzs_finding(
            "MUST", "B-044", "coverage",
            paste0(sum(bad), " coverage row(s) disagree with the evidence ",
                   "rows of their query."))
    stray <- setdiff(unique(key(ev)), key(cv))
    if (length(stray))
        f[[length(f) + 1L]] <- .mzs_finding(
            "MUST", "B-044", "evidence",
            paste0(length(stray), " queries have evidence but no coverage."))
    do.call(.mzs_findings, f)
}

#' Converter tables: known coverage states and dispositions, resolvable
#' ledger targets, and a ledger row for every wholly-null column.
#'
#' @noRd
.mzs_v_ledger <- function(path, m, read, have) {
    if (!have("conversion"))
        return(NULL)
    f <- list()
    add <- function(req, msg, table)
        f[[length(f) + 1L]] <<- .mzs_finding("MUST", req, table, msg)
    cv <- read("conversion")
    for (cm in cv$coverage_manifest) {
        st <- unlist(jsonlite::fromJSON(cm, simplifyVector = FALSE))
        bad <- setdiff(st, .MZS_COVERAGE_STATES)
        if (length(bad))
            add("R-080", paste("Unknown coverage state(s):",
                               paste(bad, collapse = ", ")), "conversion")
    }
    if (!have("loss_ledger"))
        return(do.call(.mzs_findings, c(f, list(.mzs_finding(
            "MUST", "R-082", NA_character_,
            "A converted dataset carries no loss_ledger.")))))
    ll <- read("loss_ledger")
    bad <- setdiff(ll$disposition, .MZS_DISPOSITIONS)
    if (length(bad))
        add("R-083", paste("Unknown disposition(s):",
                           paste(bad, collapse = ", ")), "loss_ledger")
    if (any(ll$disposition == "dropped_unmappable"))
        f[[length(f) + 1L]] <- .mzs_finding(
            "SHOULD", "R-085", "loss_ledger",
            "dropped_unmappable requires the lossy_conversion profile.")
    resolves <- function(ref) {
        t <- sub("\\..*$", "", ref)
        c <- if (grepl(".", ref, fixed = TRUE)) sub("^[^.]*\\.", "", ref)
             else NA
        if (!have(t))
            return(FALSE)
        is.na(c) || c %in% names(read(t))
    }
    refs <- unique(ll$target_ref[!is.na(ll$target_ref)])
    bad <- refs[!vapply(refs, resolves, logical(1))]
    if (length(bad))
        add("R-099", paste("Ledger target(s) that do not resolve:",
                           paste(bad, collapse = ", ")), "loss_ledger")
    for (n in names(m$results$tables)) {
        if (!have(n) || n == "loss_ledger" || !n %in% names(.MZS_TABLES))
            next
        df <- read(n)
        if (!nrow(df))
            next
        for (c in setdiff(names(df), "activity")) {
            v <- df[[c]]
            if (is.list(v) || !all(is.na(v) & !is.nan(v)))
                next
            if (!paste0(n, ".", c) %in% ll$target_ref)
                add("R-084", paste0("Column ", c, " is wholly null and has ",
                                    "no ledger row."), n)
        }
    }
    do.call(.mzs_findings, f)
}

#' Column naming, InChIKeys, ranks and score kinds.
#'
#' @noRd
.mzs_v_naming <- function(path, m, read, have) {
    f <- list()
    add <- function(req, msg, table)
        f[[length(f) + 1L]] <<- .mzs_finding("MUST", req, table, msg)
    units <- "_(ppm|da|dalton|seconds|second|sec|minutes|minute|min|ev)$"
    for (n in names(m$results$tables)) {
        if (!have(n))
            next
        df <- read(n)
        reg <- if (n %in% names(.MZS_TABLES)) .mzs_table_schema(n)$columns
        types <- vapply(names(df), function(c) {
            i <- match(c, reg$name)
            if (!is.na(i)) reg$type[i]
            else if (is.integer(df[[c]]) || is.double(df[[c]])) "i64" else ""
        }, "")
        bad <- .mzs_bad_names(names(df), types)
        if (length(bad))
            add("G-016", paste("Column name(s) break the naming rule:",
                               paste(bad, collapse = ", ")), n)
        own <- setdiff(names(df), reg$name)
        bad <- own[grepl(units, own, ignore.case = TRUE)]
        if (length(bad))
            add("G-014", paste("Column name(s) encode a unit:",
                               paste(bad, collapse = ", ")), n)
    }
    if (have("compound")) {
        cmp <- read("compound")
        bad <- !is.na(cmp$inchikey) &
            !grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", cmp$inchikey)
        if (any(bad))
            add("R-049", paste(sum(bad), "inchikey value(s) are not full",
                               "InChIKeys."), "compound")
        both <- !is.na(cmp$inchikey) & !is.na(cmp$inchikey_block1)
        if (any(substr(cmp$inchikey[both], 1, 14) != cmp$inchikey_block1[both]))
            add("R-049", "inchikey_block1 disagrees with inchikey.",
                "compound")
    }
    if (have("evidence") && have("evidence_score")) {
        ev <- read("evidence")
        es <- read("evidence_score")
        dpc <- es[es$score_term == .MZS_SCORE_TERMS[["dpc"]], ]
        s <- dpc$score_value[match(ev$evidence_id_, dpc$evidence_id_)]
        g <- paste(ev$evidence_input_id, ev$activity)
        ok <- !is.na(s) & !is.na(ev$rank)
        if (any(ok)) {
            want <- stats::ave(-s[ok], g[ok], FUN = function(x)
                rank(x, ties.method = "min"))
            if (any(want != ev$rank[ok]))
                add("R-051", paste("Ranks are not increasing from 1 with",
                                   "tied scores sharing a rank."), "evidence")
        }
        bad <- setdiff(es$score_kind, .MZS_SCORE_KINDS)
        if (length(bad))
            add("R-055", paste("Unknown score kind(s):",
                               paste(bad, collapse = ", ")), "evidence_score")
        if (any(es$score_kind == "rescaled" & is.na(es$rescaling_scope)))
            add("R-055", "A rescaled score has no rescaling scope.",
                "evidence_score")
    }
    if (have("x_mspurity_feature_annotation")) {
        fa <- read("x_mspurity_feature_annotation")
        if (any(fa$score_kind == "rescaled" & is.na(fa$rescaling_scope)))
            add("R-055", "A rescaled score has no rescaling scope.",
                "x_mspurity_feature_annotation")
    }
    do.call(.mzs_findings, f)
}

#' Every bound namespace is in cv_list, and every bound accession is in
#' msPurity's verified term list (inst/mzstack/cv/terms.tsv).
#'
#' @noRd
.mzs_v_cv_binding <- function(path, m, read, have) {
    idx <- .mzs_index_read(path, m)
    if (is.null(idx))
        return(NULL)
    f <- list()
    ids <- vapply(idx$cv_list %||% list(), function(c) c$id, "")
    known <- .mzs_known_terms()
    for (file in idx$files %||% list())
        for (b in file$column_mapping %||% list())
            for (acc in c(b$accession, b$unit)) {
                ns <- sub(":.*$", "", acc)
                if (!ns %in% ids)
                    f[[length(f) + 1L]] <- .mzs_finding(
                        "MUST", "R-063", file$name,
                        paste0(acc, " is in namespace ", ns, ", absent from ",
                               "cv_list."))
                else if (!acc %in% known)
                    f[[length(f) + 1L]] <- .mzs_finding(
                        "MUST", "R-063", file$name,
                        paste0(acc, " is not a term of the declared ",
                               "ontology version (as far as msPurity's ",
                               "term list records)."))
            }
    do.call(.mzs_findings, f)
}
