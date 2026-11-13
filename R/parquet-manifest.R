# The dataset manifest, sources and provenance.
#
# `mzStack.json` is read and written with MsBackendParquet's manifest API
# (readManifest(), writeManifest()), which keeps the results layer's keys
# exactly as parsed and commits atomically. The manifest is the commit point
# of every write: written last, after everything it declares is on disk.
#
# Also here: the JSON writer for msPurity's own documents (the results
# index, JSON table cells), source entries with per-run fingerprints, lazy
# source resolution, and provenance activity records.

.MZS_SPEC_VERSION <- "0.1.0"
.MZS_RESULTS_VERSION <- "1.0.0"
.MZS_RESULTS_MAJOR <- 1L
.MZS_MANIFEST <- "mzStack.json"
.MZS_COLUMN_MAP_PATH <- "results/column_map.json"
.MZS_KEY_RE <- "^[A-Za-z0-9._-]{1,128}$"

#' @noRd
.mzs_semver_major <- function(v) {
    suppressWarnings(as.integer(sub("\\..*$", "", as.character(v)[1L])))
}

#' A JSON array whatever its length: a list is never unboxed.
#'
#' @noRd
.mzs_array <- function(x) as.list(unname(x))

#' An empty JSON object.
#'
#' @noRd
.mzs_object <- function() structure(list(), names = character())

#' Current time as an RFC 3339 UTC timestamp.
#'
#' @noRd
.mzs_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")

#' Evaluate `expr` without disturbing the caller's random number stream.
#'
#' @noRd
.mzs_with_own_seed <- function(expr) {
    had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had)
        old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
        if (had)
            assign(".Random.seed", old, envir = globalenv())
        else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv())
    })
    set.seed(NULL)
    expr
}

#' A random UUIDv4, for a dataset `uid`.
#'
#' @noRd
.mzs_uid <- function() {
    b <- .mzs_with_own_seed(sample.int(256L, 16L, replace = TRUE) - 1L)
    b[7L] <- bitwOr(bitwAnd(b[7L], 0x0FL), 0x40L)
    b[9L] <- bitwOr(bitwAnd(b[9L], 0x3FL), 0x80L)
    h <- sprintf("%02x", b)
    paste(paste(h[1:4], collapse = ""), paste(h[5:6], collapse = ""),
          paste(h[7:8], collapse = ""), paste(h[9:10], collapse = ""),
          paste(h[11:16], collapse = ""), sep = "-")
}

#' A short random suffix for staging directories and part files.
#'
#' @noRd
.mzs_rand <- function(n = 8L) {
    paste(.mzs_with_own_seed(sample(c(letters, 0:9), n, replace = TRUE)),
          collapse = "")
}

# ---------------------------------------------------------------------------
# The manifest and msPurity's own JSON documents.
# ---------------------------------------------------------------------------

#' @noRd
.mzs_manifest_path <- function(path) file.path(path, .MZS_MANIFEST)

#' @noRd
.mzs_is_dataset <- function(path) {
    is.character(path) && length(path) == 1L && dir.exists(path) &&
        file.exists(.mzs_manifest_path(path))
}

#' The shortest decimal rendering of each finite double that reads back
#' exactly; `null` for NA, NaN and infinities.
#'
#' @noRd
.mzs_shortest_double <- function(x) {
    out <- rep("null", length(x))
    ok <- is.finite(x)
    if (any(ok))
        out[ok] <- as.vector(arrow::Array$create(x[ok])$cast(arrow::utf8()))
    out
}

#' Replace every double in a list by its exact JSON text.
#'
#' @noRd
.mzs_json_exact <- function(x) {
    if (is.list(x)) {
        at <- attributes(x)
        x <- lapply(x, .mzs_json_exact)
        attributes(x) <- at
        return(x)
    }
    if (!is.double(x) || inherits(x, "json"))
        return(x)
    s <- .mzs_shortest_double(as.vector(x))
    if (length(s) != 1L)
        s <- paste0("[", paste(s, collapse = ","), "]")
    structure(s, class = "json")
}

#' JSON text for a results index or parameter list.
#'
#' @noRd
.mzs_to_json <- function(x, pretty = TRUE) {
    jsonlite::toJSON(.mzs_json_exact(x), auto_unbox = TRUE, pretty = pretty,
                     null = "null", na = "null", json_verbatim = TRUE)
}

#' Write a JSON document atomically via a temporary file and rename.
#'
#' @noRd
.mzs_write_json <- function(x, fl) {
    dir.create(dirname(fl), recursive = TRUE, showWarnings = FALSE)
    tmp <- paste0(fl, ".tmp-", Sys.getpid(), "-", .mzs_rand())
    writeLines(.mzs_to_json(x), tmp, useBytes = TRUE)
    if (!file.rename(tmp, fl)) {
        unlink(tmp)
        stop("Could not write '", fl, "'.", call. = FALSE)
    }
    invisible(fl)
}

#' Commit a new version of an existing dataset's manifest.
#'
#' `m` was read at generation `expect`. `writeManifest()` serialises the
#' check and the rename with a lock; if another writer holds it or has
#' committed meanwhile, nothing is written and the `mzstack_conflict` error
#' is raised with the class `msPurity_parquet_conflict` added.
#'
#' @return the committed manifest, at generation `expect + 1`.
#'
#' @noRd
.mzs_commit_manifest <- function(path, m, expect) {
    m$generation <- as.integer(expect)
    withCallingHandlers(
        MsBackendParquet::writeManifest(path, m),
        mzstack_conflict = function(e) .mzs_conflict(e))
}

#' Re-signal a manifest conflict with msPurity's class added.
#'
#' @noRd
.mzs_conflict <- function(e) {
    class(e) <- c("msPurity_parquet_conflict", class(e))
    stop(e)
}

# ---------------------------------------------------------------------------
# Runs.
# ---------------------------------------------------------------------------

#' @noRd
.mzs_runs_frame <- function(m) {
    runs <- m$runs %||% list()
    chr <- function(f) vapply(runs, function(r)
        as.character(r[[f]] %||% NA_character_)[1L], character(1))
    int <- function(f) vapply(runs, function(r)
        as.integer(r[[f]] %||% NA_integer_)[1L], integer(1))
    data.frame(run_id = chr("run_id"), kind = chr("kind"),
               path = chr("path"), n_spectra = int("n_spectra"),
               uid_base = int("uid_base"), ingested_at = int("ingested_at"),
               stringsAsFactors = FALSE)
}

#' The directory holding a run's spectrum rows.
#'
#' A native run keeps them under its declared path; an mzpeak run's derived
#' index is under `index/spectra/`.
#'
#' @noRd
.mzs_run_dir <- function(dataset, run) {
    if (identical(run$kind, "native")) {
        p <- run$path %||% file.path("spectra", paste0("run_id=",
                                                         run$run_id))
        if (!.mzs_is_absolute(p))
            p <- file.path(dataset, p)
        return(p)
    }
    file.path(dataset, "index", "spectra", paste0("run_id=", run$run_id))
}

#' @noRd
.mzs_is_absolute <- function(p) grepl("^(/|[A-Za-z]:[/\\\\]|~)", p)

#' Every stored column of a run's spectrum rows, with their stored types, in
#' `spectrum_id_` order.
#'
#' @noRd
.mzs_read_run <- function(dataset, run) {
    ds <- arrow::open_dataset(.mzs_run_dir(dataset, run), format = "parquet",
                              unify_schemas = TRUE)
    df <- as.data.frame(dplyr::collect(ds))
    df <- df[order(df$spectrum_id_), , drop = FALSE]
    rownames(df) <- NULL
    df
}

# The Spectra variable MsBackendParquet takes and serves for each storage
# column msPurity uses; other columns (`x_mspurity_*`) keep their names.
.MZS_SPECTRA_VARS <- c(id = "spectrumId", data_origin = "dataOrigin",
                       ms_level = "msLevel", selected_ion_mz = "precursorMz",
                       scan_polarity = "polarity",
                       spectrum_representation = "centroided",
                       collision_energy = "collisionEnergy",
                       acquisition_num_ = "acquisitionNum")

#' Storage-vocabulary metadata as Spectra variables: renamed by
#' `.MZS_SPECTRA_VARS`, polarity +1 / -1 as 1 / 0, and the centroid and
#' profile CURIEs as `TRUE` / `FALSE`.
#'
#' @noRd
.mzs_spectra_vars <- function(df) {
    if (!is.null(df$scan_polarity))
        df$scan_polarity <- c(0L, 1L)[match(df$scan_polarity, c(-1L, 1L))]
    if (!is.null(df$spectrum_representation))
        df$spectrum_representation <- c(TRUE, FALSE)[match(
            df$spectrum_representation, c("MS:1000127", "MS:1000128"))]
    hit <- names(df) %in% names(.MZS_SPECTRA_VARS)
    names(df)[hit] <- .MZS_SPECTRA_VARS[names(df)[hit]]
    df
}

#' The inverse of `.mzs_spectra_vars()` for one column.
#'
#' @noRd
.mzs_storage_value <- function(column, v) {
    switch(column,
           scan_polarity = c(-1L, 1L)[match(v, c(0L, 1L))],
           spectrum_representation = c("MS:1000127", "MS:1000128")[
               match(v, c(TRUE, FALSE))],
           v)
}

#' Spectrum metadata of a dataset's runs, read through MsBackendParquet.
#'
#' Native and mzPeak-backed datasets are read alike. Columns are named and
#' encoded as stored: `scan_polarity` is +1 / -1.
#'
#' @param columns storage columns to read; any the dataset lacks is `NA`.
#'
#' @return data.frame with `run_id`, `spectrum_id_` and `columns`, in
#'     manifest run order and `spectrum_id_` order within each run.
#'
#' @noRd
.mzs_read_spectra_meta <- function(dataset, m, columns, run_ids = NULL) {
    rf <- .mzs_runs_frame(m)
    if (!is.null(run_ids))
        rf <- rf[rf$run_id %in% run_ids, , drop = FALSE]
    want <- unique(c("spectrum_id_", columns))
    ids <- unlist(lapply(seq_len(nrow(rf)), function(i)
        rf$uid_base[i] + seq_len(rf$n_spectra[i]) - 1L))
    out <- data.frame(run_id = rep(rf$run_id, rf$n_spectra),
                      stringsAsFactors = FALSE)
    if (!length(ids)) {
        for (c in want) out[[c]] <- logical()
        return(out)
    }
    be <- Spectra::backendInitialize(MsBackendParquet::MsBackendParquet(),
                                     path = dataset)
    be <- be[as.character(ids)]
    vars <- ifelse(want %in% names(.MZS_SPECTRA_VARS),
                   .MZS_SPECTRA_VARS[want], want)
    have <- vars %in% Spectra::spectraVariables(be)
    sd <- as.data.frame(Spectra::spectraData(be, columns = vars[have]),
                        optional = TRUE)
    for (k in seq_along(want))
        out[[want[k]]] <- if (have[k])
            .mzs_storage_value(want[k], sd[[vars[k]]])
        else rep(NA, length(ids))
    out
}

# ---------------------------------------------------------------------------
# Sources.
# ---------------------------------------------------------------------------

#' @noRd
.mzs_check_key <- function(key, what = "A source key") {
    if (!is.character(key) || length(key) != 1L || is.na(key) ||
        !grepl(.MZS_KEY_RE, key))
        .mzs_abort("format", what, " must match [A-Za-z0-9._-]{1,128}; ",
                   "got '", paste(key, collapse = ", "), "'.")
    if (identical(key, "self"))
        .mzs_abort("format", "'self' is reserved for this dataset and ",
                   "cannot be declared as a source key ([R-020]).")
    invisible(key)
}

#' An absolute path with symbolic links resolved, for a path that may not
#' exist yet.
#'
#' @noRd
.mzs_full_path <- function(p) {
    if (file.exists(p))
        return(normalizePath(p))
    parent <- dirname(p)
    if (identical(parent, p))
        return(p)
    file.path(.mzs_full_path(parent), basename(p))
}

#' `target` relative to the directory `from`, so a study and its results
#' can be moved together; absolute when they share only the filesystem root.
#'
#' @noRd
.mzs_relative_path <- function(target, from) {
    t <- strsplit(.mzs_full_path(target), "/", fixed = TRUE)[[1L]]
    f <- strsplit(.mzs_full_path(from), "/", fixed = TRUE)[[1L]]
    n <- 0L
    while (n < min(length(t), length(f)) && identical(t[n + 1L], f[n + 1L]))
        n <- n + 1L
    if (n <= 1L)
        return(.mzs_full_path(target))
    paste(c(rep("..", length(f) - n), t[seq_len(length(t) - n) + n]),
          collapse = "/")
}

#' Per-run fingerprints `(run_id, uid_base, n_spectra, ingested_at)` of the
#' dataset at `path`.
#'
#' @noRd
.mzs_fingerprints <- function(path, run_ids = NULL) {
    rf <- MsBackendParquet::manifestRuns(path)
    if (!is.null(run_ids))
        rf <- rf[rf$run_id %in% run_ids, , drop = FALSE]
    lapply(seq_len(nrow(rf)), function(i)
        list(run_id = rf$run_id[i], uid_base = as.integer(rf$uid_base[i]),
             n_spectra = as.integer(rf$n_spectra[i]),
             ingested_at = as.integer(rf$ingested_at[i])))
}

#' A `resolution: "dataset"` source entry.
#'
#' The source is never written: one lacking a `uid` is identified by its
#' path and fingerprints alone.
#'
#' @param dataset the results dataset's final path, for a relative `path`.
#'
#' @noRd
.mzs_dataset_source <- function(key, source, dataset, role = "study",
                                run_ids = NULL) {
    .mzs_check_key(key)
    m <- MsBackendParquet::readManifest(source)
    entry <- list(key = key, role = role, resolution = "dataset",
                  uid = m$uid, path = .mzs_relative_path(source, dataset),
                  generation = m$generation,
                  n_spectra = as.integer(sum(.mzs_runs_frame(m)$n_spectra)),
                  runs = .mzs_fingerprints(source, run_ids))
    if (is.null(m$uid))
        entry["uid"] <- list(NULL)
    if (!is.null(m$library$version))
        entry$version <- m$library$version
    if (!is.null(m$library$digest))
        entry$digest <- m$library$digest
    entry
}

#' A `resolution: "external"` source entry: no Parquet dataset exists, so
#' only native ids and USIs can be followed.
#'
#' @noRd
.mzs_external_source <- function(key, location = NULL, role = "study",
                                 version = NULL) {
    .mzs_check_key(key)
    entry <- list(key = key, role = role, resolution = "external")
    entry["uid"] <- list(NULL)
    entry["path"] <- list(NULL)
    entry["generation"] <- list(NULL)
    if (!is.null(location) && !is.na(location))
        entry$location <- as.character(location)
    if (!is.null(version) && !is.na(version))
        entry$version <- as.character(version)
    entry["runs"] <- list(NULL)
    entry
}

#' Merge source entries into a manifest's `sources`, by key.
#'
#' A key already declared must describe the same dataset; new runs are
#' added to its fingerprints, and a changed fingerprint is `Stale`.
#'
#' @noRd
.mzs_merge_sources <- function(existing, new) {
    existing <- existing %||% list()
    keys <- vapply(existing, function(s) s$key, character(1))
    for (s in new) {
        i <- match(s$key, keys)
        if (is.na(i)) {
            existing[[length(existing) + 1L]] <- s
            keys <- c(keys, s$key)
            next
        }
        old <- existing[[i]]
        if (!identical(old$resolution, s$resolution) ||
            !identical(old$uid, s$uid))
            .mzs_abort("semantic", "Source '", s$key, "' is already ",
                       "declared for a different dataset.")
        have <- vapply(old$runs %||% list(), function(r) r$run_id,
                       character(1))
        for (r in s$runs %||% list()) {
            j <- match(r$run_id, have)
            if (is.na(j)) {
                old$runs[[length(old$runs) + 1L]] <- r
            } else if (!.mzs_same_fingerprint(old$runs[[j]], r)) {
                .mzs_abort("stale", "Run '", r$run_id, "' of source '",
                           s$key, "' has changed since this dataset first ",
                           "referenced it.",
                           data = list(source = s$key, run_id = r$run_id))
            }
        }
        existing[[i]] <- old
    }
    existing
}

#' @noRd
.mzs_same_fingerprint <- function(a, b) {
    f <- c("run_id", "uid_base", "n_spectra", "ingested_at")
    identical(lapply(a[f], as.character), lapply(b[f], as.character))
}

#' @noRd
.mzs_source_entry <- function(m, key) {
    for (s in m$sources %||% list())
        if (identical(s$key, key))
            return(s)
    .mzs_abort("format", "A reference names source '", key, "', which the ",
               "manifest does not declare ([R-018]).")
}

#' Locate a `dataset` source: a caller-supplied override, then `path` as
#' given, then `path` relative to this dataset.
#'
#' @return the source's directory, or `NULL`.
#'
#' @noRd
.mzs_locate_source <- function(dataset, entry, sourcePaths = NULL) {
    cand <- c(if (!is.null(sourcePaths[[entry$key]]))
                  sourcePaths[[entry$key]],
              entry$path)
    if (!is.null(entry$path) && !.mzs_is_absolute(entry$path))
        cand <- c(cand, file.path(dataset, entry$path))
    for (p in cand)
        if (!is.null(p) && .mzs_is_dataset(p))
            return(normalizePath(p))
    NULL
}

#' Open a declared source for following references.
#'
#' An `external` source is `Capability`; a `dataset` source that cannot be
#' located or read, or has another `uid`, is `Reference`; a referenced run
#' whose fingerprint differs is `Stale`.
#'
#' @param run_ids the runs about to be followed; their fingerprints are
#'     checked.
#'
#' @return list(path, manifest).
#'
#' @noRd
.mzs_open_source <- function(dataset, m, key, sourcePaths = NULL,
                             run_ids = NULL) {
    if (identical(key, "self"))
        return(list(path = dataset, manifest = m))
    entry <- .mzs_source_entry(m, key)
    if (identical(entry$resolution, "external"))
        .mzs_abort("capability", "Source '", key, "' is external: no Parquet ",
                   "dataset holds its spectra ([R-025]).",
                   data = list(source = key))
    p <- .mzs_locate_source(dataset, entry, sourcePaths)
    if (is.null(p))
        .mzs_abort("reference", "Source '", key, "' cannot be located (",
                   "declared path '", entry$path %||% "<none>", "'). Pass ",
                   "its location in 'sourcePaths'.",
                   data = list(source = key))
    sm <- tryCatch(MsBackendParquet::readManifest(p),
                   mzstack_format = function(e)
        .mzs_abort("reference", "Source '", key, "' was located at '", p,
                   "' but is not a readable Parquet dataset: ",
                   conditionMessage(e), data = list(source = key)))
    if (!is.null(entry$uid) && !identical(entry$uid, sm$uid))
        .mzs_abort("reference", "Source '", key, "' at '", p, "' is a ",
                   "different dataset (uid ", sm$uid %||% "<none>",
                   ", expected ", entry$uid, ").", data = list(source = key))
    if (!identical(entry$generation, sm$generation)) {
        want <- entry$runs %||% list()
        if (!is.null(run_ids))
            want <- Filter(function(r) r$run_id %in% run_ids, want)
        have <- .mzs_fingerprints(p)
        names(have) <- vapply(have, function(r) r$run_id, character(1))
        for (r in want) {
            h <- have[[r$run_id]]
            if (is.null(h) || !.mzs_same_fingerprint(h, r))
                .mzs_abort("stale", "Run '", r$run_id, "' of source '", key,
                           "' no longer matches the fingerprint this ",
                           "dataset recorded ([R-029]).",
                           data = list(source = key, run_id = r$run_id))
        }
    }
    list(path = p, manifest = sm)
}

#' The resolution status of each row's spectrum reference.
#'
#' Unfollowable references are reported, not dropped. A stale source raises
#' `Stale`.
#'
#' @return character vector: `resolved`, `unresolved_external`,
#'     `unresolved_missing`.
#'
#' @noRd
.mzs_reference_status <- function(dataset, m, sources, run_ids,
                                  sourcePaths = NULL) {
    out <- rep("resolved", length(sources))
    for (k in unique(sources[!is.na(sources)])) {
        i <- which(sources == k)
        st <- tryCatch({
            .mzs_open_source(dataset, m, k, sourcePaths,
                             unique(run_ids[i]))
            "resolved"
        }, mzstack_capability = function(e) "unresolved_external",
        mzstack_reference = function(e) "unresolved_missing")
        out[i] <- st
    }
    out
}

# ---------------------------------------------------------------------------
# Provenance.
# ---------------------------------------------------------------------------

#' The next activity id, `act-0001` onwards.
#'
#' @noRd
.mzs_next_activity_id <- function(m) {
    ids <- vapply(m$provenance %||% list(), function(p)
        as.character(p$id)[1L], character(1))
    n <- suppressWarnings(max(c(0L, as.integer(sub("^act-", "", ids))),
                              na.rm = TRUE))
    sprintf("act-%04d", n + 1L)
}

#' @noRd
.mzs_environment <- function() {
    list(language = "R",
         language_version = paste(R.version$major, R.version$minor,
                                  sep = "."),
         platform = R.version$platform)
}

#' The agent recorded for an msPurity activity.
#'
#' @noRd
.mzs_agent <- function(fn) {
    list(name = "msPurity", version = .mzs_pkg_version("msPurity"),
         `function` = fn)
}

#' One activity record. Parameters are recorded verbatim under the tool's
#' own names.
#'
#' @noRd
.mzs_activity <- function(id, generation, started, action, agent, inputs,
                          outputs, parameters) {
    list(id = id, generation = as.integer(generation), started = started,
         ended = .mzs_now(), action = action, agent = agent,
         environment = .mzs_environment(), inputs = .mzs_array(inputs),
         outputs = list(runs = .mzs_array(outputs$runs %||% character()),
                        tables = .mzs_array(outputs$tables %||%
                                                character())),
         parameters = .mzs_json_safe(parameters %||% .mzs_object()))
}

#' An input that is a bare file, identified by its digest.
#'
#' @noRd
.mzs_file_input <- function(path) {
    list(path = normalizePath(path), bytes = as.numeric(file.size(path)),
         sha256 = .mzs_sha256(path))
}

#' @noRd
.mzs_sha256 <- function(path) unname(tools::sha256sum(path))

#' @noRd
.mzs_sha256_text <- function(x) {
    tmp <- tempfile()
    on.exit(unlink(tmp))
    writeBin(charToRaw(enc2utf8(x)), tmp)
    .mzs_sha256(tmp)
}

#' Make an R value serialisable as JSON without losing its meaning.
#'
#' A vector of length other than one becomes a JSON array; numbers are left
#' as is and rendered exactly when written.
#'
#' @noRd
.mzs_json_safe <- function(x) {
    if (is.null(x))
        return(NULL)
    if (isS4(x)) {
        sl <- methods::slotNames(x)
        out <- lapply(sl, function(s) .mzs_json_safe(methods::slot(x, s)))
        names(out) <- sl
        return(c(list(class = class(x)[1L]), out))
    }
    if (is.function(x))
        return(paste(deparse(x), collapse = "\n"))
    if (is.environment(x) || is.language(x))
        return(paste(deparse(x), collapse = "\n"))
    if (is.factor(x))
        x <- as.character(x)
    if (is.data.frame(x))
        x <- as.list(x)
    if (is.list(x)) {
        out <- lapply(x, .mzs_json_safe)
        if (is.null(names(x)))
            return(unname(out))
        if (!length(out))
            return(.mzs_object())
        return(out)
    }
    if (is.atomic(x)) {
        if (inherits(x, c("POSIXt", "Date")))
            x <- format(x, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
        nms <- names(x)
        x <- as.vector(x)
        if (!is.null(nms))
            return(as.list(stats::setNames(x, nms)))
        if (length(x) != 1L)
            return(as.list(x))
        return(x)
    }
    as.character(x)
}

# ---------------------------------------------------------------------------
# Reading results tables.
# ---------------------------------------------------------------------------

#' A declared results table as a data.frame, with dictionary columns as
#' character and Hive partition columns restored.
#'
#' A table declared omitted, or not declared, is `Capability`.
#'
#' @noRd
.mzs_read_table <- function(path, m, name) {
    maj <- .mzs_semver_major(m$results$version)
    if (!identical(maj, .MZS_RESULTS_MAJOR))
        .mzs_abort("unsupported", "'", path, "' holds results version ",
                   m$results$version %||% "<missing>", "; msPurity reads ",
                   .MZS_RESULTS_MAJOR, ".x. Its spectra can still be read ",
                   "(R-035).", data = list(version = m$results$version))
    e <- m$results$tables[[name]]
    if (is.null(e) || isTRUE(e$omitted))
        .mzs_abort("capability", "The results dataset '", path, "' has no ",
                   "table '", name, "'.", data = list(table = name))
    root <- file.path(path, e$path)
    files <- list.files(root, pattern = "\\.parquet$", recursive = TRUE,
                        full.names = TRUE)
    if (!length(files))
        return(.mzs_empty_table(name))
    part <- unlist(e$partitioning %||% list())
    ds <- arrow::open_dataset(
        root, format = "parquet", unify_schemas = TRUE,
        partitioning = if (length(part))
            do.call(arrow::hive_partition, stats::setNames(
                rep(list(arrow::utf8()), length(part)), part))
        else NULL)
    df <- as.data.frame(dplyr::collect(ds))
    for (c in names(df))
        if (is.factor(df[[c]]))
            df[[c]] <- as.character(df[[c]])
    df
}

#' A zero-row data.frame with a table's columns.
#'
#' @noRd
.mzs_empty_table <- function(name) {
    cols <- .mzs_table_schema(name)$columns
    df <- as.data.frame(stats::setNames(
        lapply(cols$type, function(t) .mzs_as_r_type(logical(), t)),
        cols$name), stringsAsFactors = FALSE)
    df
}
