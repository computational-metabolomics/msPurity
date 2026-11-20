# Internal batch sizes are adjustable for regression tests and benchmarks.
.mzs_batch_size <- function(spectra = FALSE) {
    n <- if (spectra) getOption("msPurity.parquet.batch_spectra", 128L)
         else getOption("msPurity.parquet.batch_rows", 65536L)
    if (length(n) != 1L || !is.numeric(n) || !is.finite(n) ||
        n < 1 || n > .Machine$integer.max || n != floor(n))
        stop("The Parquet batch size must be a positive integer.")
    as.integer(n)
}

.mzs_batches <- function(n, size, visit) {
    if (n) for (start in seq.int(1L, n, by = size))
        visit(seq.int(start, min(n, as.double(start) + size - 1)))
    invisible(NULL)
}

# Spool a table whose rows are produced alongside other outputs. Replaying
# batches uses the schema and null counts observed over the entire table.
.mzs_spool <- function(dir, name) {
    path <- file.path(dir, name)
    dir.create(path)
    parts <- 0L
    rows <- 0
    present <- character()
    prototype <- list()
    pending <- list()
    buffered <- 0L
    size <- .mzs_batch_size()
    flush <- function() {
        if (!buffered) return()
        parts <<- parts + 1L
        saveRDS(do.call(.mzs_rbind_fill, pending), file.path(path, parts), compress = FALSE)
        pending <<- list()
        buffered <<- 0L
    }
    add <- function(df) {
        if (is.null(df) || !nrow(df)) return(invisible(NULL))
        for (col in names(df))
            if (is.null(prototype[[col]]) || is.logical(prototype[[col]]) && !is.logical(df[[col]]))
                prototype[[col]] <<- df[[col]][FALSE]
        nonnull <- vapply(df, function(v) any(!is.na(v) | is.nan(v)), logical(1))
        present <<- union(present, names(df)[nonnull])
        first <- 1L
        while (first <= nrow(df)) {
            count <- min(size - buffered, nrow(df) - first + 1L)
            pending[[length(pending) + 1L]] <<- df[seq.int(first, length.out = count), , drop = FALSE]
            buffered <<- buffered + count
            first <- first + count
            if (buffered == size) flush()
        }
        rows <<- rows + nrow(df)
    }
    finish <- function() {
        flush()
        if (!parts) return(NULL)
        visit <- function(f) {
            for (i in seq_len(parts)) f(readRDS(file.path(path, i)))
        }
        structure(list(rows = rows, prototype = prototype, present = present,
                       visit = visit), class = "mzs_table_batches")
    }
    list(add = add, finish = finish)
}

.mzs_map_table <- function(dir, name, n, build) {
    if (is.null(dir) || !n) return(build(seq_len(n)))
    parts <- .mzs_spool(dir, name)
    .mzs_batches(n, .mzs_batch_size(), function(i) parts$add(build(i)))
    parts$finish()
}

.mzs_source_id_stream <- function(dir) {
    parts <- .mzs_spool(dir, "source_identifier")
    count <- 0L
    add <- function(target_table, target_key, source_table, source_column, value, type) {
        if (is.null(value)) return()
        .mzs_batches(length(target_key), .mzs_batch_size(), function(i) {
            i <- i[!is.na(value[i])]
            if (!length(i)) return()
            df <- .mzs_source_ids(target_table, target_key[i], source_table,
                                  source_column, value[i], type)
            df$source_identifier_id_ <- count + seq_along(i)
            count <<- count + length(i)
            parts$add(df)
        })
    }
    finish <- function(...) {
        for (i in seq_len(...length())) force(...elt(i))
        parts$finish()
    }
    list(add = add, finish = finish)
}
