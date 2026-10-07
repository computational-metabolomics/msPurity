# Spectral matching speed: SQLite route vs mzStack (Parquet) route.
#
# Matches the averaged spectra of the example database against a small
# library and against the full default library (Zenodo record 18700802,
# fetched through BiocFileCache), on both routes, and checks that the scores
# are identical. Also times the one-off conversions and reports disk sizes.
#
# Usage:
#
#   Rscript inst/benchmarks/spectral-matching.R [outDir] [libraryPath]
#
# outDir defaults to a temporary directory; it receives the converted
# datasets, results.csv and results.rds. libraryPath defaults to the cached
# default library, downloaded (about 480 MB) on first use. A full run takes
# about 20 minutes, almost all of it SQLite-route matching against the full
# library.

suppressPackageStartupMessages(library(msPurity))
for (p in c("arrow", "jsonlite", "BiocFileCache"))
    if (!requireNamespace(p, quietly = TRUE))
        stop("The benchmark needs the package '", p, "'.", call. = FALSE)

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args) >= 1L) args[1] else tempfile("sm-bench-")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out)

lib_full <- if (length(args) >= 2L) args[2] else {
    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
    hit <- BiocFileCache::bfcquery(bfc, "msPurity_library_spectra_db",
                                   "rname", exact = TRUE)
    if (!nrow(hit)) {
        BiocFileCache::bfcadd(bfc, "msPurity_library_spectra_db",
            "https://zenodo.org/records/18700802/files/library_spectra.db?download=1")
        hit <- BiocFileCache::bfcquery(bfc, "msPurity_library_spectra_db",
                                       "rname", exact = TRUE)
    }
    BiocFileCache::bfcrpath(bfc, rids = hit$rid[1])
}
lib_mini <- system.file("extdata", "tests", "mzstack", "library",
                        "mini_library.sqlite", package = "msPurity")
q_db <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite",
                    package = "msPurity")

log <- function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
    flush.console()
}
dir_size <- function(p)
    sum(file.info(list.files(p, recursive = TRUE, full.names = TRUE))$size)
res <- list()
add <- function(step, route, seconds = NA_real_, value = NA_real_) {
    res[[length(res) + 1L]] <<- data.frame(step = step, route = route,
                                           seconds = seconds, value = value,
                                           stringsAsFactors = FALSE)
}
timed <- function(expr) system.time(expr)[["elapsed"]]

# One-off conversions.
qd <- file.path(out, "query.mzstack")
lmd <- file.path(out, "mini-library.mzstack")
lfd <- file.path(out, "full-library.mzstack")
unlink(c(qd, lmd, lfd), recursive = TRUE)
for (x in list(list("query database", quote(convertSqliteToMzstack(q_db, qd))),
               list("mini library", quote(convertLibraryToMzstack(lib_mini, lmd))),
               list("full library", quote(convertLibraryToMzstack(lib_full, lfd))))) {
    t <- timed(eval(x[[2]]))
    log("convert", x[[1]], t)
    add(paste("convert", x[[1]]), "mzstack", seconds = t)
}
add("bytes: query", c("sqlite", "mzstack"),
    value = c(file.size(q_db), dir_size(qd)))
add("bytes: full library", c("sqlite", "mzstack"),
    value = c(file.size(lib_full), dir_size(lfd)))

# mzStack query ids back to the SQLite pids, to compare scores.
si <- as.data.frame(dplyr::collect(arrow::open_dataset(
    file.path(qd, "results", "source_identifier"))))
sp <- si[si$target_table == "spectra", ]
pid_of <- function(id) as.numeric(sp$source_value[match(id, sp$target_key)])
same_scores <- function(a, b) {
    if (is.null(a) || is.null(b))
        return(is.null(a) && is.null(b))
    b$qpid <- pid_of(b$qpid)
    ka <- paste(a$qpid, a$library_accession)
    kb <- paste(b$qpid, b$library_accession)
    if (!setequal(ka, kb))
        return(FALSE)
    k <- match(ka, kb)
    all(vapply(c("dpc", "rdpc", "cdpc", "mcount", "allcount", "mpercent"),
               function(s) identical(as.numeric(a[[s]]),
                                     as.numeric(b[[s]][k])), logical(1)))
}

scenario <- function(name, lib_sqlite, lib_mzs, types, reps_sqlite, reps_mzs) {
    sq <- mz <- NULL
    for (i in seq_len(reps_sqlite)) {
        t <- timed(sq <- suppressMessages(spectralMatching(
            q_db, lib_sqlite, q_spectraTypes = types, cores = 1)))
        log(name, "sqlite", i, t)
        add(name, "sqlite", t, NROW(sq$matchedResults))
    }
    for (i in seq_len(reps_mzs)) {
        t <- timed(mz <- suppressMessages(spectralMatching(
            qd, lib_mzs, q_spectraTypes = types, cores = 1,
            format = "mzstack")))
        log(name, "mzstack", i, t)
        add(name, "mzstack", t, NROW(mz$matchedResults))
    }
    ok <- same_scores(sq$matchedResults, mz$matchedResults)
    log(name, "identical scores:", ok)
    add(paste(name, "(identical scores)"), "both", value = as.numeric(ok))
}

scenario("mini library, av_all + inter", lib_mini, lmd,
         c("av_all", "inter"), 5, 5)
scenario("full library, av_all", lib_full, lfd, "av_all", 1, 3)

# Writing the matches: appended to the mzStack dataset.
t <- timed(suppressMessages(spectralMatching(
    qd, lfd, q_spectraTypes = "av_all", cores = 1, format = "mzstack",
    updateDb = TRUE)))
log("full library, av_all, updateDb", t)
add("full library, av_all, updateDb", "mzstack", seconds = t)

res <- do.call(rbind, res)
utils::write.csv(res, file.path(out, "results.csv"), row.names = FALSE)
saveRDS(res, file.path(out, "results.rds"))
timing <- res[!is.na(res$seconds), ]
print(stats::aggregate(seconds ~ step + route, timing, stats::median),
      row.names = FALSE)
log("results in", out)
