# The Parquet study dataset a purityA object reads its raw spectra from.
#
# When purityA() is given a Spectra on MsBackendParquet, its spectra slot
# keeps that backend: peaks are read from the dataset, never from the mzML
# files it was converted from. The dataset's path and run fingerprints are
# recorded so that a moved or rewritten dataset is reported as such.

# Whether `sp` reads from a Parquet dataset.
.pa_is_parquet <- function(sp) {
    methods::is(sp, "Spectra") &&
        methods::is(sp@backend, "MsBackendParquet")
}

# The dataset behind `sp`: its path, uid and the fingerprints of the runs
# `sp` holds. NULL unless `sp` reads from a Parquet dataset.
.pa_study_identity <- function(sp) {
    if (!.pa_is_parquet(sp) || !length(sp))
        return(NULL)
    path <- unique(Spectra::dataStorage(sp))
    if (length(path) != 1L)
        stop("The spectra come from more than one Parquet dataset.")
    path <- normalizePath(path, mustWork = TRUE)
    m <- MsBackendParquet::readManifest(path)
    list(path = path, uid = m$uid %||% NA_character_,
         runs = .mzs_fingerprints(path, unique(as.character(sp$run_id))))
}

# Each run of a Parquet input must be one source file: msPurity numbers files
# by dataOrigin, and two runs converted from the same path would merge.
.pa_check_parquet_runs <- function(sp) {
    if (!"run_id" %in% Spectra::spectraVariables(sp))
        return(invisible(TRUE))
    pairs <- unique(data.frame(run = as.character(sp$run_id),
                               origin = as.character(Spectra::dataOrigin(sp))))
    if (anyDuplicated(pairs$run) || anyDuplicated(pairs$origin))
        stop("Each run of the Parquet dataset must come from one source ",
             "file, and each source file be one run; select the runs to ",
             "process with filterDataOrigin().")
    invisible(TRUE)
}

# The study recorded by purityA(), or NULL.
.pa_study <- function(pa) {
    pa@params$purityA$study
}

# Stop unless the study purityA() read is still where it was and its runs
# are unchanged: a missing dataset is a Reference error, a changed run Stale.
.pa_study_check <- function(pa) {
    study <- .pa_study(pa)
    if (is.null(study))
        return(invisible(NULL))
    m <- tryCatch(MsBackendParquet::readManifest(study$path),
                  error = function(e) NULL)
    if (is.null(m))
        .mzs_abort("reference", "purityA() read its spectra from the ",
                   "Parquet dataset '", study$path, "', which is no longer ",
                   "readable there. Restore it, or run purityA() again on ",
                   "the dataset where it is now.",
                   data = list(path = study$path))
    if (!identical(m$uid %||% NA_character_, study$uid))
        .mzs_abort("reference", "The Parquet dataset at '", study$path,
                   "' is not the one purityA() read (uid ",
                   m$uid %||% "<none>", ", expected ", study$uid, ").",
                   data = list(path = study$path))
    have <- .mzs_fingerprints(study$path)
    names(have) <- vapply(have, `[[`, "", "run_id")
    for (r in study$runs) {
        h <- have[[r$run_id]]
        if (is.null(h) || !.mzs_same_fingerprint(h, r))
            .mzs_abort("stale", "Run '", r$run_id, "' of the Parquet ",
                       "dataset '", study$path, "' changed after purityA() ",
                       "read it.", data = list(path = study$path,
                                               run_id = r$run_id))
    }
    invisible(study)
}

# The study createDatabase(format = "parquet") links scans to: `study` when
# given, none when NA, and by default the dataset purityA() read, if any.
.pa_study_resolve <- function(pa, study) {
    if (length(study) == 1L && is.na(study))
        return(NULL)
    if (!is.null(study))
        return(normalizePath(study, mustWork = TRUE))
    read <- .pa_study(pa)
    if (is.null(read))
        return(NULL)
    .pa_study_check(pa)
    message("Linking scans to the Parquet dataset purityA() read: ",
            read$path)
    read$path
}
