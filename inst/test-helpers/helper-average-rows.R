.av_rows <- function(path) {
    runs <- lapply(.manifest(path)$runs, function(r) {
        d <- .read_run(path, r$run_id)
        d$run_id <- r$run_id
        d[!vapply(d, is.list, TRUE)]
    })
    do.call(msPurity:::.mzs_rbind_fill, runs)
}
