.px_validate_parameters <- function(values) {
    modes <- c("purityFWHMmedian", "purityFWmedian", "purityFWHMmean",
        "purityFWmean")
    if (!is.character(values$purityType) || length(values$purityType) != 1L ||
            anyNA(values$purityType) || !values$purityType %in% modes)
        stop("purityType must name one documented purity prediction mode.")
    offsets <- values$offsets
    if (!is.numeric(offsets) || length(offsets) != 2L ||
            any(!is.finite(offsets)) || any(offsets < 0) || sum(offsets) <= 0)
        stop("offsets must contain nonnegative upper and lower offsets, ",
            "with positive total width.")
    for (name in c("iwNorm", "isotopes")) {
        value <- values[[name]]
        if (!is.logical(value) || length(value) != 1L || anyNA(value))
            stop(name, " must be TRUE or FALSE.")
    }
    if (!is.null(values$iwNormFun) && !is.function(values$iwNormFun))
        stop("iwNormFun must be a function or NULL.")
    ilim <- values$ilim
    if (!is.numeric(ilim) || length(ilim) != 1L || !is.finite(ilim) ||
            ilim < 0 || ilim > 1)
        stop("ilim must be between zero and one.")
    im <- values$im
    if (!is.null(im) && (!is.matrix(im) || !is.numeric(im) ||
                nrow(im) < 1L || ncol(im) != 8L || any(!is.finite(im))))
        stop("im must be NULL or a finite numeric matrix with eight columns.")
    invisible(values)
}
