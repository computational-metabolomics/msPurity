.blink_scores <- c("blinkCosine", "blinkMatchingIonCount")
.blink_revision <- "dfdf0c13672b9f5e5bbf59187ea4f6cd04627c64"

.sm_scoring <- function(param) {
    blink <- is(param, "BlinkMatchingParam")
    list(method = if (blink) "BLINK" else "msPurity",
        implementation = if (blink) "msPurity/Matrix" else "msPurity",
        version = .mzs_pkg_version("msPurity"),
        referenceRevision = if (blink) .blink_revision else NULL,
        settings = if (blink) as.list(param) else c(as.list(param), list(raW = 0.5, mzW = 2)),
        rankScore = if (blink) "blinkCosine" else "dpc")
}

# References form the columns of a bounded sparse block. Compress the mass
# axis to occupied reference bins so memory is independent of the m/z span.
.sm_blink_block <- function(query, references, param) {
    normalize <- function(p) {
        p <- p[p$i > 0, , drop = FALSE]
        # Scale in log space so very unequal finite intensities survive a
        # fractional power even when their unpowered ratio underflows.
        logIntensity <- log(p$i)
        w <- exp((logIntensity - max(logIntensity)) * param@intensityPower)
        list(bin = round(p$mz / param@binWidth), weight = w / sqrt(sum(w * w)))
    }
    q <- normalize(query)
    ref <- lapply(references, normalize)
    refBins <- lapply(ref, `[[`, "bin")
    bin <- unlist(refBins, use.names = FALSE)
    bins <- sort(unique(bin))
    # Upstream kernel length truncates 2*tolerance/binWidth - 1. An
    # even-length kernel has one more positive offset than negative offsets.
    width <- trunc(2 * (param@blurTolerance / param@binWidth) - 1)
    if (!is.finite(width) || any(!is.finite(c(bins, q$bin))))
        stop("BLINK bin indices overflow. Increase binWidth.")
    lo <- findInterval(q$bin + 1 - ceiling(width / 2), bins, left.open = TRUE) + 1L
    hi <- findInterval(q$bin + floor(width / 2), bins)
    use <- which(lo <= hi)
    n <- length(ref)
    if (!length(use))
        return(data.frame(blinkCosine = numeric(n), blinkMatchingIonCount = numeric(n)))
    rows <- sequence(hi[use] - lo[use] + 1L, from = lo[use])
    ions <- rep(use, hi[use] - lo[use] + 1L)
    dims <- c(length(bins), n)
    ri <- match(bin, bins)
    rj <- rep(seq_along(ref), lengths(refBins))
    rweight <- unlist(lapply(ref, `[[`, "weight"), use.names = FALSE)
    dot <- function(qx, rx) {
        qmat <- Matrix::sparseMatrix(i = rows, j = rep.int(1L, length(rows)),
            x = qx, dims = c(length(bins), 1L))
        rmat <- Matrix::sparseMatrix(i = ri, j = rj, x = rx, dims = dims)
        as.numeric(Matrix::crossprod(qmat, rmat))
    }
    data.frame(blinkCosine = dot(q$weight[ions], rweight),
        blinkMatchingIonCount = dot(rep(1, length(ions)), rep(1, length(bin))))
}

.sm_fragment_scores <- function(query, refs, param) {
    if (is(param, "BlinkMatchingParam")) return(.sm_blink_block(query, refs, param))
    as.data.frame(do.call(rbind, lapply(refs, function(p)
        .sm_score_peaks(query, p, param@q_ppmProd, param@l_ppmProd,
            raW = 0.5, mzW = 2, q_daProd = param@q_daProd, l_daProd = param@l_daProd))))
}
