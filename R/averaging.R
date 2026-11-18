# msPurity R package for processing MS/MS data - Copyright (C)
#
# This file is part of msPurity.
#
# msPurity is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msPurity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msPurity.  If not, see <https://www.gnu.org/licenses/>.

group1d <- function(v, thr){

  v <- sort(v)
  a <- ''
  cl <-  1
  c <- 1

  for(i in 1:(length(v)-1)){
    a[i] <- (1e6*(v[i+1]-v[i]))/v[i]

    if(as.numeric(a[i])<thr){
      cl[i+1] <- c
    }else{
      c <- c+1
      cl[i+1] <- c
    }

  }
  return(cl)
}

# Clustering and batch boundaries share the same distance and cutoff.
.hc_distance <- function(difference, mean_mz, ppm, tolerance) {
  if (tolerance > 0) {
    # Absolute and relative tolerances add, with the numerical margin of
    # MsCoreUtils::closest().
    allowed <- tolerance + mean_mz * ppm * 0.000001 +
      sqrt(.Machine$double.eps)
    return(list(distance = difference / allowed, height = 1))
  }
  list(distance = difference / (mean_mz * 0.000001), height = ppm)
}

performHc <- function(mzs, ppm, tolerance = 0){
  if (length(mzs) == 1L)
    return(1L)
  span <- .hc_distance(mzs[length(mzs)] - mzs[1L],
    (mzs[1L] + mzs[length(mzs)]) / 2, ppm, tolerance)
  if (span$distance <= span$height)
    return(rep.int(1L, length(mzs)))
  distance <- .hc_distance(stats::dist(mzs),
    stats::as.dist(outer(mzs, mzs, "+") / 2), ppm, tolerance)
  clh <- fastcluster::hclust(distance$distance, method = "complete")
  stats::cutree(clh, h = distance$height)
}

# Sorted positive masses have increasing distance across a gap. Thus a gap
# above the cutoff separates every pair on opposite sides, including ties.
# Connected runs only locate safe boundaries; clustering remains complete-linkage.
.hc_blocks <- function(mz, ppm, tolerance, max_size = 5000L) {
  n <- length(mz)
  if (n <= max_size)
    return(data.frame(first = 1L, last = n, one_cluster = FALSE))

  adjacent <- .hc_distance(diff(mz),
    (mz[-n] + mz[-1L]) / 2, ppm, tolerance)
  ends <- c(which(adjacent$distance > adjacent$height), n)
  # One vectorized lookup avoids repeatedly scanning all safe boundaries.
  end_indices <- findInterval(c(0L, ends[-length(ends)]) + max_size, ends)
  blocks <- list()
  first <- 1L
  first_run <- 1L
  while (first <= n) {
    end_index <- end_indices[first_run]
    one_cluster <- FALSE
    if (end_index >= first_run) {
      last <- ends[end_index]
    } else {
      end_index <- first_run
      last <- ends[end_index]
      span <- .hc_distance(mz[last] - mz[first],
        (mz[first] + mz[last]) / 2, ppm, tolerance)
      if (span$distance > span$height)
        stop("A tolerance-connected group of ", last - first + 1L,
          " peaks exceeds the ", max_size,
          "-peak clustering limit and cannot be split safely.",
          call. = FALSE)
      one_cluster <- TRUE
    }
    blocks[[length(blocks) + 1L]] <- data.frame(first = first,
      last = last, one_cluster = one_cluster)
    first <- last + 1L
    first_run <- end_index + 1L
  }
  do.call(rbind, blocks)
}

clustering <- function(mz, ppm=1.5, clustType="hc", cores=2, tolerance=0){
  if (clustType == "hc") {
    if (!is.numeric(ppm) || length(ppm) != 1L || !is.finite(ppm) || ppm < 0)
      stop("'ppm' must be a finite, non-negative numeric scalar.", call. = FALSE)
    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        !is.finite(tolerance) || tolerance < 0)
      stop("'tolerance' must be a finite, non-negative numeric scalar.",
        call. = FALSE)
    if (!is.numeric(mz) || any(!is.finite(mz)) || any(mz <= 0) ||
        is.unsorted(mz))
      stop("Hierarchical clustering requires sorted, positive, finite m/z values.",
        call. = FALSE)
    if (!length(mz))
      return(integer())

    blocks <- .hc_blocks(mz, ppm, tolerance)
    clustlist <- vector("list", nrow(blocks))
    for (i in which(blocks$one_cluster))
      clustlist[[i]] <- rep.int(1L, blocks$last[i] - blocks$first[i] + 1L)

    regular <- which(!blocks$one_cluster)
    if (length(regular)) {
      mzsplit <- lapply(regular, function(i)
        mz[seq.int(blocks$first[i], blocks$last[i])])
      workers <- min(cores, length(regular))
      if (workers > 1L) {
        clust <- parallel::makeCluster(workers)
        on.exit(parallel::stopCluster(clust), add = TRUE)
        clustlist[regular] <- parallel::parLapply(clust, mzsplit, performHc,
          ppm = ppm, tolerance = tolerance)
      } else {
        clustlist[regular] <- lapply(mzsplit, performHc,
          ppm = ppm, tolerance = tolerance)
      }
    }

    offset <- 0L
    for (i in seq_along(clustlist)) {
      clustlist[[i]] <- clustlist[[i]] + offset
      offset <- max(clustlist[[i]])
    }
    return(unlist(clustlist, use.names = FALSE))
  }else if (clustType=="simple"){
    # order by mz
    #print("group1d")
    return(group1d(mz, ppm))
  }


}

stde <- function(x) sd(x)/sqrt(length(x))
rsde <- function(x) (sd(x)/mean(x))*100


averageCluster <- function(x, av="median", minnum=1,
                           missingV="ignore", totalScans, normTIC, sumI=FALSE){
  # Filter out any that do not match the minimum
  if(nrow(x)<minnum){
    return(NULL)
  }


  # Get rsd (if normalise flagged then use the normalised intensity)
  if(normTIC){
    rsdRes <- (sd(x$inorm)/mean(x$inorm))*100
  }else{
    rsdRes <- rsde(x$i)

  }

  # Calc average first of the mz. We don't want to mess around with
  # missing values for this stage
  if(av=="median"){
    mz <- stats::median(x$mz)
  }else{
    mz <- base::mean(x$mz)
  }

  # Zero any missing values for intensity and SNR (this is the thermo approach,
  # but might not be best approach)
  if(missingV=="zero"){
    toadd <- totalScans-nrow(x)
    if(toadd>0){
      for(i in 1:toadd){
        x <- rbind(x, c(mz,0,0,unique(x$cl)))
      }

    }
  }

  # Finally calculate the averaged SNR and intensity
  if(av=="median"){
    i <- median(x$i)
    snr <- median(x$snr)
    inorm <- median(x$inorm)
    purity <- median(x$inPurity)
    ra <- median(x$ra)
  }else{
    i <- mean(x$i)
    snr <- mean(x$snr)
    inorm <- mean(x$inorm)
    purity <- mean(x$inPurity)
    ra <- mean(x$ra)
  }

  if (sumI){
    i <- sum(x$i)
  }



  return(c("mz" = mz, "i" = i, "snr" = snr, "rsd" = rsdRes,
           "inorm" = inorm, 'count' =length(x$i), 'total'=totalScans, 'inPurity'=purity, 'ra'=ra))

}

snrFilter <- function(x, snthr, snMeth){
  # If snMethod is "precalc" then the SNR does not need to
  # be calculated required
  if (snMeth=="median"){
    x$snr <- x$i/median(x$i)
  }else if(snMeth=="mean"){
    x$snr <- x$i/mean(x$i)
  }

  x <- x[x$snr>snthr, ]

  return(x)

}
