# Measurements are supplied by the assessment, so history does not repeat its
# selection rules or retain peak arrays.
.pa_ms1_observations <- function(purity, metadata, used, method) {
    count <- lengths(used)
    method <- rep(method, length.out = nrow(purity))
    unavailable <- count == 0L & method != "invalid_isolation_window"
    method[unavailable] <- "no_usable_ms1"
    method[count == 1L] <- "single_ms1"
    result <- data.frame(seqNum = purity$seqNum,
        acquisitionNum = purity$acquisitionNum, method = method)
    result$ms1_scanIndex <- lapply(used, function(rows) metadata$seqNum[rows])
    result$ms1_acquisitionNum <- lapply(used,
        function(rows) metadata$acquisitionNum[rows])
    result
}
