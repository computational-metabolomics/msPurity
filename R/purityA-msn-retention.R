# Retain fragmentation measurements independently of MS1-based MS2 purity.
.pa_msn_retained_spectra <- function(raw, fileid, puritydf) {
    level <- Spectra::msLevel(raw)
    selected <- which(level %in% 2:4)
    if (!length(selected))
        return(Spectra::Spectra(Spectra::MsBackendMemory()))
    sp <- raw[selected]
    native <- as.data.frame(Spectra::spectraData(sp,
            columns = intersect(c("scanIndex", "acquisitionNum", "rtime",
                    "precursorMz", "precScanNum", "precursorIntensity",
                    "dataOrigin"),
                Spectra::spectraVariables(sp))), optional = TRUE)
    rows <- if (nrow(puritydf)) match(
        paste(fileid[selected], native$scanIndex),
        paste(puritydf$fileid, puritydf$seqNum)) else
        rep(NA_integer_, length(selected))
    if (anyNA(rows[level[selected] == 2L]))
        stop("MS2 purity results do not cover the retained spectrum positions.")
    columns <- names(puritydf)
    df <- if (length(columns)) puritydf[rows, , drop = FALSE] else
        data.frame(row.names = seq_along(selected))
    higher <- which(level[selected] > 2L)
    if (!"pid" %in% names(df)) df$pid <- rep(NA_integer_, nrow(df))
    next_pid <- if (nrow(puritydf)) max(puritydf$pid) else 0L
    df$pid[higher] <- next_pid + seq_along(higher)
    if (!"fileid" %in% names(df)) df$fileid <- rep(NA_integer_, nrow(df))
    df$fileid[higher] <- fileid[selected[higher]]
    mapping <- c(seqNum = "scanIndex", acquisitionNum = "acquisitionNum",
        retentionTime = "rtime", precursorMZ = "precursorMz",
        precursorScanNum = "precScanNum",
        precursorIntensity = "precursorIntensity")
    for (field in names(mapping)) {
        value <- native[[mapping[[field]]]]
        if (is.null(value)) value <- rep(NA_real_, length(sp))
        if (!field %in% names(df)) df[[field]] <- value
        else df[[field]][higher] <- value[higher]
    }
    if ("filename" %in% names(df))
        df$filename[higher] <- basename(native$dataOrigin[higher])
    if ("id" %in% names(df)) df$id[higher] <- native$scanIndex[higher]
    for (field in names(df)) sp[[field]] <- df[[field]]
    offsets <- .pa_raw_peak_offsets(raw)
    sp$msPurityRawPeakOffset <- offsets$offset[selected]
    sp$msPurityRawPeakCount <- offsets$count[selected]
    sp$msPurityRawPeakIdBasis <- rep("source_reference", length(sp))
    purity <- if ("inPurity" %in% Spectra::spectraVariables(sp))
        sp$inPurity else rep(NA_real_, length(sp))
    sp$msPurityPurityStatus <- ifelse(level[selected] != 2L,
        "not_assessed", ifelse(is.na(purity), "unavailable", "assessed"))
    meta <- .pa_spectra_metadata(sp)
    meta$puritydfColumns <- columns
    if (all(Spectra::dataStorage(sp) == "<memory>"))
        meta$msPurityRawPeaksAvailable <- TRUE
    .pa_spectra_metadata(sp) <- meta
    sp
}
