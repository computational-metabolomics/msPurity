# Collection search shares the pair scorer with scoreSpectrumPair(). Row
# indices always refer to the supplied objects, even after metadata filters.
.sm_rank_scores <- function(scores) {
    rank(-scores, ties.method = "min", na.last = "keep")
}

.sm_validate_topn <- function(topn) {
    if (length(topn) != 1L || !is.na(topn) &&
            (!is.numeric(topn) || !is.finite(topn) ||
                topn < 0 || topn != trunc(topn)))
        stop("topn must be a nonnegative integer or NA.")
}

.sm_search_param <- function(a, supplied) {
    ppm <- c("q_ppmPrec", "l_ppmPrec", "q_ppmProd", "l_ppmProd")
    if (is.null(a$param)) return(do.call(SpectralMatchingParam, a[ppm]))
    if (!methods::is(a$param, "SpectralMatchingParam"))
        stop("param must be a SpectralMatchingParam object.")
    methods::validObject(a$param)
    for (name in intersect(ppm, supplied)) {
        x <- a[[name]]
        if (is(a$param, "BlinkMatchingParam") && endsWith(name, "Prod") && !isTRUE(x == 0))
            stop("BLINK uses blurTolerance in Da. Set ", name, " to zero or omit it.")
        if (!is.numeric(x) || length(x) != 1L ||
                !isTRUE(x == methods::slot(a$param, name)))
            stop(name, " conflicts with param.")
    }
    a$param
}

.sm_reference <- function(x, role) {
    if (!is.null(role))
        role <- match.arg(role, c("experimental", "library"))
    if (methods::is(x, "Spectra"))
        return(list(spectra = x, role = role %||% "experimental"))
    if (!.mzs_is_dataset(x))
        stop("l_dbPth must be a Spectra object or a Parquet library or results path.")
    dataset <- MsBackendParquet::readManifest(x)$role
    if (!isTRUE(dataset %in% c("library", "results")))
        stop("The reference path must identify a Parquet library or results dataset.")
    inferred <- if (dataset == "library") "library" else "experimental"
    if (!is.null(role) && role != inferred)
        stop("referenceRole conflicts with the reference dataset role.")
    sp <- if (dataset == "results") readProcessedSpectra(x) else
        Spectra::Spectra(Spectra::backendInitialize(
            MsBackendParquet::MsBackendParquet(), path = x))
    list(spectra = sp, role = inferred)
}

.sm_search_metadata <- function(sp) {
    vars <- c("msLevel", "precursorMz", "rtime", "polarity", "run_id", "dataOrigin",
        "fileid", "spectrumId", "nativeId", "acquisitionNum", "spectrum_id_",
        "x_mspurity_spectrum_key", "accession", "x_mspurity_splash", "splash",
        "x_mspurity_retention_time", "av_level")
    d <- Spectra::spectraData(sp, columns = intersect(vars, Spectra::spectraVariables(sp)))
    value <- function(...) .sm_metadata_value(d, ...)
    data.frame(level = d$msLevel, mz = d$precursorMz,
        rt = .sm_metadata_value(d, "x_mspurity_retention_time", "rtime", numeric = TRUE),
        polarity = unname(c("0" = "negative", "1" = "positive")[as.character(d$polarity)]),
        run = value("run_id", "dataOrigin", "fileid"),
        native = value("x_mspurity_spectrum_key", "spectrumId", "nativeId", "acquisitionNum"),
        id = value("spectrum_id_"),
        accession = value("accession", "spectrumId", "x_mspurity_splash", "splash"),
        consensus = !is.na(value("av_level", "x_mspurity_spectrum_key")),
        origin = value("dataOrigin"), runId = value("run_id"),
        nativeId = value("spectrumId", "nativeId"), scan = value("acquisitionNum"))
}

.sm_query_status <- function(meta, rows, type) {
    n <- nrow(meta)
    status <- rep("filtered", n)
    reason <- rep("query_filter", n)
    reason[!is.finite(meta$mz) | meta$mz <= 0] <- "missing_precursor"
    reason[!meta$level %in% 2:4] <- "ms_level_unsupported"
    reason[is.na(meta$level)] <- "ms_level_missing"
    status[reason != "query_filter"] <- "invalid_query"
    status[rows] <- "no_eligible_reference"
    reason[rows] <- "no_compatible_reference"
    data.frame(queryIndex = seq_len(n), query_run_id = meta$run,
        query_native_id = meta$native, query_spectrum_id_ = meta$id,
        queryType = type, msLevel = meta$level, status = status, reason = reason,
        nEligible = integer(n), nScored = integer(n), nRetained = integer(n))
}

.sm_search_rows <- function(meta, polarity, rt) {
    keep <- meta$level %in% 2:4 & is.finite(meta$mz) & meta$mz > 0
    if (length(polarity) != 1L || !is.na(polarity) &&
            !tolower(polarity) %in% c("positive", "negative"))
        stop("Polarity must be positive, negative or NA.")
    if (!is.na(polarity)) keep <- keep & meta$polarity %in% tolower(polarity)
    if (!all(is.na(rt))) {
        if (!is.numeric(rt) || length(rt) != 2L || anyNA(rt) || rt[1L] > rt[2L])
            stop("Retention-time ranges require ordered lower and upper bounds.")
        keep <- keep & meta$rt >= rt[1L] & meta$rt <= rt[2L]
    }
    which(keep)
}

.sm_candidate_index <- function(meta, rows, param) {
    by_level <- split(rows, meta$level[rows])
    windows <- lapply(by_level, function(ids) {
        ids <- ids[order(meta$mz[ids])]
        mz <- meta$mz[ids]
        width <- .sm_mass_width(mz, param@l_ppmPrec, param@l_daPrec)
        # Monotone bounds on the actual endpoints also work above 1e6 ppm
        # and at floating-point boundaries. The scorer checks exact overlap.
        list(ids = ids, upper = cummax(mz + width),
            lower = rev(cummin(rev(mz - width))))
    })
    function(level, mz) {
        w <- windows[[as.character(level)]]
        if (is.null(w)) return(integer())
        width <- .sm_mass_width(mz, param@q_ppmPrec, param@q_daPrec)
        first <- findInterval(mz - width, w$upper, left.open = TRUE, checkSorted = FALSE) + 1L
        last <- findInterval(mz + width, w$lower, checkSorted = FALSE)
        if (first > last) return(integer())
        sort.int(w$ids[seq.int(first, last)])
    }
}

.sm_stage_reader <- function(sp, parentData = .sm_parent_data(sp)) {
    cache <- new.env(parent = emptyenv())
    function(i, level) {
        key <- as.character(i)
        if (is.null(cache[[key]]))
            cache[[key]] <- .sm_stages(sp, i, level, parentData)
        cache[[key]]
    }
}

.sm_search_spectra <- function(a, supplied, queryRows = NULL, referenceRows = NULL,
        queryPeaks = NULL, referencePeaks = NULL, queryRt = NULL, referenceRt = NULL) {
    supported <- c("q_dbPth", "l_dbPth", "param", "referenceRole", "q_ppmPrec",
        "l_ppmPrec", "q_ppmProd", "l_ppmProd", "q_pol", "l_pol", "q_rtrange",
        "l_rtrange", "rttol", "l_accessions", "topn", "usePrecursors", "format",
        "updateDb", "copyDb")
    # These explicit settings describe the serial D04 scorer and the caller's
    # unfiltered peak selection. Other values require a different operation.
    if (is(a$param, "BlinkMatchingParam") && any(c("raW", "mzW") %in% supplied))
        stop("Set BLINK intensityPower in BlinkMatchingParam().")
    fixed <- c(cores = 1, raW = 0.5, mzW = 2)
    for (name in intersect(names(fixed), supplied)) {
        if (!is.numeric(a[[name]]) || !isTRUE(a[[name]] == fixed[[name]]))
            stop(name, " must be ", fixed[[name]], " for a Spectra query.")
    }
    flags <- intersect(c("q_spectraFilter", "l_spectraFilter"), supplied)
    for (name in flags)
        if (!identical(a[[name]], FALSE))
            stop(name, " is unsupported for a Spectra query. Filter the supplied peaks first.")
    unsupported <- setdiff(supplied, c(supported, names(fixed), flags))
    if (length(unsupported))
        stop("Unsupported for a Spectra query: ", paste(unsupported, collapse = ", "), ".")
    if (!identical(a$updateDb, FALSE) || !identical(a$copyDb, FALSE))
        stop("updateDb and copyDb must be FALSE for an in-memory Spectra search.")
    if ("format" %in% supplied && !identical(a$format, "parquet"))
        stop("A Spectra query cannot be searched with a SQLite format.")
    if (!isTRUE(a$usePrecursors))
        stop("usePrecursors must be TRUE for a Spectra query.")
    .sm_validate_topn(a$topn)
    if (length(a$rttol) != 1L || !is.na(a$rttol) &&
            (!is.numeric(a$rttol) || a$rttol < 0))
        stop("rttol must be nonnegative or NA.")
    param <- .sm_search_param(a, supplied)
    ref <- .sm_reference(a$l_dbPth, a$referenceRole)
    query <- a$q_dbPth
    reference <- ref$spectra
    accessions <- as.character(stats::na.omit(a$l_accessions))
    if (length(accessions) && ref$role != "library")
        stop("l_accessions requires a library reference.")
    qm <- .sm_search_metadata(query)
    lm <- .sm_search_metadata(reference)
    if (!is.null(queryRt)) qm$rt <- queryRt
    if (!is.null(referenceRt)) lm$rt <- referenceRt
    if (is.null(queryPeaks)) queryPeaks <- function(i) .sm_pair_peaks(query[i])
    if (is.null(referencePeaks)) referencePeaks <- function(ids)
        lapply(.msp_peaks(reference[ids]), .sm_clean_peaks)
    qrows <- .sm_search_rows(qm, a$q_pol, a$q_rtrange)
    lrows <- .sm_search_rows(lm, if ("l_pol" %in% supplied) a$l_pol else NA, a$l_rtrange)
    if (!is.null(queryRows)) qrows <- intersect(qrows, queryRows)
    if (!is.null(referenceRows)) lrows <- intersect(lrows, referenceRows)
    if (length(accessions)) {
        keep <- lm$accession %in% accessions |
            (is.na(lm$accession) & any(startsWith(accessions, "sha256:")))
        lrows <- lrows[keep[lrows]]
    }
    candidates <- .sm_candidate_index(lm, lrows, param)
    qstages <- .sm_stage_reader(query)
    lstages <- .sm_stage_reader(reference)
    qsource <- .sm_source_reader(query, qm)
    lsource <- .sm_source_reader(reference, lm, ref$role)
    status <- .sm_query_status(qm, qrows, qsource$type)

    scoring <- .sm_scoring(param)
    batchSize <- .mzs_batch_size(spectra = TRUE)
    match_query <- function(i, qp = queryPeaks(i)) {
        level <- qm$level[i]
        ids <- candidates(level, qm$mz[i])
        if (!is.na(a$rttol))
            ids <- ids[which(abs(lm$rt[ids] - qm$rt[i]) < a$rttol)]
        reasons <- character()
        parts <- list()
        eligible <- 0L
        .mzs_batches(length(ids), batchSize, function(ix) {
            j <- ids[ix]
            decision <- vapply(j, function(k) .sm_eligibility(level, lm$level[k],
                qm$mz[i], lm$mz[k], qstages(i, level), lstages(k, level), param)$reason, "")
            ok <- decision %in% c("compatible", "compatible_partial")
            reasons <<- union(reasons, decision[!ok])
            j <- j[ok]
            decision <- decision[ok]
            if (!length(j)) return()
            accession <- lm$accession[j]
            signal <- any(qp$i > 0)
            if (!signal && (!length(accessions) || !anyNA(accession))) {
                eligible <<- eligible + length(j)
                return()
            }
            peaks <- referencePeaks(j)
            usable <- signal & vapply(peaks, function(p) any(p$i > 0), logical(1))
            anonymous <- which(ref$role == "library" & is.na(accession) &
                (usable | length(accessions) > 0L))
            for (k in anonymous)
                accession[k] <- .mzs_anonymous_accession(lm$mz[j[k]], lm$polarity[j[k]], peaks[[k]])
            selected <- !length(accessions) | accession %in% accessions
            eligible <<- eligible + sum(selected)
            keep <- which(selected & usable)
            if (!length(keep)) return()
            scores <- .sm_fragment_scores(qp, peaks[keep], param)
            j <- j[keep]
            parts[[length(parts) + 1L]] <<- data.frame(queryIndex = i, referenceIndex = j,
                msLevel = level, query_run_id = qm$run[i], query_native_id = qm$native[i],
                query_spectrum_id_ = qm$id[i], reference_run_id = lm$run[j],
                reference_native_id = lm$native[j], reference_spectrum_id_ = lm$id[j],
                library_accession = if (ref$role == "library") accession[keep] else NA_character_,
                reason = decision[keep], scores, row.names = NULL)
        })
        hits <- do.call(rbind, parts)
        scored <- NROW(hits)
        if (scored) {
            hits$rank <- .sm_rank_scores(hits[[scoring$rankScore]])
            if (!is.na(a$topn))
                hits <- hits[which(hits$rank <= a$topn), , drop = FALSE]
        }
        outcome <- if (scored) "scored" else if (eligible) "unscorable" else
            "no_eligible_reference"
        reason <- if (scored) {
            if (NROW(hits)) "hit_retained" else "rank_cutoff"
        } else if (eligible) "no_usable_peaks" else if (length(unique(reasons)) == 1L)
            reasons[1L] else "no_compatible_reference"
        list(hits = hits, status = data.frame(status = outcome, reason = reason,
            nEligible = eligible, nScored = scored, nRetained = NROW(hits)))
    }
    results <- lapply(qrows, match_query)
    if (length(results)) {
        counts <- do.call(rbind, lapply(results, `[[`, "status"))
        status[qrows, names(counts)] <- counts
    }
    hits <- do.call(rbind, lapply(results, `[[`, "hits"))
    if (!NROW(hits)) hits <- NULL
    sources <- .sm_hit_sources(hits, qsource, lsource)
    list(matchedResults = sources$hits, query = query, reference = reference,
        referenceRole = ref$role, queryStatus = status, sourceMeasurements = sources$table,
        scoring = scoring)
}
