# Keep spectra in their original backends. A results dataset can contain
# consensus rows of its own and measured rows in several external studies.
.blink_results_blocks <- function(path, manifest, types, filter, sourcePaths) {
    reader <- .mzs_reader(path, sourcePaths, manifest)
    blocks <- list()
    methods <- if (!all(is.na(types))) sub("^av_", "", setdiff(types, "scan"))
    stored <- .mzs_processed_selection(reader, NULL, 2:4, methods, NULL)
    if (length(stored$spectra)) {
        sp <- stored$spectra
        info <- stored$info
        peaks <- function(rows) {
            flags <- if (filter) paste0("x_mspurity_", c("snr", "minnum", "minfrac", "ra"), "_pass_flag")
            lapply(.mzs_matching_peaks(sp[rows], flags), function(p) {
                if (filter) {
                    have <- intersect(flags, names(p))
                    keep <- if (length(have)) Reduce(`&`, lapply(p[have], `%in%`, TRUE)) else FALSE
                    p <- p[keep, , drop = FALSE]
                }
                .sm_clean_peaks(cbind(mz = p$mz, intensity = p$i))
            })
        }
        blocks$self <- list(sp = sp, peaks = peaks, path = path,
            meta = data.frame(query_source = "self", query_run_id = info$run_id,
                query_spectrum_id_ = info$spectrum_id_, query_native_id = info$native_id,
                feature_id_ = info$feature_id_, grpid = info$grpid, rt = info$retention_time,
                precursor_mz = info$precursor_mz,
                purity = info$x_mspurity_precursor_purity))
    }
    if (!"scan" %in% types) return(blocks)
    scans <- reader$table("x_mspurity_scan")
    refs <- .mzs_read_scan_refs(scans)
    raw <- .mzs_raw_spectra(reader, refs, nativeOnly = TRUE)
    links <- reader$table("spectrum_feature")
    feature <- links$feature_id_[match(.mzs_ref_key(refs), .mzs_ref_key(.mzs_read_scan_refs(links)))]
    msn <- reader$table("x_mspurity_msn_feature")
    missing <- is.na(feature)
    feature[missing] <- msn$feature_id_[match(scans$scan_annotation_id_[missing], msn$scan_annotation_id_)]
    groups <- .mzs_original_ids(reader, "feature", "grpid")
    blocks <- c(blocks, lapply(names(raw), function(key) {
        sp <- raw[[key]]
        rows <- which(refs$source == key)
        rows <- rows[match(sp$spectrum_id_, refs$spectrum_id_[rows])]
        s <- scans[rows, , drop = FALSE]
        peaks <- function(i) {
            if (!filter) return(lapply(.msp_peaks(sp[i]), .sm_clean_peaks))
            ann <- reader$table("x_mspurity_scan_peak", "scan_annotation_id_", s$scan_annotation_id_[i])
            lapply(s$scan_annotation_id_[i], function(id) {
                p <- ann[ann$scan_annotation_id_ == id & ann$purity_pass_flag %in% TRUE &
                    ann$intensity_pass_flag %in% TRUE & ann$ra_pass_flag %in% TRUE &
                    ann$snr_pass_flag %in% TRUE, , drop = FALSE]
                p <- p[order(p$peak_rank), , drop = FALSE]
                .sm_clean_peaks(cbind(mz = p$mz, intensity = p$intensity))
            })
        }
        fid <- feature[rows]
        list(sp = sp, peaks = peaks, path = unique(Spectra::dataStorage(sp)),
            meta = data.frame(query_source = key, query_run_id = s$ms2_run_id,
                query_spectrum_id_ = s$ms2_spectrum_id_, query_native_id = s$ms2_native_id,
                feature_id_ = fid, grpid = as.numeric(groups[as.character(fid)]),
                rt = s$retention_time, precursor_mz = s$precursor_mz, purity = s$in_purity))
    }))
    blocks
}

.blink_selected <- function(meta, groups, purity) {
    if (!NROW(meta)) return(integer())
    keep <- rep(TRUE, nrow(meta))
    if (!all(is.na(groups))) keep <- meta$grpid %in% groups
    if (!is.na(purity)) keep <- keep & meta$purity > purity
    which(keep)
}

.blink_parquet_search <- function(path, manifest, a, lm, key) {
    queries <- .blink_results_blocks(path, manifest, a$q_spectraTypes, a$q_spectraFilter, a$sourcePaths)
    experimental <- identical(lm$role, "results")
    references <- if (experimental)
        .blink_results_blocks(a$l_dbPth, lm, a$l_spectraTypes, a$l_spectraFilter, a$sourcePaths)
    else list(list(sp = .sm_reference(a$l_dbPth, "library")$spectra))
    rawSources <- list()
    # An empty reference selection must still produce a status for every query.
    if (!length(references)) references <- list(list(sp = Spectra::Spectra()))
    offset <- 0L
    for (k in seq_along(references)) {
        r <- references[[k]]
        r$key <- key
        r$offset <- offset
        offset <- offset + length(r$sp)
        if (experimental && length(r$sp) && !all(r$meta$query_source == "self")) {
            # Reference measured scans belong to their original study dataset.
            sm <- MsBackendParquet::readManifest(r$path)
            known <- manifest
            known$sources <- c(known$sources, list(list(key = key)), rawSources)
            # Older native studies have no UID. Reuse the source key when
            # both manifests resolve to the same study directory.
            samePath <- Filter(function(src) identical(src$role, "study") &&
                identical(.mzs_locate_source(path, src, a$sourcePaths), normalizePath(r$path)), manifest$sources)
            r$key <- if (length(samePath)) samePath[[1L]]$key else .mzs_library_key(known, sm)
            rawSources[[r$key]] <- list(key = r$key, path = r$path, uid = sm$uid,
                runs = union(rawSources[[r$key]]$runs, r$meta$query_run_id))
        }
        filters <- a
        if (experimental) {
            filters$l_purity <- filters$l_rtrange <- NA_real_
        }
        metadata <- .mzs_library_metadata(r$sp, filters, r$key)
        parts <- list()
        .mzs_batches(length(r$sp), .mzs_batch_size(TRUE), function(i)
            parts[[length(parts) + 1L]] <<- metadata(i))
        r$library <- if (length(parts)) do.call(rbind, parts) else metadata(integer())
        r$rows <- r$library$reference_spectrum_id_
        if (experimental) r$rows <- intersect(r$rows,
            .blink_selected(r$meta, a$l_xcmsGroups, a$l_purity))
        references[[k]] <- r
    }
    direct <- a
    direct$updateDb <- direct$copyDb <- FALSE
    direct$referenceRole <- if (experimental) "experimental" else "library"
    direct$topn <- NA_integer_  # rank across all reference sources below
    matches <- attempted <- statuses <- identities <- list()
    offset <- 0L
    for (q in queries) {
        rows <- .blink_selected(q$meta, a$q_xcmsGroups, a$q_purity)
        direct$q_dbPth <- q$sp
        status <- NULL
        for (r in references) {
            direct$l_dbPth <- r$sp
            result <- .sm_search_spectra(direct, "l_pol", rows, r$rows,
                queryPeaks = function(i) q$peaks(i)[[1L]], referencePeaks = r$peaks,
                queryRt = q$meta$rt, referenceRt = r$meta$rt)
            s <- result$queryStatus
            if (is.null(status)) status <- s else {
                for (name in c("nEligible", "nScored")) status[[name]] <- status[[name]] + s[[name]]
                disagree <- status$reason != s$reason
                status$reason[disagree] <- "no_compatible_reference"
            }
            h <- result$matchedResults
            if (is.null(h)) next
            sources <- result$sourceMeasurements
            sources$inputIndex <- sources$inputIndex + ifelse(sources$side == "query", offset, r$offset)
            identities[[length(identities) + 1L]] <- sources
            i <- h$queryIndex
            j <- h$referenceIndex
            lib <- r$library[match(j, r$library$reference_spectrum_id_), , drop = FALSE]
            qmeta <- q$meta[i, , drop = FALSE]
            h$queryIndex <- i + offset
            h$referenceIndex <- j + r$offset
            h$qkey <- as.character(h$queryIndex)
            h$order <- h$referenceIndex
            h$query_source <- qmeta$query_source
            h$query_run_id <- qmeta$query_run_id
            h$query_spectrum_id_ <- qmeta$query_spectrum_id_
            h$query_native_id <- qmeta$query_native_id
            h$feature_id_ <- qmeta$feature_id_
            h$reference_source <- r$key
            h$reference_spectrum_id_ <- if (is(a$l_dbPth, "Spectra")) j else as.numeric(h$reference_spectrum_id_)
            if (experimental) {
                h$reference_run_id <- r$meta$query_run_id[j]
                h$reference_native_id <- r$meta$query_native_id[j]
                h$library_accession <- h$reference_native_id
            } else if (is(a$l_dbPth, "Spectra")) {
                h$reference_native_id <- h$library_accession
                h$reference_run_id <- NA_character_
            }
            h$query_rt <- qmeta$rt
            h$library_rt <- if (experimental) r$meta$rt[j] else lib$retention_time
            h$query_precursor_mz <- qmeta$precursor_mz
            h$library_precursor_mz <- lib$precursor_mz
            h$query_precursor_ion_purity <- qmeta$purity
            h$library_precursor_ion_purity <- if (experimental) r$meta$purity[j] else lib$inPurity
            columns <- c(library_precursor_type = "precursor_type", library_entry_name = "name",
                inchikey = "inchikey", library_source_name = "source_name",
                library_compound_name = "compound_name", library_formula = "formula",
                library_smiles = "smiles", library_inchi = "inchi", library_exact_mass = "exact_mass")
            h[names(columns)] <- lib[unname(columns)]
            matches[[length(matches) + 1L]] <- h
        }
        status$queryIndex <- status$queryIndex + offset
        identity <- c("query_source", "query_run_id", "query_spectrum_id_", "query_native_id")
        status[identity] <- q$meta[identity]
        status$status[status$nEligible > 0L] <- "unscorable"
        status$reason[status$nEligible > 0L] <- "no_usable_peaks"
        status$status[status$nScored > 0L] <- "scored"
        status$reason[status$nScored > 0L] <- "hit_retained"
        statuses[[length(statuses) + 1L]] <- status
        att <- q$meta
        att$qkey <- as.character(status$queryIndex)
        attempted[[length(attempted) + 1L]] <- att[!status$status %in% c("invalid_query", "filtered"), ]
        offset <- offset + nrow(att)
    }
    status <- do.call(rbind, statuses)
    if (is.null(status)) {
        meta <- data.frame(level = integer(), mz = numeric(), run = character(), native = character(), id = numeric())
        status <- .sm_query_status(meta, integer(), character())
        status$query_source <- character()
    }
    hits <- do.call(rbind, matches)
    if (is.null(hits)) {
        hits <- data.frame(qkey = character(), order = integer(),
            blinkCosine = numeric(), blinkMatchingIonCount = numeric())
    } else {
        hits <- hits[order(hits$queryIndex, hits$referenceIndex), , drop = FALSE]
        rownames(hits) <- NULL
    }
    list(matches = hits,
        attempted = do.call(rbind, attempted) %||% data.frame(qkey = character(),
            query_source = character(), query_run_id = character(), query_spectrum_id_ = numeric(),
            query_native_id = character()),
        queryStatus = status,
        sourceMeasurements = unique(do.call(rbind, identities)) %||% .sm_hit_sources(NULL, NULL, NULL)$table,
        rawSources = unname(rawSources), scoring = .sm_scoring(a$param))
}
