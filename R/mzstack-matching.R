# Spectral matching on mzStack datasets, and identification evidence.
#
# Scoring is the same as the SQLite route (queryVlibrarySingle()); library
# metadata is read once and peaks are read only for candidates.
# Every attempted query gets a coverage row, so "not attempted" (no row),
# "nothing matched" (no hits) and "list cut" (truncated, with boundary score)
# can be told apart.

# Scores spectralMatching() computes, with their kinds.
.MZS_SCORE_KINDS_SM <- c(dpc = "raw", rdpc = "raw", cdpc = "raw",
                         mcount = "count", allcount = "count",
                         mpercent = "raw")

# ---------------------------------------------------------------------------
# Compounds.
# ---------------------------------------------------------------------------

#' A full InChIKey, or NA.
#'
#' @noRd
.mzs_inchikey <- function(x) {
    x <- toupper(trimws(as.character(x)))
    ok <- !is.na(x) & grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", x)
    ifelse(ok, x, NA_character_)
}

#' Compound rows for a set of matches, reusing existing compounds with the
#' same full InChIKey.
#'
#' @param cmp data.frame: inchikey, chemical_name, chemical_formula, smiles,
#'     inchi, theoretical_neutral_mass, average_molecular_weight,
#'     compound_class, plus optional xref columns pubchem / chemspider.
#'
#' @param existing the dataset's compound table, or NULL.
#'
#' @return list(compound = new rows, xref = new xref rows, id = function
#'     mapping an InChIKey to its compound_id_).
#'
#' @noRd
.mzs_compound_rows <- function(cmp, existing = NULL) {
    raw <- .mzs_na_empty(cmp$inchikey)
    valid <- .mzs_inchikey(raw)
    ## A partial InChIKey still gets a compound, with a null inchikey; the
    ## raw value is kept as a source identifier.
    cmp$key <- ifelse(!is.na(valid), valid, ifelse(is.na(raw), NA,
                                                    paste0("raw:", raw)))
    cmp$inchikey <- valid
    cmp$raw <- raw
    cmp <- cmp[!is.na(cmp$key) & !duplicated(cmp$key), , drop = FALSE]
    have <- if (!is.null(existing) && nrow(existing))
        data.frame(compound_id_ = existing$compound_id_,
                   key = existing$inchikey, stringsAsFactors = FALSE)
    else data.frame(compound_id_ = numeric(), key = character())
    have <- have[!is.na(have$key), , drop = FALSE]
    new <- cmp[!cmp$key %in% have$key, , drop = FALSE]
    offset <- if (!is.null(existing) && nrow(existing))
        max(existing$compound_id_) else 0
    compound <- NULL
    xref <- NULL
    col <- function(c, f = as.character)
        if (is.null(new[[c]])) rep(f(NA), nrow(new)) else f(new[[c]])
    if (nrow(new)) {
        compound <- data.frame(
            compound_id_ = offset + seq_len(nrow(new)),
            inchikey = new$inchikey,
            inchikey_block1 = substr(new$inchikey, 1L, 14L),
            inchi = col("inchi"), smiles = col("smiles"),
            chemical_formula = col("chemical_formula"),
            chemical_name = col("chemical_name"),
            theoretical_neutral_mass = suppressWarnings(
                col("theoretical_neutral_mass", as.numeric)),
            average_molecular_weight = suppressWarnings(
                col("average_molecular_weight", as.numeric)),
            compound_class = col("compound_class"),
            stringsAsFactors = FALSE)
        dbs <- c(pubchem = "PubChem Compound", chemspider = "ChemSpider")
        xref <- do.call(.mzs_rbind_fill, lapply(names(dbs), function(c) {
            if (is.null(new[[c]]))
                return(NULL)
            v <- .mzs_na_empty(new[[c]])
            keep <- !is.na(v)
            if (!any(keep))
                return(NULL)
            data.frame(compound_id_ = compound$compound_id_[keep],
                       database = dbs[[c]], identifier = v[keep],
                       stringsAsFactors = FALSE)
        }))
        have <- rbind(have, data.frame(compound_id_ = compound$compound_id_,
                                       key = new$key,
                                       stringsAsFactors = FALSE))
    }
    keyof <- function(ik) {
        r <- .mzs_na_empty(ik)
        v <- .mzs_inchikey(r)
        ifelse(!is.na(v), v, ifelse(is.na(r), NA, paste0("raw:", r)))
    }
    list(compound = compound, xref = xref, raw = new$raw,
         id = function(ik) have$compound_id_[match(keyof(ik), have$key)])
}

#' Compound match level: structure, formula or mass only.
#'
#' @noRd
.mzs_match_level <- function(inchikey, formula = NA_character_) {
    ifelse(!is.na(.mzs_inchikey(inchikey)), "exact_structure",
           ifelse(!is.na(.mzs_na_empty(formula)), "formula", "mass_only"))
}

# ---------------------------------------------------------------------------
# Evidence.
# ---------------------------------------------------------------------------

#' Rank each query's candidates by score (ties share a rank) and apply a
#' top-n cut that keeps every tie at the boundary.
#'
#' @return `matches` with `rank` and `keep`, and `coverage` for every
#'     attempted query.
#'
#' @noRd
.mzs_rank_matches <- function(matches, attempted, topn = NA) {
    if (nrow(matches)) {
        matches$rank <- as.integer(stats::ave(
            -matches$dpc, matches$qkey,
            FUN = function(s) rank(s, ties.method = "min", na.last = "keep")))
        matches$keep <- is.na(topn) | (!is.na(matches$rank) &
                                       matches$rank <= topn)
    } else {
        matches$rank <- integer()
        matches$keep <- logical()
    }
    n_cand <- table(factor(matches$qkey, levels = attempted$qkey))
    n_keep <- table(factor(matches$qkey[matches$keep],
                           levels = attempted$qkey))
    boundary <- vapply(attempted$qkey, function(k) {
        s <- matches$dpc[matches$qkey == k & matches$keep]
        if (length(s)) min(s, na.rm = TRUE) else NA_real_
    }, numeric(1))
    attempted$candidates_considered <- as.numeric(n_cand)
    attempted$hits_retained <- as.integer(n_keep)
    attempted$truncated <- as.integer(n_keep) < as.numeric(n_cand)
    attempted$boundary_score <- ifelse(attempted$truncated, boundary,
                                       NA_real_)
    list(matches = matches, coverage = attempted)
}

#' evidence, evidence_score and coverage rows for a set of matches.
#'
#' @param matches data.frame, one row per scored query-library pair, with
#'     qkey, query_* and reference_* columns, feature_id_, compound_id_,
#'     compound_match_level, the six scores and library/query attributes.
#'
#' @param coverage data.frame, one row per attempted query, from
#'     `.mzs_rank_matches()`.
#'
#' @param software_id the software_id_ of the activity.
#'
#' @noRd
.mzs_evidence_rows <- function(matches, coverage, software_id) {
    ev <- matches[matches$keep, , drop = FALSE]
    ev <- ev[order(ev$qkey, ev$rank, ev$order), , drop = FALSE]
    n <- nrow(ev)
    evidence <- data.frame(
        evidence_id_ = seq_len(n),
        evidence_input_id = match(ev$qkey, unique(ev$qkey)),
        query_source = ev$query_source, query_run_id = ev$query_run_id,
        query_spectrum_id_ = ev$query_spectrum_id_,
        query_native_id = ev$query_native_id,
        reference_source = ev$reference_source,
        reference_run_id = ev$reference_run_id,
        reference_spectrum_id_ = ev$reference_spectrum_id_,
        reference_native_id = ev$reference_native_id,
        feature_id_ = ev$feature_id_,
        compound_id_ = ev$compound_id_,
        compound_match_level = ev$compound_match_level,
        identification_method = rep(.MZS_LIBRARY_SEARCH, n),
        ms_level = rep(2L, n),
        adduct_ion = .mzs_na_empty(ev$library_precursor_type),
        exp_mass_to_charge = ev$query_precursor_mz,
        theoretical_mass_to_charge = ev$library_precursor_mz,
        rank = ev$rank,
        database_identifier = ev$library_accession,
        chemical_name = ifelse(is.na(ev$compound_id_),
                               .mzs_na_empty(ev$library_compound_name),
                               NA_character_),
        software_id_ = rep(as.integer(software_id), n),
        x_mspurity_library_entry_name = .mzs_na_empty(ev$library_entry_name),
        x_mspurity_library_source_name = .mzs_na_empty(
            ev$library_source_name),
        x_mspurity_library_rt = as.numeric(ev$library_rt),
        x_mspurity_query_rt = as.numeric(ev$query_rt),
        x_mspurity_rtdiff = as.numeric(ev$library_rt) -
            as.numeric(ev$query_rt),
        x_mspurity_query_precursor_ion_purity = as.numeric(
            ev$query_precursor_ion_purity),
        x_mspurity_library_precursor_ion_purity = as.numeric(
            ev$library_precursor_ion_purity),
        stringsAsFactors = FALSE)
    if (!is.null(ev$x_mspurity_mid))
        evidence$x_mspurity_mid <- ev$x_mspurity_mid
    score <- do.call(rbind, lapply(names(.MZS_SCORE_KINDS_SM), function(t)
        data.frame(evidence_id_ = evidence$evidence_id_,
                   score_term = .MZS_SCORE_TERMS[[t]],
                   score_value = as.numeric(ev[[t]]),
                   score_kind = .MZS_SCORE_KINDS_SM[[t]],
                   stringsAsFactors = FALSE)))
    cov <- data.frame(
        query_source = coverage$query_source,
        query_run_id = coverage$query_run_id,
        query_spectrum_id_ = coverage$query_spectrum_id_,
        query_native_id = coverage$query_native_id,
        candidates_considered = coverage$candidates_considered,
        hits_retained = coverage$hits_retained,
        truncated = coverage$truncated,
        boundary_score = coverage$boundary_score,
        stringsAsFactors = FALSE)
    list(evidence = evidence, evidence_score = score, coverage = cov)
}

#' The selection declaration of an evidence table.
#'
#' @noRd
.mzs_selection <- function(topn) {
    if (is.na(topn))
        return(list(policy = "complete"))
    list(policy = "top_n", n = as.integer(topn),
         score_column = .MZS_SCORE_TERMS[["dpc"]], order = "desc",
         ties = "all")
}

# ---------------------------------------------------------------------------
# Reading spectra to match.
# ---------------------------------------------------------------------------

#' Peaks of a dataset's spectra, read only for `ids`.
#'
#' Read as a Spectra object through MsBackendParquet when it is installed,
#' otherwise with arrow directly.
#'
#' @return named list (by spectrum_id_) of data.frame(mz, i, ...flags).
#'
#' @noRd
.mzs_read_peaks <- function(path, m, ids, extra = character()) {
    if (requireNamespace("MsBackendParquet", quietly = TRUE))
        return(.mzs_read_peaks_spectra(path, ids, extra))
    .mzs_read_peaks_arrow(path, m, ids, extra)
}

#' @noRd
.mzs_read_peaks_spectra <- function(path, ids, extra = character()) {
    ids <- unique(ids[!is.na(ids)])
    if (!length(ids))
        return(list())
    sp <- Spectra::Spectra(Spectra::backendInitialize(
        MsBackendParquet::MsBackendParquet(), path = path))
    k <- match(ids, sp$spectrum_id_)
    sp <- sp[sort(k[!is.na(k)])]
    if (!length(sp))
        return(list())
    cols <- intersect(c("mz", "intensity", extra), Spectra::peaksVariables(sp))
    pk <- Spectra::peaksData(sp, columns = cols)
    out <- lapply(pk, function(x) {
        p <- data.frame(mz = unname(x[, "mz"]), i = unname(x[, "intensity"]))
        for (c in setdiff(cols, c("mz", "intensity")))
            p[[c]] <- as.logical(x[, c])
        p
    })
    names(out) <- as.character(sp$spectrum_id_)
    out
}

#' @noRd
.mzs_read_peaks_arrow <- function(path, m, ids, extra = character()) {
    ids <- unique(ids[!is.na(ids)])
    out <- list()
    if (!length(ids))
        return(out)
    rf <- .mzs_runs_frame(m)
    for (i in seq_len(nrow(rf))) {
        lo <- rf$uid_base[i]
        hi <- lo + rf$n_spectra[i] - 1L
        want <- ids[ids >= lo & ids <= hi]
        if (!length(want))
            next
        ds <- arrow::open_dataset(.mzs_run_dir(path, m$runs[[i]]),
                                  format = "parquet")
        cols <- intersect(c("spectrum_id_", "mz", "intensity", extra),
                          names(ds))
        df <- as.data.frame(dplyr::collect(dplyr::select(dplyr::filter(
            ds, spectrum_id_ %in% want), dplyr::all_of(cols))))
        for (k in seq_len(nrow(df))) {
            p <- data.frame(mz = df$mz[[k]], i = df$intensity[[k]])
            for (c in setdiff(cols, c("spectrum_id_", "mz", "intensity"))) {
                v <- df[[c]][[k]]
                p[[c]] <- if (is.null(v)) rep(NA, nrow(p)) else v
            }
            out[[as.character(df$spectrum_id_[k])]] <- p
        }
    }
    out
}

#' Query spectra of a results dataset: averaged spectra by method and, for
#' "scan", the referenced MS2 scans in a native study.
#'
#' @return list(meta, peaks): meta one row per query (qkey, query_* columns,
#'     feature_id_, grpid, precursor_mz, retention_time, inPurity,
#'     polarity); peaks a list by qkey of data.frame(mz, i).
#'
#' @noRd
.mzs_query_spectra <- function(path, m, a, sourcePaths) {
    tab <- function(n) .mzs_read_table(path, m, n)
    sf <- tab("spectrum_feature")
    ft <- tab("feature")
    si <- tab("source_identifier")
    grp_of <- function(fid) as.numeric(si$source_value[match(
        paste("feature", "grpid", fid),
        paste(si$target_table, si$source_column, si$target_key))])
    types <- a$q_spectraTypes
    methods <- unique(sub("^av_", "", types[types != "scan"]))
    meta <- list()
    peaks <- list()
    rr <- m$results$runs
    runs <- names(rr)[vapply(rr, function(r) r$method %in% methods,
                             logical(1))]
    if (length(runs)) {
        sm <- .mzs_read_spectra_meta(path, m, c(
            "id", "selected_ion_mz", "scan_polarity",
            "x_mspurity_precursor_purity"), runs)
        self <- sf[sf$ms2_source == "self", ]
        fid <- self$feature_id_[match(sm$spectrum_id_,
                                      self$ms2_spectrum_id_)]
        meta[[1]] <- data.frame(
            query_source = "self", query_run_id = sm$run_id,
            query_spectrum_id_ = as.numeric(sm$spectrum_id_),
            query_native_id = sm$id, feature_id_ = fid, grpid = grp_of(fid),
            precursor_mz = sm$selected_ion_mz,
            retention_time = ft$retention_time_in_seconds[match(
                fid, ft$feature_id_)],
            inPurity = sm$x_mspurity_precursor_purity,
            polarity = c("-1" = "negative", "1" = "positive")[
                as.character(sm$scan_polarity)],
            stringsAsFactors = FALSE)
        flags <- c("x_mspurity_snr_pass_flag", "x_mspurity_minnum_pass_flag",
                   "x_mspurity_minfrac_pass_flag", "x_mspurity_ra_pass_flag")
        pk <- .mzs_read_peaks(path, m, sm$spectrum_id_, flags)
        if (isTRUE(a$q_spectraFilter))
            pk <- lapply(pk, function(p) {
                have <- intersect(flags, names(p))
                if (!length(have))
                    return(p[0, c("mz", "i")])
                ok <- Reduce(`&`, lapply(have, function(f) p[[f]] %in% TRUE))
                p[ok, c("mz", "i"), drop = FALSE]
            })
        else
            pk <- lapply(pk, function(p) p[, c("mz", "i"), drop = FALSE])
        names(pk) <- paste0("self:", names(pk))
        peaks <- c(peaks, pk)
    }
    if ("scan" %in% types) {
        sc <- tab("x_mspurity_scan")
        keys <- unique(sc$ms2_source)
        for (k in keys) {
            s <- sc[sc$ms2_source == k, ]
            src <- .mzs_open_source(path, m, k, sourcePaths,
                                    unique(s$ms2_run_id))
            kind <- unique(vapply(src$manifest$runs, function(r) r$kind, ""))
            if (!identical(kind, "native"))
                .mzs_abort("unsupported", "Query scans are in source '", k,
                           "', whose runs are of kind ", paste(kind,
                                                               collapse = ", "),
                           "; spectralMatching(format = \"mzstack\") reads ",
                           "scan peaks from native runs only.")
            lk <- sf[sf$ms2_source == k, ]
            fid <- lk$feature_id_[match(paste(s$ms2_run_id, s$scan_number),
                                        paste(lk$ms2_run_id,
                                              lk$x_mspurity_scan_number))]
            meta[[length(meta) + 1L]] <- data.frame(
                query_source = k, query_run_id = s$ms2_run_id,
                query_spectrum_id_ = s$ms2_spectrum_id_,
                query_native_id = s$ms2_native_id, feature_id_ = fid,
                grpid = grp_of(fid), precursor_mz = s$precursor_mz,
                retention_time = s$retention_time, inPurity = s$in_purity,
                polarity = NA_character_, stringsAsFactors = FALSE)
            pk <- .mzs_read_peaks(src$path, src$manifest, s$ms2_spectrum_id_)
            if (isTRUE(a$q_spectraFilter)) {
                ## Without filterFragSpectra() scan flags no peaks pass,
                ## as in the SQLite route.
                ann <- if (!is.null(m$results$tables$x_mspurity_scan_peak))
                    tab("x_mspurity_scan_peak")
                pk <- stats::setNames(lapply(names(pk), function(id) {
                    if (is.null(ann))
                        return(pk[[id]][0, , drop = FALSE])
                    a_id <- s$scan_annotation_id_[match(as.numeric(id),
                                                        s$ms2_spectrum_id_)]
                    p <- ann[ann$scan_annotation_id_ == a_id, ]
                    p <- p[order(p$peak_rank), ]
                    ok <- p$purity_pass_flag %in% TRUE &
                        p$intensity_pass_flag %in% TRUE &
                        p$ra_pass_flag %in% TRUE & p$snr_pass_flag %in% TRUE
                    data.frame(mz = p$mz[ok], i = p$intensity[ok])
                }), names(pk))
            }
            names(pk) <- paste0(k, ":", names(pk))
            peaks <- c(peaks, pk)
        }
    }
    meta <- do.call(rbind, meta)
    if (is.null(meta))
        return(list(meta = data.frame(), peaks = list()))
    meta$qkey <- paste0(meta$query_source, ":", meta$query_spectrum_id_)
    keep <- rep(TRUE, nrow(meta))
    if (!is.na(a$q_purity))
        keep <- keep & !is.na(meta$inPurity) & meta$inPurity > a$q_purity
    if (!is.na(a$q_pol))
        keep <- keep & !is.na(meta$polarity) &
            tolower(meta$polarity) == tolower(a$q_pol)
    if (!all(is.na(a$q_rtrange)))
        keep <- keep & !is.na(meta$retention_time) &
            meta$retention_time >= a$q_rtrange[1] &
            meta$retention_time <= a$q_rtrange[2]
    if (!all(is.na(a$q_xcmsGroups)))
        keep <- keep & meta$grpid %in% a$q_xcmsGroups
    meta <- meta[keep, , drop = FALSE]
    list(meta = meta, peaks = peaks[meta$qkey])
}

#' Library spectra: a library dataset, or a results dataset's averaged
#' spectra.
#'
#' @noRd
.mzs_library_spectra <- function(lpath, lm, a, key) {
    cols <- c("id", "selected_ion_mz", "scan_polarity", "x_mspurity_name",
              "x_mspurity_precursor_type", "x_mspurity_inchikey",
              "x_mspurity_compound_name", "x_mspurity_formula",
              "x_mspurity_smiles", "x_mspurity_inchi",
              "x_mspurity_exact_mass", "x_mspurity_source_name",
              "x_mspurity_instrument_type", "x_mspurity_instrument",
              "x_mspurity_retention_time", "x_mspurity_precursor_purity")
    if (identical(lm$role, "results")) {
        q <- .mzs_query_spectra(lpath, lm, list(
            q_spectraTypes = a$l_spectraTypes %||% c("av_all", "inter"),
            q_spectraFilter = a$l_spectraFilter, q_purity = a$l_purity,
            q_pol = a$l_pol, q_rtrange = a$l_rtrange,
            q_xcmsGroups = a$l_xcmsGroups), NULL)
        mt <- q$meta
        meta <- data.frame(
            reference_source = key, reference_run_id = mt$query_run_id,
            reference_spectrum_id_ = mt$query_spectrum_id_,
            reference_native_id = mt$query_native_id,
            precursor_mz = mt$precursor_mz,
            retention_time = mt$retention_time,
            inPurity = mt$inPurity, accession = mt$query_native_id,
            name = NA_character_, precursor_type = NA_character_,
            inchikey = NA_character_, compound_name = NA_character_,
            formula = NA_character_, smiles = NA_character_,
            inchi = NA_character_, exact_mass = NA_real_,
            source_name = NA_character_, stringsAsFactors = FALSE)
        peaks <- q$peaks
        names(peaks) <- as.character(mt$query_spectrum_id_)
        return(list(meta = meta, peaks = function(ids)
            peaks[as.character(ids)]))
    }
    sm <- .mzs_read_spectra_meta(lpath, lm, cols)
    pol <- c("-1" = "negative", "1" = "positive")[as.character(
        sm$scan_polarity)]
    meta <- data.frame(
        reference_source = key, reference_run_id = sm$run_id,
        reference_spectrum_id_ = as.numeric(sm$spectrum_id_),
        reference_native_id = sm$id, precursor_mz = sm$selected_ion_mz,
        retention_time = as.numeric(sm$x_mspurity_retention_time),
        inPurity = as.numeric(sm$x_mspurity_precursor_purity),
        accession = sm$id, name = sm$x_mspurity_name,
        precursor_type = sm$x_mspurity_precursor_type,
        inchikey = sm$x_mspurity_inchikey,
        compound_name = sm$x_mspurity_compound_name,
        formula = sm$x_mspurity_formula, smiles = sm$x_mspurity_smiles,
        inchi = sm$x_mspurity_inchi,
        exact_mass = as.numeric(sm$x_mspurity_exact_mass),
        source_name = sm$x_mspurity_source_name,
        polarity = unname(pol),
        instrument_type = sm$x_mspurity_instrument_type,
        instrument = sm$x_mspurity_instrument, stringsAsFactors = FALSE)
    keep <- rep(TRUE, nrow(meta))
    if (!is.na(a$l_pol))
        keep <- keep & !is.na(meta$polarity) &
            tolower(meta$polarity) == tolower(a$l_pol)
    if (!all(is.na(a$l_accessions)))
        keep <- keep & meta$accession %in% a$l_accessions
    if (!all(is.na(a$l_instrumentTypes)))
        keep <- keep & meta$instrument_type %in% a$l_instrumentTypes
    if (!all(is.na(a$l_instruments)))
        keep <- keep & meta$instrument %in% a$l_instruments
    if (!all(is.na(a$l_sources)))
        keep <- keep & meta$source_name %in% a$l_sources
    if (!all(is.na(a$l_rtrange)))
        keep <- keep & !is.na(meta$retention_time) &
            meta$retention_time >= a$l_rtrange[1] &
            meta$retention_time <= a$l_rtrange[2]
    if (!is.na(a$l_purity))
        keep <- keep & !is.na(meta$inPurity) & meta$inPurity > a$l_purity
    meta <- meta[keep, , drop = FALSE]
    list(meta = meta, peaks = function(ids)
        .mzs_read_peaks(lpath, lm, ids))
}

#' Score every query against the library spectra in its precursor window.
#'
#' Candidates are those whose ppm range overlaps the query's: a coarse m/z
#' prefilter, then the exact SQLite-route rule.
#'
#' @noRd
.mzs_score <- function(query, library, a) {
    qm <- query$meta
    lm <- library$meta
    attempted <- list()
    out <- list()
    for (i in seq_len(nrow(qm))) {
        q <- qm[i, ]
        qp <- query$peaks[[q$qkey]]
        if (is.null(qp) || !nrow(qp))
            next
        attempted[[length(attempted) + 1L]] <- q
        cand <- lm
        if (isTRUE(a$usePrecursors)) {
            e <- a$l_ppmPrec * 1e-6
            lo <- q$precursor_mz - q$precursor_mz * 1e-6 * a$q_ppmPrec
            hi <- q$precursor_mz + q$precursor_mz * 1e-6 * a$q_ppmPrec
            near <- !is.na(cand$precursor_mz) &
                cand$precursor_mz >= lo / (1 + e) * (1 - 1e-12) &
                cand$precursor_mz <= hi / (1 - e) * (1 + 1e-12)
            cand <- cand[near, , drop = FALSE]
            cand <- cand[(hi >= cand$precursor_mz -
                              ((cand$precursor_mz * 0.000001) * a$l_ppmPrec)) &
                         (cand$precursor_mz +
                              ((cand$precursor_mz * 0.000001) * a$l_ppmPrec) >=
                              lo), , drop = FALSE]
        }
        if (!is.na(a$rttol))
            cand <- cand[!is.na(cand$retention_time) &
                         abs(cand$retention_time - q$retention_time) < a$rttol,
                         , drop = FALSE]
        if (!nrow(cand))
            next
        lpk <- library$peaks(cand$reference_spectrum_id_)
        scored <- lapply(seq_len(nrow(cand)), function(j) {
            lp <- lpk[[as.character(cand$reference_spectrum_id_[j])]]
            if (is.null(lp) || !nrow(lp))
                return(NULL)
            l_speaks <- data.frame(pid = 1L, mz = lp$mz, i = lp$i)
            l_meta <- data.frame(pid = 1L,
                                 retention_time = cand$retention_time[j],
                                 accession = cand$accession[j],
                                 precursor_mz = cand$precursor_mz[j],
                                 precursor_type = cand$precursor_type[j],
                                 name = cand$name[j],
                                 inchikey_id = cand$inchikey[j],
                                 inPurity = cand$inPurity[j],
                                 stringsAsFactors = FALSE)
            r <- queryVlibrarySingle(1L, q_speaksi = qp, l_speakmeta = l_meta,
                                     l_speaks = l_speaks,
                                     q_ppmProd = a$q_ppmProd,
                                     l_ppmProd = a$l_ppmProd, raW = a$raW,
                                     mzW = a$mzW)
            as.numeric(r[c("dpc", "rdpc", "cdpc", "mcount", "allcount",
                           "mpercent")])
        })
        ok <- !vapply(scored, is.null, logical(1))
        if (!any(ok))
            next
        sc <- do.call(rbind, scored[ok])
        colnames(sc) <- c("dpc", "rdpc", "cdpc", "mcount", "allcount",
                          "mpercent")
        cj <- cand[ok, , drop = FALSE]
        out[[length(out) + 1L]] <- data.frame(
            qkey = q$qkey, order = seq_len(nrow(cj)),
            query_source = q$query_source, query_run_id = q$query_run_id,
            query_spectrum_id_ = q$query_spectrum_id_,
            query_native_id = q$query_native_id, feature_id_ = q$feature_id_,
            reference_source = cj$reference_source,
            reference_run_id = cj$reference_run_id,
            reference_spectrum_id_ = cj$reference_spectrum_id_,
            reference_native_id = cj$reference_native_id,
            sc, library_rt = cj$retention_time,
            query_rt = q$retention_time,
            library_precursor_mz = cj$precursor_mz,
            query_precursor_mz = q$precursor_mz,
            library_precursor_ion_purity = cj$inPurity,
            query_precursor_ion_purity = q$inPurity,
            library_accession = cj$accession,
            library_precursor_type = cj$precursor_type,
            library_entry_name = cj$name, inchikey = cj$inchikey,
            library_source_name = cj$source_name,
            library_compound_name = cj$compound_name,
            library_formula = cj$formula, library_smiles = cj$smiles,
            library_inchi = cj$inchi, library_exact_mass = cj$exact_mass,
            stringsAsFactors = FALSE)
    }
    att <- do.call(rbind, attempted)
    if (is.null(att))
        att <- qm[0, , drop = FALSE]
    list(matches = do.call(rbind, out) %||% data.frame(
             qkey = character(), dpc = numeric()),
         attempted = att)
}

# ---------------------------------------------------------------------------
# spectralMatching(format = "mzstack")
# ---------------------------------------------------------------------------

# SQLite-only arguments; refused, not ignored, when set.
.MZS_SM_UNMAPPED <- c("q_raThres", "q_instrumentTypes", "q_instruments",
                      "q_sources", "q_pids", "q_accessions", "l_raThres",
                      "l_pids", "q_dbName", "q_dbHost", "q_dbUser",
                      "q_dbPass", "q_dbPort", "l_dbName", "l_dbHost",
                      "l_dbUser", "l_dbPass", "l_dbPort")

#' A source key for a library, unique in the manifest.
#'
#' @noRd
.mzs_library_key <- function(m, lm) {
    base <- gsub("[^A-Za-z0-9._-]+", "_", lm$library$version %||%
                     if (identical(lm$role, "results")) "query_dataset"
                     else "library")
    base <- substr(sub("^_+", "", base), 1L, 100L)
    if (!nzchar(base))
        base <- "library"
    keys <- vapply(m$sources %||% list(), function(s) s$key, "")
    for (s in m$sources %||% list())
        if (!is.null(lm$uid) && identical(s$uid, lm$uid))
            return(s$key)
    key <- base
    i <- 1L
    while (key %in% keys) {
        i <- i + 1L
        key <- paste0(base, "_", i)
    }
    key
}

#' Copy a results dataset to a new location with a new uid.
#'
#' @noRd
.mzs_fork <- function(from, to, overwrite) {
    .mzs_check_destination(to, overwrite)
    stage <- .mzs_stage(to)
    on.exit(unlink(stage, recursive = TRUE), add = TRUE)
    file.copy(list.files(from, full.names = TRUE, all.files = TRUE,
                         no.. = TRUE), stage, recursive = TRUE)
    unlink(list.files(stage, pattern = "\\.lock$", full.names = TRUE),
           recursive = TRUE)
    m <- .mzs_read_manifest(stage)
    old <- m$uid
    m$uid <- .mzs_uid()
    ## Re-resolve relative source paths against the new location.
    for (i in seq_along(m$sources)) {
        p <- m$sources[[i]]$path
        if (!is.null(p) && !.mzs_is_absolute(p) &&
            .mzs_is_dataset(file.path(from, p)))
            m$sources[[i]]$path <- .mzs_relative_path(
                normalizePath(file.path(from, p)), to)
    }
    .mzs_write_manifest(stage, m)
    .mzs_finalise(stage, to, overwrite)
    old
}

#' spectralMatching(format = "mzstack").
#'
#' @noRd
.spectralMatching_mzstack <- function(a) {
    .mzs_require("spectralMatching(format = \"mzstack\")")
    started <- .mzs_now()
    set <- .MZS_SM_UNMAPPED[vapply(.MZS_SM_UNMAPPED, function(n)
        !all(is.na(a[[n]])), logical(1))]
    if (length(set))
        .mzs_abort("unsupported", "spectralMatching(format = \"mzstack\") ",
                   "does not apply ", paste(set, collapse = ", "), ".",
                   data = list(constructs = set))
    if (!identical(a$q_dbType, "sqlite") || !identical(a$l_dbType, "sqlite"))
        .mzs_abort("unsupported", "q_dbType and l_dbType do not apply to ",
                   "format = \"mzstack\".")
    if (length(a$l_dbPth) != 1L || is.na(a$l_dbPth))
        .mzs_abort("unsupported", "spectralMatching(format = \"mzstack\") ",
                   "needs an mzStack library dataset in 'l_dbPth'; the ",
                   "default library is a SQLite database. Convert it with ",
                   "convertLibraryToMzstack().")
    if (!.mzs_is_dataset(a$q_dbPth))
        .mzs_abort("format", "'", a$q_dbPth, "' is not an mzStack dataset.")
    q <- normalizePath(a$q_dbPth)
    m <- .mzs_read_manifest(q)
    if (!identical(m$role, "results"))
        .mzs_abort("unsupported", "'", q, "' is not a results dataset.")
    lpath <- normalizePath(a$l_dbPth, mustWork = TRUE)
    lm <- .mzs_read_manifest(lpath)
    if (!m$role %in% "results" || !lm$role %in% c("library", "results"))
        .mzs_abort("unsupported", "'", lpath, "' is neither a library nor ",
                   "a results dataset.")
    if (!identical(lm$role, "results") &&
        any(!is.na(c(a$l_spectraTypes, a$l_xcmsGroups))))
        .mzs_abort("unsupported", "l_spectraTypes and l_xcmsGroups apply ",
                   "only when the library is a results dataset.")
    key <- .mzs_library_key(m, lm)

    message("Running msPurity spectral matching (mzStack)")
    query <- .mzs_query_spectra(q, m, a, a$sourcePaths)
    library <- .mzs_library_spectra(lpath, lm, a, key)
    scored <- .mzs_score(query, library, a)
    ranked <- .mzs_rank_matches(scored$matches, scored$attempted, a$topn)
    matches <- ranked$matches

    existing_cmp <- if (!is.null(m$results$tables$compound))
        .mzs_read_table(q, m, "compound")
    cmp <- .mzs_compound_rows(data.frame(
        inchikey = matches$inchikey %||% character(),
        chemical_name = matches$library_compound_name %||% character(),
        chemical_formula = matches$library_formula %||% character(),
        smiles = matches$library_smiles %||% character(),
        inchi = matches$library_inchi %||% character(),
        theoretical_neutral_mass = matches$library_exact_mass %||% numeric(),
        stringsAsFactors = FALSE), existing_cmp)
    if (nrow(matches)) {
        matches$compound_id_ <- cmp$id(matches$inchikey)
        matches$compound_match_level <- .mzs_match_level(
            matches$inchikey, matches$library_formula)
    }
    software_id <- 1L + if (!is.null(m$results$tables$software))
        as.integer(m$results$tables$software$rows) else 0L
    rows <- .mzs_evidence_rows(matches, ranked$coverage, software_id)
    params <- a[setdiff(names(a), c("q_dbPth", "l_dbPth", "sourcePaths"))]

    target <- q
    if (isTRUE(a$updateDb)) {
        forked <- NULL
        if (isTRUE(a$copyDb)) {
            forked <- .mzs_fork(q, a$outPth, overwrite = FALSE)
            target <- normalizePath(a$outPth)
        }
        append <- list(evidence = rows$evidence,
                       evidence_score = rows$evidence_score,
                       coverage = rows$coverage,
                       software = data.frame(
                           software_id_ = 1L, name = "msPurity",
                           version = .mzs_pkg_version("msPurity"),
                           parameters = as.character(.mzs_to_json(
                               .mzs_json_safe(params), pretty = FALSE)),
                           stringsAsFactors = FALSE))
        if (!is.null(cmp$compound))
            append$compound <- cmp$compound
        if (!is.null(cmp$xref))
            append$compound_xref <- cmp$xref
        lsrc <- .mzs_dataset_source(key, lpath, target,
                                    role = lm$role,
                                    run_ids = unique(matches$reference_run_id))
        .mzs_commit_update(
            target, append = append,
            selection = list(evidence = .mzs_selection(a$topn)),
            sources = list(lsrc),
            activity = list(action = "annotate", started = started,
                            fn = "spectralMatching",
                            inputs = list(list(source = key)),
                            parameters = c(params, list(
                                forked_from = forked,
                                candidate_space = paste(
                                    "Library spectra whose precursor ppm",
                                    "window (l_ppmPrec) overlaps the query's",
                                    "(q_ppmPrec), after the l_* filters"),
                                rank_score = "dpc"))))
    }

    ## Same shape as the SQLite route's return value.
    ev <- matches[matches$keep, , drop = FALSE]
    mr <- if (nrow(ev)) data.frame(
        qpid = ev$query_spectrum_id_, query_source = ev$query_source,
        query_run_id = ev$query_run_id, mid = seq_len(nrow(ev)),
        dpc = ev$dpc, rdpc = ev$rdpc, cdpc = ev$cdpc, mcount = ev$mcount,
        allcount = ev$allcount, mpercent = ev$mpercent,
        lpid = ev$reference_spectrum_id_, library_rt = ev$library_rt,
        query_rt = ev$query_rt, rtdiff = ev$library_rt - ev$query_rt,
        library_precursor_mz = ev$library_precursor_mz,
        query_precursor_mz = ev$query_precursor_mz,
        library_precursor_ion_purity = ev$library_precursor_ion_purity,
        query_precursor_ion_purity = ev$query_precursor_ion_purity,
        library_accession = ev$library_accession,
        library_precursor_type = ev$library_precursor_type,
        library_entry_name = ev$library_entry_name, inchikey = ev$inchikey,
        library_source_name = ev$library_source_name,
        library_compound_name = ev$library_compound_name,
        rank = ev$rank, feature_id_ = ev$feature_id_,
        stringsAsFactors = FALSE)
    else NULL
    xr <- NULL
    if (!is.null(mr)) {
        ft <- .mzs_read_table(q, m, "feature")
        f <- ft[match(mr$feature_id_, ft$feature_id_), c(
            "feature_id_", "feature_name", "exp_mass_to_charge",
            "retention_time_in_seconds")]
        xr <- cbind(f, mr[, setdiff(names(mr), "feature_id_")])
        xr <- xr[!is.na(xr$feature_id_), , drop = FALSE]
        xr <- xr[order(xr$feature_id_, -xr$dpc), , drop = FALSE]
        rownames(xr) <- NULL
    }
    list(q_dbPth = target, matchedResults = mr, xcmsMatchedResults = xr)
}
