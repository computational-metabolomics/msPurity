# Spectral libraries as Parquet library datasets: native runs with role
# "library", one run per source library file.
#
# Sources: an msp2db SQLite library (one run per library source), or MSP
# files in a caller-declared dialect. Unknown MSP keys are refused, not
# skipped.
#
# Spectra in a run are sorted by polarity and precursor m/z so a precursor
# window reads few row groups. Peaks keep source order and exact values.

# Library spectrum metadata columns and their storage types.
.MZS_LIBRARY_TYPES <- c(
    ms_level = "i32", time = "f64", scan_polarity = "i32",
    spectrum_representation = "str", id = "str", data_origin = "str",
    selected_ion_mz = "f64", collision_energy = "f64",
    x_mspurity_name = "str", x_mspurity_synonyms = "str",
    x_mspurity_precursor_type = "str", x_mspurity_instrument_type = "str",
    x_mspurity_instrument = "str", x_mspurity_fragmentation_type = "str",
    x_mspurity_resolution = "str", x_mspurity_collision_energy_text = "str",
    x_mspurity_copyright = "str", x_mspurity_column = "str",
    x_mspurity_origin = "str", x_mspurity_splash = "str",
    x_mspurity_mass_accuracy = "f64", x_mspurity_mass_error = "f64",
    x_mspurity_retention_index = "f64", x_mspurity_retention_time = "f64",
    x_mspurity_retention_time_text = "str", x_mspurity_inchikey = "str",
    x_mspurity_compound_name = "str", x_mspurity_formula = "str",
    x_mspurity_smiles = "str", x_mspurity_inchi = "str",
    x_mspurity_exact_mass = "f64", x_mspurity_molecular_weight = "f64",
    x_mspurity_compound_class = "str", x_mspurity_pubchem = "str",
    x_mspurity_chemspider = "str", x_mspurity_other_names = "str",
    x_mspurity_source_name = "str", x_mspurity_source_row = "i32",
    x_mspurity_record_title = "str", x_mspurity_comment = "str",
    x_mspurity_annotation = "str", x_mspurity_grpid = "i32")

#' Value of the first of `keys` present in a named vector, or NA.
#'
#' @noRd
.mzs_lookup <- function(x, keys) {
    for (k in keys)
        if (!is.null(x) && k %in% names(x))
            return(unname(as.character(x[[k]])))
    NA_character_
}

#' Numeric collision energy for `35` or `35 eV`; anything else gives NA.
#'
#' @noRd
.mzs_collision_energy <- function(x) {
    x <- as.character(x)
    ok <- !is.na(x) & grepl("^\\s*[0-9]+(\\.[0-9]+)?\\s*(eV)?\\s*$", x,
                            ignore.case = TRUE)
    out <- rep(NA_real_, length(x))
    out[ok] <- as.numeric(sub("^\\s*([0-9.]+).*$", "\\1", x[ok]))
    out
}

#' Polarity text as +1 / -1 / NA.
#'
#' @noRd
.mzs_library_polarity <- function(x) {
    x <- tolower(trimws(as.character(x)))
    out <- rep(NA_integer_, length(x))
    out[x %in% c("positive", "pos", "p", "+", "1")] <- 1L
    out[x %in% c("negative", "neg", "n", "-", "-1")] <- -1L
    out
}

#' @noRd
.mzs_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

#' Write a library dataset.
#'
#' @param spectra data.frame: `run`, `order_key` (source order) and columns
#'     of `.MZS_LIBRARY_TYPES`.
#'
#' @param peaks list(mz, intensity): each a list with one numeric vector per
#'     spectra row.
#'
#' @param runs data.frame: run (as in `spectra$run`), run_id, ordinal.
#'
#' @noRd
.mzs_write_library <- function(path, spectra, peaks, runs, library,
                               rtUnit, overwrite, activity) {
    if (!is.na(rtUnit)) {
        f <- switch(rtUnit, s = , sec = , second = , seconds = 1 / 60,
                    min = , minute = , minutes = 1,
                    stop("rtUnit must be \"s\" or \"min\".", call. = FALSE))
        spectra$time <- spectra$x_mspurity_retention_time * f
    } else {
        spectra$time <- NA_real_
    }
    out <- list()
    for (i in order(runs$ordinal)) {
        sel <- which(spectra$run == runs$run[i])
        sel <- sel[order(spectra$scan_polarity[sel], spectra$selected_ion_mz[sel],
                         spectra$order_key[sel], na.last = TRUE)]
        meta <- spectra[sel, intersect(names(.MZS_LIBRARY_TYPES),
                                       names(spectra)), drop = FALSE]
        meta <- meta[, vapply(meta, function(v) !all(is.na(v)), logical(1)) |
                         names(meta) %in% c("ms_level", "scan_polarity", "id",
                                            "data_origin", "selected_ion_mz"),
                     drop = FALSE]
        rownames(meta) <- NULL
        out[[length(out) + 1L]] <- list(
            run_id = runs$run_id[i], meta = meta,
            peaks = list(mz = .mzs_peak_list(peaks$mz[sel]),
                         intensity = .mzs_peak_list(peaks$intensity[sel])))
    }
    .mzs_commit_new(path, list(
        role = "library", library = library, runs = out,
        run_types = .MZS_LIBRARY_TYPES, tables = list(),
        activity = activity), overwrite)
}

# ---------------------------------------------------------------------------
# msp2db SQLite libraries.
# ---------------------------------------------------------------------------

.MZS_MSP2DB_TABLES <- c("library_spectra_meta", "library_spectra",
                        "library_spectra_source", "metab_compound")

#' @noRd
.mzs_from_msp2db <- function(db, path, version, rtUnit, runIds, overwrite,
                             started) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
    on.exit(DBI::dbDisconnect(con))
    have <- DBI::dbListTables(con)
    miss <- setdiff(.MZS_MSP2DB_TABLES, have)
    if (length(miss))
        .mzs_abort("format", "'", db, "' is not an msp2db library: it has no ",
                   "table(s) ", paste(miss, collapse = ", "), ".")
    q <- function(sql) suppressWarnings(DBI::dbGetQuery(con, sql))
    extra <- setdiff(have, .MZS_MSP2DB_TABLES)
    n_extra <- vapply(extra, function(t)
        q(sprintf('SELECT COUNT(*) AS n FROM "%s"', t))$n, numeric(1))
    if (any(n_extra > 0))
        .mzs_abort("unsupported", "'", db, "' holds table(s) this converter ",
                   "does not map: ", paste(extra[n_extra > 0],
                                           collapse = ", "), " (R-081).",
                   data = list(constructs = extra[n_extra > 0]))
    src <- q("SELECT * FROM library_spectra_source ORDER BY id")
    meta <- q("SELECT * FROM library_spectra_meta ORDER BY id")
    cmp <- q("SELECT * FROM metab_compound")
    n_other <- q("SELECT COUNT(*) AS n FROM library_spectra WHERE
                  TRIM(CAST(other AS TEXT),
                       ' ' || char(9) || char(10) || char(13)) <> ''")$n
    if (n_other > 0)
        .mzs_abort("unsupported", "library_spectra.other holds peak ",
                   "annotations this converter does not map (R-081).")
    pk <- q("SELECT mz, i, library_spectra_meta_id FROM library_spectra
             ORDER BY id")
    checksum <- .mzs_sha256(db)
    runs <- data.frame(run = src$id, name = src$name,
                       ordinal = seq_len(nrow(src)), stringsAsFactors = FALSE)
    runs$run_id <- vapply(seq_len(nrow(runs)), function(i) {
        r <- .mzs_lookup(runIds, runs$name[i])
        if (!is.na(r)) r else .mzs_mint_run_id(checksum, runs$ordinal[i])
    }, character(1))
    k <- match(meta$inchikey_id, cmp$inchikey_id)
    sname <- src$name[match(meta$library_spectra_source_id, src$id)]
    spectra <- data.frame(
        run = meta$library_spectra_source_id, order_key = meta$id,
        ms_level = as.integer(.mzs_num(meta$ms_level)),
        scan_polarity = .mzs_library_polarity(meta$polarity),
        spectrum_representation = "MS:1000127",
        id = .mzs_na_empty(meta$accession),
        data_origin = paste0("msp2db:", basename(db), ":", sname),
        selected_ion_mz = .mzs_num(meta$precursor_mz),
        collision_energy = .mzs_collision_energy(meta$collision_energy),
        x_mspurity_name = .mzs_na_empty(meta$name),
        x_mspurity_precursor_type = .mzs_na_empty(meta$precursor_type),
        x_mspurity_instrument_type = .mzs_na_empty(meta$instrument_type),
        x_mspurity_instrument = .mzs_na_empty(meta$instrument),
        x_mspurity_fragmentation_type = .mzs_na_empty(meta$fragmentation_type),
        x_mspurity_resolution = .mzs_na_empty(meta$resolution),
        x_mspurity_collision_energy_text = .mzs_na_empty(meta$collision_energy),
        x_mspurity_copyright = .mzs_na_empty(meta$copyright),
        x_mspurity_column = .mzs_na_empty(meta$column),
        x_mspurity_origin = .mzs_na_empty(meta$origin),
        x_mspurity_splash = .mzs_na_empty(meta$splash),
        x_mspurity_mass_accuracy = .mzs_num(meta$mass_accuracy),
        x_mspurity_mass_error = .mzs_num(meta$mass_error),
        x_mspurity_retention_index = .mzs_num(meta$retention_index),
        x_mspurity_retention_time = .mzs_num(meta$retention_time),
        x_mspurity_inchikey = .mzs_na_empty(meta$inchikey_id),
        x_mspurity_compound_name = .mzs_na_empty(cmp$name[k]),
        x_mspurity_formula = .mzs_na_empty(cmp$molecular_formula[k]),
        x_mspurity_smiles = .mzs_na_empty(cmp$smiles[k]),
        x_mspurity_exact_mass = .mzs_num(cmp$exact_mass[k]),
        x_mspurity_molecular_weight = .mzs_num(cmp$molecular_weight[k]),
        x_mspurity_compound_class = .mzs_na_empty(cmp$compound_class[k]),
        x_mspurity_pubchem = .mzs_na_empty(cmp$pubchem_id[k]),
        x_mspurity_chemspider = .mzs_na_empty(cmp$chemspider_id[k]),
        x_mspurity_other_names = .mzs_na_empty(cmp$other_names[k]),
        x_mspurity_source_name = sname,
        x_mspurity_source_row = as.integer(meta$id),
        stringsAsFactors = FALSE)
    by <- factor(pk$library_spectra_meta_id, levels = meta$id)
    peaks <- list(mz = unname(split(pk$mz, by)),
                  intensity = unname(split(pk$i, by)))
    n_peaks <- nrow(pk)
    pk <- by <- NULL
    version <- version %||% paste0(basename(db), " (",
                                   paste(unique(src$parsing_software),
                                         collapse = ", "), ")")
    .mzs_write_library(
        path, spectra, peaks, runs,
        library = list(version = version, digest = paste0("sha256:",
                                                          checksum)),
        rtUnit = rtUnit, overwrite = overwrite,
        activity = list(
            action = "convert", started = started,
            fn = "convertLibraryToParquet",
            inputs = list(.mzs_file_input(db)),
            parameters = list(
                format = "msp2db", rtUnit = rtUnit,
                runIds = if (length(runIds)) as.list(runIds),
                run_names = as.list(stats::setNames(runs$name, runs$run_id)),
                coverage_manifest = list(
                    library_spectra_meta = "mapped",
                    library_spectra = "mapped",
                    library_spectra_source = "mapped",
                    metab_compound = "mapped_partial",
                    library_spectra_annotation = "not_present_in_source"),
                losses = list(
                    list(construct = "metab_compound.created_at, updated_at",
                         disposition = "no_such_concept_in_target",
                         reason = "Row timestamps of the library database."),
                    list(construct = "metab_compound (unmatched rows)",
                         disposition = "dropped_by_configuration",
                         affected_rows = sum(!cmp$inchikey_id %in%
                                             meta$inchikey_id),
                         reason = "Compounds no library spectrum names."),
                    list(construct = "library_spectra.id",
                         disposition = "dropped_by_configuration",
                         affected_rows = n_peaks,
                         reason = "A positional peak row number."),
                    list(construct = "SPLASH",
                         disposition = "not_present_in_source",
                         reason = paste("Copied where the library has one;",
                                        "not computed (R-093)."))))))
}

# ---------------------------------------------------------------------------
# MSP files.
# ---------------------------------------------------------------------------

# Dialect key -> field. MassBank keys include their subtag
# ("AC$MASS_SPECTROMETRY: ION_MODE").
.MZS_MSP_DIALECTS <- local({
    massbank <- c(
        "ACCESSION" = "id", "RECORD_TITLE" = "x_mspurity_record_title",
        "DATE" = "ignore", "AUTHORS" = "ignore", "LICENSE" = "ignore",
        "PUBLICATION" = "ignore", "PROJECT" = "ignore",
        "COPYRIGHT" = "x_mspurity_copyright", "COMMENT" = "x_mspurity_comment",
        "CH$NAME" = "x_mspurity_name",
        "CH$COMPOUND_CLASS" = "x_mspurity_compound_class",
        "CH$FORMULA" = "x_mspurity_formula",
        "CH$EXACT_MASS" = "x_mspurity_exact_mass",
        "CH$SMILES" = "x_mspurity_smiles", "CH$IUPAC" = "x_mspurity_inchi",
        "CH$LINK: INCHIKEY" = "x_mspurity_inchikey",
        "CH$LINK: PUBCHEM" = "x_mspurity_pubchem",
        "CH$LINK: CHEMSPIDER" = "x_mspurity_chemspider",
        "AC$INSTRUMENT" = "x_mspurity_instrument",
        "AC$INSTRUMENT_TYPE" = "x_mspurity_instrument_type",
        "AC$MASS_SPECTROMETRY: MS_TYPE" = "ms_level",
        "AC$MASS_SPECTROMETRY: ION_MODE" = "scan_polarity",
        "AC$MASS_SPECTROMETRY: COLLISION_ENERGY" =
            "x_mspurity_collision_energy_text",
        "AC$MASS_SPECTROMETRY: FRAGMENTATION_MODE" =
            "x_mspurity_fragmentation_type",
        "AC$MASS_SPECTROMETRY: RESOLUTION" = "x_mspurity_resolution",
        "AC$CHROMATOGRAPHY: RETENTION_TIME" = "x_mspurity_retention_time_text",
        "AC$CHROMATOGRAPHY: COLUMN_NAME" = "x_mspurity_column",
        "MS$FOCUSED_ION: PRECURSOR_M/Z" = "selected_ion_mz",
        "MS$FOCUSED_ION: PRECURSOR_TYPE" = "x_mspurity_precursor_type",
        "PK$SPLASH" = "x_mspurity_splash",
        "PK$NUM_PEAK" = "num_peaks", "PK$PEAK" = "peaks")
    mona <- c(
        "NAME" = "x_mspurity_name", "SYNON" = "x_mspurity_synonyms",
        "SYNONYM" = "x_mspurity_synonyms", "DB#" = "id",
        "INCHIKEY" = "x_mspurity_inchikey", "INCHI" = "x_mspurity_inchi",
        "SMILES" = "x_mspurity_smiles", "FORMULA" = "x_mspurity_formula",
        "MW" = "x_mspurity_molecular_weight",
        "EXACTMASS" = "x_mspurity_exact_mass",
        "EXACT_MASS" = "x_mspurity_exact_mass",
        "PRECURSORMZ" = "selected_ion_mz",
        "PRECURSOR_TYPE" = "x_mspurity_precursor_type",
        "SPECTRUM_TYPE" = "ms_level",
        "INSTRUMENT_TYPE" = "x_mspurity_instrument_type",
        "INSTRUMENT" = "x_mspurity_instrument",
        "ION_MODE" = "scan_polarity",
        "COLLISION_ENERGY" = "x_mspurity_collision_energy_text",
        "RETENTIONTIME" = "x_mspurity_retention_time_text",
        "RETENTION_TIME" = "x_mspurity_retention_time_text",
        "SPLASH" = "x_mspurity_splash", "COMMENTS" = "x_mspurity_comment",
        "COMMENT" = "x_mspurity_comment", "NOTES" = "x_mspurity_comment",
        "NUM PEAKS" = "num_peaks")
    list(massbank = massbank, mona = mona,
         mspurity = c(massbank, "XCMS GROUPID (GRPID)" = "x_mspurity_grpid"))
})

#' Key and value of an MSP line in a dialect.
#'
#' @noRd
.mzs_msp_key <- function(line, dialect) {
    keys <- names(.MZS_MSP_DIALECTS[[dialect]])
    if (dialect == "mona") {
        k <- toupper(trimws(sub(":.*$", "", line)))
        return(list(key = k, value = trimws(sub("^[^:]*:", "", line))))
    }
    k <- toupper(trimws(sub(":.*$", "", line)))
    v <- trimws(sub("^[^:]*:", "", line))
    if (k %in% c("CH$LINK", "AC$MASS_SPECTROMETRY", "AC$CHROMATOGRAPHY",
                 "MS$FOCUSED_ION")) {
        ## Subtag may be followed by ":" or whitespace.
        sub_tag <- toupper(sub(":$", "", sub("\\s.*$", "", v)))
        return(list(key = paste0(k, ": ", sub_tag),
                    value = trimws(sub("^\\S+\\s*", "", v))))
    }
    list(key = k, value = v)
}

#' Parse one MSP file in a declared dialect.
#'
#' @return list of records: `fields`, `mz`, `intensity` (literals, file
#'     order) and `line`.
#'
#' @noRd
.mzs_parse_msp <- function(file, dialect) {
    ## Lines may end in \n, \r\n or \r\r\n; readLines() would split \r\r\n.
    txt <- rawToChar(readBin(file, "raw", file.size(file)))
    Encoding(txt) <- "UTF-8"
    lines <- strsplit(gsub("\r", "\n", gsub("\r+\n", "\n", txt)),
                      "\n", fixed = TRUE)[[1]]
    map <- .MZS_MSP_DIALECTS[[dialect]]
    recs <- list()
    cur <- NULL
    in_peaks <- FALSE
    flush <- function() {
        if (!is.null(cur) && length(cur$fields))
            recs[[length(recs) + 1L]] <<- cur
        cur <<- NULL
        in_peaks <<- FALSE
    }
    for (i in seq_along(lines)) {
        ln <- lines[i]
        if (!nzchar(trimws(ln)) || trimws(ln) == "//") {
            flush()
            next
        }
        if (is.null(cur))
            cur <- list(fields = list(), mz = character(),
                        intensity = character(), line = i)
        if (in_peaks && !grepl("^[A-Za-z][A-Za-z$#_ ]*:", ln)) {
            pairs <- strsplit(trimws(ln), ";")[[1]]
            for (p in pairs) {
                tok <- strsplit(trimws(p), "[[:space:],]+")[[1]]
                tok <- tok[nzchar(tok)]
                if (length(tok) < 2L)
                    next
                cur$mz <- c(cur$mz, tok[1])
                cur$intensity <- c(cur$intensity, tok[2])
            }
            next
        }
        kv <- .mzs_msp_key(ln, dialect)
        field <- if (kv$key %in% names(map)) map[[kv$key]] else NULL
        if (is.null(field))
            .mzs_abort("unsupported", "MSP key '", kv$key, "' at line ", i,
                       " of '", file, "' is not a key of the ", dialect,
                       " dialect. Refusing to skip it (R-081).",
                       data = list(constructs = kv$key, line = i))
        if (field == "peaks" || (dialect == "mona" && field == "num_peaks")) {
            in_peaks <- TRUE
            if (field == "num_peaks")
                cur$fields$num_peaks <- kv$value
            next
        }
        if (field == "ignore")
            next
        old <- cur$fields[[field]]
        cur$fields[[field]] <- if (is.null(old)) kv$value
                               else paste(old, kv$value, sep = "\n")
    }
    flush()
    recs
}

#' MSP records as library spectra.
#'
#' @noRd
.mzs_msp_spectra <- function(recs, run, file) {
    f <- function(name) vapply(recs, function(r)
        as.character(r$fields[[name]] %||% NA_character_)[1L], "")
    n <- length(recs)
    lv <- toupper(trimws(f("ms_level")))
    spectra <- data.frame(
        run = rep(run, n), order_key = seq_len(n),
        ms_level = suppressWarnings(as.integer(sub("^MS", "", lv))),
        scan_polarity = .mzs_library_polarity(f("scan_polarity")),
        spectrum_representation = rep("MS:1000127", n),
        id = .mzs_na_empty(ifelse(is.na(f("id")), f("x_mspurity_record_title"),
                                  f("id"))),
        data_origin = rep(normalizePath(file), n),
        selected_ion_mz = .mzs_num(f("selected_ion_mz")),
        stringsAsFactors = FALSE)
    for (c in setdiff(names(.MZS_LIBRARY_TYPES), names(spectra))) {
        v <- f(c)
        if (all(is.na(v)))
            next
        spectra[[c]] <- if (.MZS_LIBRARY_TYPES[[c]] == "f64") .mzs_num(v)
                        else if (.MZS_LIBRARY_TYPES[[c]] == "i32")
                            suppressWarnings(as.integer(v))
                        else v
    }
    if (!is.null(spectra$x_mspurity_collision_energy_text))
        spectra$collision_energy <- .mzs_collision_energy(
            spectra$x_mspurity_collision_energy_text)
    if (!is.null(spectra$x_mspurity_retention_time_text))
        spectra$x_mspurity_retention_time <- .mzs_num(
            sub("\\s*(s|sec|min)\\s*$", "", spectra$x_mspurity_retention_time_text))
    np <- .mzs_num(f("num_peaks"))
    mz <- lapply(recs, function(r) as.numeric(r$mz))
    intensity <- lapply(recs, function(r) as.numeric(r$intensity))
    bad <- which(!is.na(np) & np != lengths(mz))
    if (length(bad))
        .mzs_abort("semantic", "Record at line ", recs[[bad[1]]]$line, " of '",
                   file, "' declares ", np[bad[1]], " peaks and lists ",
                   length(mz[[bad[1]]]), ".")
    list(spectra = spectra, peaks = list(mz = mz, intensity = intensity))
}

#' @noRd
.mzs_from_msp <- function(files, path, dialect, version, rtUnit, runIds,
                          overwrite, started) {
    if (is.null(dialect) || !dialect %in% names(.MZS_MSP_DIALECTS))
        .mzs_abort("unsupported", "MSP has no specification and its dialects ",
                   "disagree on every key: declare 'dialect' as one of ",
                   paste(names(.MZS_MSP_DIALECTS), collapse = ", "), ".")
    files <- normalizePath(files, mustWork = TRUE)
    sums <- vapply(files, .mzs_sha256, "")
    parts <- lapply(seq_along(files), function(i)
        .mzs_msp_spectra(.mzs_parse_msp(files[i], dialect), i, files[i]))
    spectra <- do.call(.mzs_rbind_fill, lapply(parts, `[[`, "spectra"))
    peaks <- list(mz = do.call(c, lapply(parts, function(p) p$peaks$mz)),
                  intensity = do.call(c, lapply(parts, function(p)
                      p$peaks$intensity)))
    runs <- data.frame(run = seq_along(files), ordinal = seq_along(files),
                       stringsAsFactors = FALSE)
    runs$run_id <- vapply(seq_along(files), function(i) {
        r <- .mzs_lookup(runIds, c(files[i], basename(files[i])))
        if (!is.na(r)) r else .mzs_mint_run_id(sums[i], i)
    }, character(1))
    digest <- if (length(files) == 1L) sums[[1]]
              else .mzs_sha256_text(paste(sums, collapse = "\n"))
    .mzs_write_library(
        path, spectra, peaks, runs,
        library = list(version = version %||% paste(basename(files),
                                                    collapse = ", "),
                       digest = paste0("sha256:", digest)),
        rtUnit = rtUnit, overwrite = overwrite,
        activity = list(
            action = "convert", started = started,
            fn = "convertLibraryToParquet",
            inputs = lapply(files, .mzs_file_input),
            parameters = list(
                format = "msp", dialect = dialect, rtUnit = rtUnit,
                runIds = if (length(runIds)) as.list(runIds),
                coverage_manifest = list(
                    records = "mapped", peaks = "mapped",
                    unknown_keys = "refused"),
                losses = list(list(
                    construct = "SPLASH", disposition = "not_present_in_source",
                    reason = paste("Copied where the file has one; not",
                                   "computed (R-093)."))))))
}

# ---------------------------------------------------------------------------
# Spectra objects (MSP through MsBackendMsp, MassBank, MGF and others).
# ---------------------------------------------------------------------------

# Library field -> spectra variable, for the variables Spectra and its
# common library backends (MsBackendMsp, MsBackendMassbank) use.
.MZS_SPECTRA_LIBRARY_MAP <- c(
    x_mspurity_name = "name", x_mspurity_compound_name = "name",
    x_mspurity_synonyms = "synonym", x_mspurity_precursor_type = "adduct",
    x_mspurity_inchikey = "inchikey", x_mspurity_inchi = "inchi",
    x_mspurity_smiles = "smiles", x_mspurity_formula = "formula",
    x_mspurity_exact_mass = "exactmass",
    x_mspurity_instrument = "instrument",
    x_mspurity_instrument_type = "instrument_type",
    x_mspurity_splash = "splash", x_mspurity_comment = "comment")

# Spectra variables the conversion reads besides the mapped ones.
.MZS_SPECTRA_CORE_READ <- c("msLevel", "polarity", "precursorMz", "rtime",
                            "collisionEnergy", "dataOrigin", "accession",
                            "spectrumId")

#' Library spectra from a Spectra object, one run per data origin.
#'
#' @noRd
.mzs_spectra_library <- function(x, mapping) {
    sv <- Spectra::spectraVariables(x)
    d <- Spectra::spectraData(x, columns = sv)
    n <- length(x)
    col <- function(v) {
        if (!v %in% sv) return(rep(NA, n))
        val <- d[[v]]
        # List variables, such as synonyms, are kept as one string.
        if (is.list(val) || is(val, "List"))
            val <- vapply(as.list(val), function(e)
                if (!length(e) || all(is.na(e))) NA_character_
                else paste(as.character(e), collapse = "; "), "")
        val
    }
    origin <- as.character(col("dataOrigin"))
    origin[is.na(origin)] <- "Spectra"
    runs <- unique(origin)
    id <- as.character(col("accession"))
    if (all(is.na(id)))
        id <- as.character(col("spectrumId"))
    id[is.na(id)] <- paste0("spectrum_", which(is.na(id)))
    pol <- suppressWarnings(as.integer(col("polarity")))
    spectra <- data.frame(
        run = match(origin, runs), order_key = seq_len(n),
        ms_level = as.integer(col("msLevel")),
        scan_polarity = ifelse(pol %in% 1L, 1L,
                               ifelse(pol %in% 0L, -1L, NA_integer_)),
        spectrum_representation = rep("MS:1000127", n),
        id = id, data_origin = origin,
        selected_ion_mz = as.numeric(col("precursorMz")),
        collision_energy = .mzs_num(col("collisionEnergy")),
        x_mspurity_retention_time = as.numeric(col("rtime")),
        stringsAsFactors = FALSE)
    for (f in names(mapping)) {
        v <- col(mapping[[f]])
        if (all(is.na(v)))
            next
        spectra[[f]] <- switch(.MZS_LIBRARY_TYPES[[f]],
                               f64 = .mzs_num(v),
                               i32 = suppressWarnings(as.integer(v)),
                               as.character(v))
    }
    pk <- Spectra::peaksData(x, columns = c("mz", "intensity"))
    list(spectra = spectra,
         peaks = list(mz = lapply(pk, function(p) unname(p[, "mz"])),
                      intensity = lapply(pk, function(p)
                          unname(p[, "intensity"]))),
         runs = runs,
         unmapped = setdiff(sv, c(.MZS_SPECTRA_CORE_READ, mapping,
                                  Spectra::coreSpectraVariables())))
}

#' @noRd
.mzs_from_spectra <- function(x, path, mapping, version, rtUnit, runIds,
                              overwrite, started) {
    bad <- setdiff(names(mapping), names(.MZS_LIBRARY_TYPES))
    if (length(bad))
        stop("'mapping' names fields that a library does not have: ",
             paste(bad, collapse = ", "), ".", call. = FALSE)
    lib <- .mzs_spectra_library(x, mapping)
    runs <- data.frame(run = seq_along(lib$runs),
                       ordinal = seq_along(lib$runs),
                       stringsAsFactors = FALSE)
    num <- function(v) format(v, digits = 17)
    digest <- .mzs_sha256_text(paste(c(
        lib$spectra$id, num(lib$spectra$selected_ion_mz),
        vapply(seq_along(lib$peaks$mz), function(i) paste(
            num(lib$peaks$mz[[i]]), num(lib$peaks$intensity[[i]]),
            collapse = " "), "")), collapse = "\n"))
    runs$run_id <- vapply(seq_along(lib$runs), function(i) {
        r <- .mzs_lookup(runIds, lib$runs[i])
        if (!is.na(r)) r else .mzs_mint_run_id(
            .mzs_sha256_text(paste(digest, lib$runs[i])), i)
    }, character(1))
    .mzs_write_library(
        path, lib$spectra, lib$peaks, runs,
        library = list(version = version %||% paste(basename(lib$runs),
                                                    collapse = ", "),
                       digest = paste0("sha256:", digest)),
        rtUnit = rtUnit, overwrite = overwrite,
        activity = list(
            action = "convert", started = started,
            fn = "convertLibraryToParquet", inputs = list(),
            parameters = list(
                format = "spectra", rtUnit = rtUnit,
                mapping = as.list(mapping),
                runIds = if (length(runIds)) as.list(runIds),
                coverage_manifest = list(
                    spectra = "mapped", peaks = "mapped",
                    unmapped_variables = "dropped"),
                losses = lapply(lib$unmapped, function(v) list(
                    construct = v, disposition = "dropped_by_configuration",
                    reason = "A spectra variable with no library field.")))))
}

#' Convert a spectral library to a Parquet library dataset
#'
#' @description
#'
#' Writes a spectral library as a Parquet library dataset (one native run
#' per source library file) for use with
#' `spectralMatching(format = "parquet")`.
#'
#'  * `format = "msp2db"`: an msp2db SQLite library, such as the default
#'    library of `spectralMatching()`; one run per library source.
#'  * `format = "msp"`: MSP files, one run per file. `dialect` must be
#'    declared: `"massbank"`, `"mona"` (NIST/MoNA) or `"mspurity"` (files
#'    written by [createMSP()]). Unknown keys are an error.
#'  * `format = "spectra"`, the default when `x` is a `Spectra` object: a
#'    library read with any Spectra backend, such as MsBackendMsp,
#'    MsBackendMassbank or MsBackendMgf; one run per data origin. Spectra
#'    variables are mapped to library fields by `mapping`, and those with no
#'    field are dropped and recorded as such in the dataset's provenance.
#'    Retention times are taken as seconds, the Spectra convention, unless
#'    `rtUnit` says otherwise.
#'
#' Peaks keep their source order and exact values. Collision energy is kept
#' as written, with a numeric value only for plain numbers (`35`, `35 eV`).
#' Retention time is converted to minutes only when `rtUnit` is given.
#' Licence and attribution fields are kept.
#'
#' Requires the suggested packages arrow and jsonlite.
#'
#' @param x `character`: path of the msp2db database, or of the MSP files;
#'     or a `Spectra` object.
#'
#' @param path `character(1)`, the library dataset to create.
#'
#' @param format `"msp2db"`, `"msp"` or `"spectra"`.
#'
#' @param dialect `character(1)`, required for MSP: `"massbank"`, `"mona"`
#'     or `"mspurity"`.
#'
#' @param version `character(1)`, the library release; derived from the
#'     source by default.
#'
#' @param rtUnit `character(1)`, retention time unit of the source, `"s"`
#'     or `"min"`, or `NA` if unknown.
#'
#' @param runIds named `character`: library source name (msp2db) or file
#'     (MSP) -> run id. By default minted from checksum and ordinal.
#'
#' @param overwrite `logical(1)`, whether to replace an existing dataset.
#'
#' @param mapping named `character`, for `format = "spectra"`: library field
#'     (such as `x_mspurity_inchikey`) -> spectra variable. `NULL`, the
#'     default, maps the variable names used by MsBackendMsp and
#'     MsBackendMassbank: name, synonym, adduct, inchikey, inchi, smiles,
#'     formula, exactmass, instrument, instrument_type, splash and comment.
#'
#' @return The path of the dataset, invisibly.
#'
#' @seealso [spectralMatching()], [validateParquet()]
#'
#' @examples
#' if (requireNamespace("arrow", quietly = TRUE) &&
#'     requireNamespace("jsonlite", quietly = TRUE)) {
#'     msp <- system.file("extdata", "tests", "msp", "av_all.msp",
#'                        package = "msPurity")
#'     out <- file.path(tempdir(), "library.parquet")
#'     convertLibraryToParquet(msp, out, format = "msp", dialect = "mspurity",
#'                             overwrite = TRUE)
#'     validateParquet(out)
#'
#'     if (requireNamespace("MsBackendMsp", quietly = TRUE)) {
#'         mona <- system.file("extdata", "tests", "library",
#'                             "mini_mona.msp", package = "msPurity")
#'         sp <- Spectra::Spectra(mona, source = MsBackendMsp::MsBackendMsp(),
#'                                mapping = c(name = "Name", accession = "DB#",
#'                                            precursorMz = "PrecursorMZ",
#'                                            adduct = "Precursor_type",
#'                                            inchikey = "InChIKey",
#'                                            formula = "Formula",
#'                                            polarity = "Ion_mode"))
#'         out <- file.path(tempdir(), "library-spectra.parquet")
#'         convertLibraryToParquet(sp, out, overwrite = TRUE)
#'     }
#' }
#' @export
convertLibraryToParquet <- function(x, path,
                                    format = c("msp2db", "msp", "spectra"),
                                    dialect = NULL, version = NULL,
                                    rtUnit = NA_character_, runIds = NULL,
                                    overwrite = FALSE, mapping = NULL) {
    .mzs_require("convertLibraryToParquet()")
    started <- .mzs_now()
    if (is(x, "Spectra")) {
        if (!missing(format) && !identical(format, "spectra"))
            stop("'x' is a Spectra object; use format = \"spectra\".",
                 call. = FALSE)
        .mzs_check_destination(path, overwrite)
        if (is.na(rtUnit))
            rtUnit <- "s"
        .mzs_from_spectra(x, path, mapping %||% .MZS_SPECTRA_LIBRARY_MAP,
                          version, rtUnit, runIds, overwrite, started)
        return(invisible(normalizePath(path)))
    }
    format <- match.arg(format)
    if (format == "spectra")
        stop("format = \"spectra\" needs a Spectra object as 'x'.",
             call. = FALSE)
    if (!length(x) || !all(file.exists(x)))
        stop("'x' must name existing file(s).", call. = FALSE)
    .mzs_check_destination(path, overwrite)
    if (format == "msp2db") {
        if (length(x) != 1L)
            stop("An msp2db library is a single database.", call. = FALSE)
        .mzs_from_msp2db(x, path, version, rtUnit, runIds, overwrite, started)
    } else {
        .mzs_from_msp(x, path, dialect, version, rtUnit, runIds, overwrite,
                      started)
    }
    invisible(normalizePath(path))
}

#' Write a library dataset as MSP (MassBank or MoNA dialect), for
#' round-trip tests. Numbers use their shortest exact form.
#'
#' @noRd
.mzs_write_msp <- function(path, file, dialect = c("massbank", "mona")) {
    dialect <- match.arg(dialect)
    m <- .mzs_read_manifest(path)
    num <- function(v) ifelse(is.na(v), NA_character_, .mzs_shortest_double(v))
    out <- character()
    for (r in m$runs) {
        sp <- as.data.frame(arrow::read_parquet(file.path(
            .mzs_run_dir(path, r), "part-0.parquet")))
        sp <- sp[order(sp$spectrum_index), , drop = FALSE]
        col <- function(c) if (c %in% names(sp)) sp[[c]] else
            rep(NA, nrow(sp))
        pol <- c("-1" = "NEGATIVE", "1" = "POSITIVE")[as.character(
            col("scan_polarity"))]
        for (i in seq_len(nrow(sp))) {
            kv <- if (dialect == "massbank") c(
                "ACCESSION" = col("id")[i],
                "CH$NAME" = col("x_mspurity_name")[i],
                "CH$FORMULA" = col("x_mspurity_formula")[i],
                "CH$LINK: INCHIKEY" = col("x_mspurity_inchikey")[i],
                "AC$INSTRUMENT_TYPE" = col("x_mspurity_instrument_type")[i],
                "AC$MASS_SPECTROMETRY: MS_TYPE" = paste0("MS", col("ms_level")[i]),
                "AC$MASS_SPECTROMETRY: ION_MODE" = unname(pol[[i]]),
                "AC$MASS_SPECTROMETRY: COLLISION_ENERGY" =
                    col("x_mspurity_collision_energy_text")[i],
                "MS$FOCUSED_ION: PRECURSOR_M/Z" = num(col("selected_ion_mz")[i]),
                "MS$FOCUSED_ION: PRECURSOR_TYPE" =
                    col("x_mspurity_precursor_type")[i],
                "PK$NUM_PEAK" = length(sp$mz[[i]]))
            else c(
                "Name" = col("x_mspurity_name")[i], "DB#" = col("id")[i],
                "InChIKey" = col("x_mspurity_inchikey")[i],
                "Formula" = col("x_mspurity_formula")[i],
                "PrecursorMZ" = num(col("selected_ion_mz")[i]),
                "Precursor_type" = col("x_mspurity_precursor_type")[i],
                "Spectrum_type" = paste0("MS", col("ms_level")[i]),
                "Instrument_type" = col("x_mspurity_instrument_type")[i],
                "Ion_mode" = unname(c(NEGATIVE = "N", POSITIVE = "P")[pol[[i]]]),
                "Collision_energy" = col("x_mspurity_collision_energy_text")[i],
                "Num Peaks" = length(sp$mz[[i]]))
            kv <- kv[!is.na(kv)]
            sep <- if (dialect == "massbank") ": " else ": "
            out <- c(out, paste0(names(kv), sep, kv))
            if (dialect == "massbank")
                out <- c(out, "PK$PEAK: m/z int. rel.int.")
            out <- c(out, paste0(if (dialect == "massbank") "  " else "",
                                 num(sp$mz[[i]]), " ",
                                 num(sp$intensity[[i]])),
                     if (dialect == "massbank") "//" else "")
        }
    }
    writeLines(out, file, useBytes = TRUE)
    invisible(file)
}
