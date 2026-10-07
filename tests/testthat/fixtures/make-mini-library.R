# Build the mini library fixtures used by the mzStack library tests.
#
# Eight spectra of the msp2db library spectralMatching() downloads (Zenodo
# record 18700802), written three ways that must convert to the same
# spectra:
#
#   inst/extdata/tests/mzstack/library/mini_library.sqlite   msp2db schema
#   inst/extdata/tests/mzstack/library/mini_massbank.msp     MassBank keys
#   inst/extdata/tests/mzstack/library/mini_mona.msp         MoNA/NIST keys
#
# The MSP files are written here, independently of msPurity's converters, with
# every number in a form that reads back to the database's own double.
#
# Run from the package root, with the library in BiocFileCache:
#
#   Rscript tests/testthat/fixtures/make-mini-library.R

accessions <- c("PR100407", "ML005101", "CCMSLIB00003740024",
                "CCMSLIB00000479720", "KNA00052", "CE000616",
                "CCMSLIB00000577898", "AU101851")

bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
hit <- BiocFileCache::bfcquery(bfc, "msPurity_library_spectra_db",
                               "rname", exact = TRUE)
stopifnot(nrow(hit) == 1L)
src_db <- BiocFileCache::bfcrpath(bfc, rids = hit$rid)

out <- file.path("inst", "extdata", "tests", "mzstack", "library")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
db <- file.path(out, "mini_library.sqlite")
unlink(db)

src <- DBI::dbConnect(RSQLite::SQLite(), src_db, flags = RSQLite::SQLITE_RO)
dst <- DBI::dbConnect(RSQLite::SQLite(), db)
q <- function(sql) suppressWarnings(DBI::dbGetQuery(src, sql))
for (sql in q("SELECT sql FROM sqlite_master WHERE type = 'table'")$sql)
    DBI::dbExecute(dst, sql)
inlist <- paste0("'", accessions, "'", collapse = ",")
meta <- q(sprintf("SELECT * FROM library_spectra_meta WHERE accession IN (%s)
                   ORDER BY id", inlist))
stopifnot(nrow(meta) == length(accessions))
peaks <- q(sprintf("SELECT * FROM library_spectra WHERE
                    library_spectra_meta_id IN (%s) ORDER BY id",
                   paste(meta$id, collapse = ",")))
sources <- q(sprintf("SELECT * FROM library_spectra_source WHERE id IN (%s)",
                     paste(unique(meta$library_spectra_source_id),
                           collapse = ",")))
cmp <- q(sprintf("SELECT * FROM metab_compound WHERE inchikey_id IN (%s)",
                 paste0("'", unique(meta$inchikey_id), "'", collapse = ",")))
DBI::dbWriteTable(dst, "library_spectra_source", sources, append = TRUE)
DBI::dbWriteTable(dst, "library_spectra_meta", meta, append = TRUE)
DBI::dbWriteTable(dst, "library_spectra", peaks, append = TRUE)
DBI::dbWriteTable(dst, "metab_compound", cmp, append = TRUE)
DBI::dbDisconnect(dst)
DBI::dbDisconnect(src)

num <- function(x) ifelse(is.na(x), NA, sprintf("%.17g", x))
k <- match(meta$inchikey_id, cmp$inchikey_id)
pol <- toupper(meta$polarity)

mb <- character()
mona <- character()
for (i in seq_len(nrow(meta))) {
    p <- peaks[peaks$library_spectra_meta_id == meta$id[i], ]
    kv <- c("ACCESSION" = meta$accession[i],
            "RECORD_TITLE" = meta$name[i],
            "CH$NAME" = cmp$name[k[i]],
            "CH$FORMULA" = cmp$molecular_formula[k[i]],
            "CH$LINK: INCHIKEY" = meta$inchikey_id[i],
            "AC$INSTRUMENT_TYPE" = meta$instrument_type[i],
            "AC$MASS_SPECTROMETRY: MS_TYPE" = "MS2",
            "AC$MASS_SPECTROMETRY: ION_MODE" = pol[i],
            "AC$MASS_SPECTROMETRY: COLLISION_ENERGY" = meta$collision_energy[i],
            "MS$FOCUSED_ION: PRECURSOR_M/Z" = num(meta$precursor_mz[i]),
            "MS$FOCUSED_ION: PRECURSOR_TYPE" = meta$precursor_type[i],
            "PK$NUM_PEAK" = nrow(p))
    kv <- kv[!is.na(kv) & nzchar(kv)]
    ## MassBank writes a subtag and its value separated by a space.
    sep <- ifelse(grepl(": ", names(kv)), " ", ": ")
    mb <- c(mb, paste0(names(kv), sep, kv), "PK$PEAK: m/z int. rel.int.",
            paste0("  ", num(p$mz), " ", num(p$i), " ",
                   round(p$i / max(p$i) * 999)), "//")
    kv <- c("Name" = cmp$name[k[i]], "DB#" = meta$accession[i],
            "InChIKey" = meta$inchikey_id[i],
            "Formula" = cmp$molecular_formula[k[i]],
            "Instrument_type" = meta$instrument_type[i],
            "Spectrum_type" = "MS2",
            "Ion_mode" = substr(pol[i], 1, 1),
            "Collision_energy" = meta$collision_energy[i],
            "PrecursorMZ" = num(meta$precursor_mz[i]),
            "Precursor_type" = meta$precursor_type[i],
            "Num Peaks" = nrow(p))
    kv <- kv[!is.na(kv) & nzchar(kv)]
    mona <- c(mona, paste0(names(kv), ": ", kv),
              paste0(num(p$mz), " ", num(p$i)), "")
}
writeLines(mb, file.path(out, "mini_massbank.msp"))
writeLines(mona, file.path(out, "mini_mona.msp"))
message("Wrote ", nrow(meta), " spectra to ", out)
