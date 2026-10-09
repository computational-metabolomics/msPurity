context("purityA with the legacy slots switched off")

suppressPackageStartupMessages(library(xcms))

test_that("results are unchanged when only the Spectra slots are filled", {
    ref <- .parity_ref("purityA")
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    res <- .parity_purityA(via = "accessors")
    for (stage in names(ref)) {
        r <- ref[[stage]]
        if (is.list(r) && !is.null(r$all_frag_scans))
            r <- .parity_fix_allfrag(r)
        expect_identical(res[[stage]], r, label = stage)
    }
})

# The full workflow up to createDatabase(), with allfrag = TRUE so that the
# frozen all_frag_scans table is rebuilt when the legacy slot is empty.
.workflow_db <- function(format) {
    q <- .parity_quiet
    pa <- q(purityA(unname(.lcmsms_paths())))
    xcmsObj <- .golden_xcms("msms_only_xcmsnexp.rds")
    pa <- q(frag4feature(pa, xcmsObj))
    pa <- q(filterFragSpectra(pa, plim = 0.7, snr = 3, allfrag = TRUE))
    pa <- q(averageAllFragSpectra(q(averageInterFragSpectra(
        q(averageIntraFragSpectra(pa))))))
    out <- tempfile("legacy-off-")
    dir.create(out)
    list(pa = pa, path = q(createDatabase(pa, xcmsObj, outDir = out,
                                          dbName = paste0("out.", format),
                                          format = format)))
}

test_that("the SQLite output does not depend on the legacy slots", {
    on <- .workflow_db("sqlite")
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    off <- .workflow_db("sqlite")
    expect_identical(length(off$pa@grped_ms2), 0L)
    expect_identical(nrow(off$pa@all_frag_scans), 0L)
    expect_identical(.sqlite_dump(off$path), .sqlite_dump(on$path))
})

test_that("the mzStack output does not depend on the legacy slots", {
    skip_if_no_mzstack()
    on <- .workflow_db("mzstack")
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    off <- .workflow_db("mzstack")
    tables <- names(.manifest(on$path)$results$tables)
    expect_setequal(names(.manifest(off$path)$results$tables), tables)
    # Tables that record when and where the dataset was written differ.
    for (t in setdiff(tables, c("activity", "activity_input", "software",
                                "conversion")))
        expect_identical(.read_table(off$path, t), .read_table(on$path, t),
                         label = t)
    v <- validateMzstack(off$path)
    expect_identical(sum(v$level == "MUST"), 0L)
})
