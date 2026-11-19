context("purityA with the legacy slots switched off")

suppressPackageStartupMessages(library(xcms))

test_that("results are unchanged when only the Spectra slots are filled", {
    ref <- .parity_ref("purityA")
    old <- options(msPurity.legacySlots = FALSE)
    on.exit(options(old))
    res <- .parity_purityA(via = "accessors")
    for (stage in names(ref))
        expect_equal(res[[stage]], ref[[stage]],
                     tolerance = 1e-12, label = stage)
})

# Use allfrag = TRUE so the export rebuilds all_frag_scans when its slot is empty.
.workflow_inputs <- function() {
    q <- .parity_quiet
    pa <- q(purityA(unname(.lcmsms_paths())))
    xcmsObj <- .golden_xcms("msms_only_xcmsnexp.rds")
    pa <- q(frag4feature(pa, xcmsObj))
    pa <- q(filterFragSpectra(pa, plim = 0.7, snr = 3, allfrag = TRUE))
    pa <- q(averageAllFragSpectra(q(averageInterFragSpectra(
        q(averageIntraFragSpectra(pa))))))
    list(pa = pa, xcmsObj = xcmsObj)
}

# Each format gets a fresh destination and the same analysis for its slot setting.
.workflow_db <- local({
    inputs <- list()
    function(format, legacy) {
        old <- options(msPurity.legacySlots = legacy)
        on.exit(options(old))
        mode <- if (legacy) "on" else "off"
        if (is.null(inputs[[mode]]))
            inputs[[mode]] <<- .workflow_inputs()
        input <- inputs[[mode]]
        out <- tempfile("legacy-off-")
        dir.create(out)
        path <- .parity_quiet(createDatabase(input$pa, input$xcmsObj, outDir = out,
                                             dbName = paste0("out.", format),
                                             format = format))
        list(pa = input$pa, path = path)
    }
})

test_that("the SQLite output does not depend on the legacy slots", {
    on <- .workflow_db("sqlite", legacy = TRUE)
    off <- .workflow_db("sqlite", legacy = FALSE)
    expect_identical(length(off$pa@grped_ms2), 0L)
    expect_identical(nrow(off$pa@all_frag_scans), 0L)
    expect_identical(.sqlite_dump(off$path), .sqlite_dump(on$path))
})

test_that("the Parquet output does not depend on the legacy slots", {
    skip_if_no_arrow()
    on <- .workflow_db("parquet", legacy = TRUE)
    off <- .workflow_db("parquet", legacy = FALSE)
    tables <- names(.manifest(on$path)$results$tables)
    expect_setequal(names(.manifest(off$path)$results$tables), tables)
    # Tables that record when and where the dataset was written differ, as
    # do the times and ids of the processing history.
    read <- function(path, t) {
        d <- .read_table(path, t)
        if (!is.null(d$x_mspurity_history))
            d$x_mspurity_history <- lapply(d$x_mspurity_history, function(h) {
                if (is.na(h)) return(NULL)
                lapply(jsonlite::fromJSON(h, simplifyVector = FALSE),
                       function(r) {
                           r$time <- NULL
                           r$step_id <- NULL
                           r$input$previous_steps <- NULL
                           r
                       })
            })
        d
    }
    for (t in setdiff(tables, c("activity", "activity_input", "software",
                                "conversion")))
        expect_identical(read(off$path, t), read(on$path, t), label = t)
    v <- validateParquet(off$path)
    expect_identical(sum(v$level == "MUST"), 0L)
})
