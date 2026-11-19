# Run from the msPurity checkout. Prepare input once, then export the same
# fixture in separate R processes using the candidate and baseline checkouts.
# Rscript inst/benchmarks/parquet-export.R prepare DIR N
# Rscript inst/benchmarks/parquet-export.R bridge DIR N BASELINE_CHECKOUT
# Rscript inst/benchmarks/parquet-export.R direct DIR N CANDIDATE_CHECKOUT
# Rscript inst/benchmarks/parquet-export.R compare DIR N
args <- commandArgs(TRUE)
action <- match.arg(args[1L], c("prepare", "bridge", "direct", "compare"))
n <- as.integer(args[3L])
stopifnot(n >= 2L)
dir <- file.path(args[2L], paste0("runs-", n))
dir.create(dir, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(if (length(args) >= 4L) args[4L] else ".", helpers = FALSE)
suppressPackageStartupMessages(library(xcms))
fixture <- file.path(dir, "input.rds")
if (action == "prepare") {
    src <- list.files(system.file("extdata/lcms/mzML", package = "msPurityData"),
                      pattern = "MSMS.*mzML$", full.names = TRUE)
    stopifnot(length(src) == 2L)
    files <- file.path(normalizePath(dir), sprintf("sample-%03d.mzML", seq_len(n)))
    stopifnot(all(file.copy(rep(src, length.out = n), files, overwrite = TRUE)))
    serial <- BiocParallel::SerialParam()
    xe <- MsExperiment::readMsExperiment(files, BPPARAM = serial)
    xe <- xcms::findChromPeaks(xe, param = xcms::CentWaveParam(
        ppm = 15, peakwidth = c(3, 30), prefilter = c(3, 100)),
        msLevel = 1L, BPPARAM = serial)
    xe <- xcms::groupChromPeaks(xe, param = xcms::PeakDensityParam(
        sampleGroups = rep(1L, n), minFraction = 0.5))
    pa <- purityA(files, cores = 1L, plotP = FALSE)
    pa <- frag4feature(pa, xe, convert2RawRT = FALSE)
    pa <- filterFragSpectra(pa, allfrag = TRUE)
    pa <- averageIntraFragSpectra(pa, cores = 1L)
    pa <- averageInterFragSpectra(pa, cores = 1L)
    pa <- averageAllFragSpectra(pa, cores = 1L)
    saveRDS(list(pa = pa, xe = xe), fixture, compress = FALSE)
} else if (action == "compare") {
    source("tests/testthat/helper-sqlite-golden.R")
    source("tests/testthat/helper-parquet.R")
    source("tests/testthat/helper-parquet-golden.R")
    a <- .parquet_canon(file.path(dir, "direct"))
    b <- .parquet_canon(file.path(dir, "bridge"))
    # This fixture retains failed output peaks (rmp = FALSE). Its only
    # logical addition is the ledger of attempts which produced NULL.
    ll <- a$tables$loss_ledger
    attempt <- grepl("^av_(intra|inter|all): feature ", ll$source_construct)
    # All three methods were run. Independently enumerate their NULL
    # results in the original, nested feature/sample representation.
    legacy <- readRDS(fixture)$pa@av_spectra
    expected <- unlist(lapply(names(legacy), function(group) {
        av <- legacy[[group]]
        empty <- vapply(av[c("av_all", "av_inter")], is.null, logical(1))
        refs <- paste0(names(empty)[empty], ": feature ", group, recycle0 = TRUE)
        intra <- av$av_intra
        if (!length(intra)) return(c(refs, paste0("av_intra: feature ", group)))
        samples <- names(intra)[vapply(intra, is.null, logical(1))]
        c(refs, paste0("av_intra: feature ", group, ", sample ", samples, recycle0 = TRUE))
    }), use.names = FALSE)
    stopifnot(identical(sort(ll$source_construct[attempt]), sort(as.character(expected))),
              all(ll$disposition[attempt] == "not_present_in_source"),
              all(ll$affected_rows[attempt] == 1L))
    ll <- ll[!attempt, , drop = FALSE]
    ll$loss_ledger_id_ <- seq_len(nrow(ll))
    a$tables$loss_ledger <- ll
    a$manifest$results$tables$loss_ledger$rows <- nrow(ll)
    for (name in names(b$tables)) {
        by <- unlist(b$manifest$results$tables[[name]]$sorted_by)
        sort_rows <- function(x) {
            at <- do.call(order, c(unname(as.list(x[by])), list(method = "radix")))
            x <- x[at, , drop = FALSE]
            rownames(x) <- NULL
            x
        }
        expected <- sort_rows(b$tables[[name]])
        actual <- sort_rows(a$tables[[name]])
        if (!identical(expected, actual))
            stop(name, ": ", paste(all.equal(expected, actual), collapse = "; "))
        a$manifest$results$tables[[name]]$parts <- NULL
        b$manifest$results$tables[[name]]$parts <- NULL
    }
    stopifnot(identical(names(a$tables), names(b$tables)),
              identical(a$runs, b$runs), identical(a$manifest, b$manifest),
              identical(a$json, b$json))
    print(vapply(a$tables, NROW, integer(1)))
} else {
    input <- readRDS(fixture)
    if (action == "bridge") {
        stopifnot(length(args) >= 4L)
        options(msPurity.parquet.route = "bridge")
    }
    gc()
    print(system.time(createDatabase(input$pa, input$xe, outDir = dir,
        dbName = action, format = "parquet", overwrite = TRUE)))
    print(sessionInfo())
}
