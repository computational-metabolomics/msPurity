# Parity references for the move to Spectra containers.
#
# Each workflow below runs msPurity from the raw msPurityData files and
# returns the results in their legacy shapes (puritydf, grped_df, grped_ms2,
# all_frag_scans, av_spectra, MSP text and the DIMS peak lists). The
# references in msPurityData were written by the code before the move,
# and later results must agree within floating-point rounding. Explicit
# identity and saved-object checks remain exact.

# Feature-link and provenance columns appended after the legacy columns.
.f4f_added_cols <- c("precurMtchRTadj", "rtMapStatus", "precurMzSource",
                     "is_filled", "rtminWin", "rtmaxWin", "rtApexDist",
                     "matchRank", "linkId")

# grped_df without the appended columns, in the legacy column order.
.parity_grped_df <- function(g) {
    added <- intersect(names(g), .f4f_added_cols)
    if (!length(added))
        return(g)
    g[, setdiff(names(g), added), drop = FALSE]
}

# The legacy view of a purityA object, read from the legacy slots or, once
# they exist, rebuilt from the Spectra slots by the accessors. The same
# references check both.
.parity_view <- function(pa, via = c("slots", "accessors")) {
    via <- match.arg(via)
    acc <- function(f, slot) {
        if (via == "accessors")
            get(f, envir = asNamespace("msPurity"))(pa, legacy = TRUE)
        else methods::slot(pa, slot)
    }
    list(puritydf = acc("purityTable", "puritydf"),
         grped_df = .parity_grped_df(pa@grped_df),
         grped_ms2 = acc("groupedSpectra", "grped_ms2"),
         all_frag_scans = acc("allFragSpectra", "all_frag_scans"),
         av_spectra = acc("averagedSpectra", "av_spectra"))
}

# plyr leaves split attributes on some lists, filtering leaves gaps in data
# frame row names, and cbind() names some matrix rows; none of these is part
# of the result, so they are dropped before comparing.
.parity_strip <- function(x) {
    if (is.data.frame(x) && !nrow(x) && !ncol(x))
        return(data.frame())
    if (is.data.frame(x)) {
        attr(x, "split_type") <- NULL
        attr(x, "split_labels") <- NULL
        rownames(x) <- NULL
        return(x)
    }
    if (is.matrix(x)) {
        rownames(x) <- NULL
        return(x)
    }
    if (is.list(x)) {
        a <- attributes(x)
        x <- lapply(x, .parity_strip)
        keep <- setdiff(names(a), c("split_type", "split_labels"))
        attributes(x)[keep] <- a[keep]
        return(x)
    }
    x
}

.msp_text <- function(pth) {
    s <- readChar(pth, file.info(pth)$size)
    # The legacy writer adds CRLF on Windows through a text connection,
    # which can add a second CR. Preserve line boundaries and all content.
    s <- gsub("\r{1,2}\n", "\n", s)
    gsub("msPurity version:\\d+\\.\\d+\\.\\d+", "", s)
}

.parity_msp <- function(pa) {
    metadata <- data.frame("grpid" = c(162, 187),
                           "MS$FOCUSED_ION: PRECURSOR_TYPE" = c("[M+H]+", "[M+H]+"),
                           "AC$MASS_SPECTROMETRY: ION_MODE" = c("POSITIVE", "POSITIVE"),
                           "CH$NAME:" = c("Unknown", "Methionine"),
                           check.names = FALSE, stringsAsFactors = FALSE)
    methods <- c("all", "max", "av_inter", "av_intra", "av_all")
    run <- function(m, ...) {
        f <- tempfile(fileext = ".msp")
        createMSP(pa, msp_file_pth = f, method = m, ...)
        .msp_text(f)
    }
    # method = "max" over all groups failed in the reference code, so it is
    # only covered for grpid 162 and 187; test.2_purityA_fixes.R covers the
    # rest.
    every <- setdiff(methods, "max")
    out <- c(lapply(methods, run, metadata = metadata, xcms_groupids = c(162, 187)),
             lapply(every, run))
    names(out) <- c(paste0(methods, "_metadata"), paste0(every, "_all_groups"))
    out
}

# Add a link from one scan to a second feature. The scan of the first row of
# grpid 187 is also linked to grpid 162, so the scan belongs to two features.
.parity_two_features <- function(pa) {
    g <- pa@grped_df
    src <- which(as.character(g$grpid) == "187")[1]
    row <- g[src, ]
    row$grpid <- factor("162", levels = levels(g$grpid))
    dst <- max(which(as.character(g$grpid) == "162"))
    g <- rbind(g[seq_len(dst), ], row, g[-seq_len(dst), ])
    rownames(g) <- NULL
    ms2 <- pa@grped_ms2
    k <- match(src, which(as.character(pa@grped_df$grpid) == "187"))
    ms2[["162"]] <- c(ms2[["162"]], list(ms2[["187"]][[k]]))
    if ("linkId" %in% names(g)) {
        # This hand-built association has its own occurrences and source.
        input <- msPurity:::.pa_processing_input(pa, "frag4feature",
            source = pa@fragSpectra)
        pa@grped_df <- g
        pa <- msPurity:::.pa_graph_mint_links(pa)
        pa <- msPurity:::.pa_processing_record(pa, "frag4feature",
            list(association_input = "test_fixture", added_pid = row$pid,
                added_feature = "162"), Sys.time(), input)
    } else {
        pa@grped_df <- g
    }
    pa@grped_ms2 <- ms2
    pa
}

.parity_quiet <- function(expr) {
    utils::capture.output(res <- suppressWarnings(suppressMessages(expr)))
    res
}

# Keep the stage objects so both views can be checked against the same run.
.parity_stages <- function() {
    q <- .parity_quiet
    p <- unname(.lcmsms_paths())
    out <- list()

    # These references used a fixed window around precursorMz.
    # Recorded event centres are covered by the window tests.
    pa1 <- q(purityA(p, offsets = c(0.5, 0.5)))
    out$purityA <- pa1

    pa2 <- q(frag4feature(pa1, .golden_xcms("msms_only_xcmsnexp.rds")))
    out$frag4feature_xcmsnexp <- pa2
    out$frag4feature_xset <- q(frag4feature(pa1, .golden_xcms("msms_only_xset.rds")))
    out$frag4feature_group <- q(frag4feature(
        pa1, .golden_xcms("msms_only_xcmsnexp.rds"), useGroup = TRUE))

    pa3 <- q(filterFragSpectra(pa2, plim = 0.7, snr = 3, allfrag = TRUE))
    out$filter_allfrag <- pa3
    pa3r <- q(filterFragSpectra(pa2, plim = 0.7, snr = 3, rmp = TRUE))
    out$filter_rmp <- pa3r

    pa4 <- q(averageIntraFragSpectra(pa3))
    pa5 <- q(averageInterFragSpectra(pa4))
    pa6 <- q(averageAllFragSpectra(pa5))
    out$average_filtered <- pa6
    # Averaging after filterFragSpectra(rmp = TRUE) failed in the reference
    # code when every spectrum of a feature lost all its peaks, so it has no
    # reference; test.2_purityA_fixes.R covers it.
    out$average_unfiltered <- q(averageAllFragSpectra(pa2))
    out$msp <- q(.parity_msp(pa6))

    pa2d <- .parity_two_features(pa2)
    pa6d <- q(averageAllFragSpectra(q(averageInterFragSpectra(
        q(averageIntraFragSpectra(
            q(filterFragSpectra(pa2d, plim = 0.7, snr = 3))))))))
    out$two_features <- pa6d
    out$two_features_msp <- q(.parity_msp(pa6d))
    out
}

.parity_purityA <- function(via = "slots", stages = .parity_stages()) {
    lapply(stages, function(stage) {
        if (methods::is(stage, "purityA")) .parity_view(stage, via) else stage
    })
}

.dims_paths <- function() {
    system.file("extdata", "dims", "mzML", package = "msPurityData")
}

# The retired DIMS workflow, returning the peak lists by stage.
.parity_purityD <- function() .with_dims_deprecation({
    q <- .parity_quiet
    inDF <- Getfiles(.dims_paths(), pattern = ".mzML", check = FALSE)
    pd <- q(purityD(inDF, mzML = TRUE))
    pd <- q(averageSpectra(pd, snMeth = "median", snthr = 5))
    out <- list(averageSpectra = pd@avPeaks)
    pd <- q(filterp(pd, thr = 5000, rsd = 10))
    out$filterp <- pd@avPeaks
    pd <- q(subtract(pd))
    out$subtract <- pd@avPeaks
    pd <- q(dimsPredictPurity(pd))
    out$dimsPredictPurity <- pd@avPeaks
    pd <- q(groupPeaks(pd))
    out$groupPeaks <- pd@groupedPeaks
    out
})

# Matrices saved before R 4.2 have no column names in grped_ms2.
.label_ms2 <- function(ms2) {
    lapply(ms2, lapply, function(m) {
        if (!is.null(m) && ncol(m) == 2) colnames(m) <- c("mz", "intensity")
        m
    })
}

# Compare a purityA object with one saved before the Spectra slots existed.
# Every slot of the saved object is compared, which is what comparing the
# whole objects did before. The accessors of the current object must also
# rebuild its own legacy slots.
#
# The saved filtered objects were made from grped_ms2 matrices without column
# names (xcms before R 4.2). setFlagMatrix() names the row of a one-peak
# matrix only when its input has column names, so for grped_ms2 the row names
# of matrices are not compared; every value is.
expect_legacy_equal <- function(object, expected) {
    unrow <- function(ms2) lapply(ms2, lapply, function(m) {
        if (is.matrix(m)) rownames(m) <- NULL
        m
    })
    for (s in setdiff(names(attributes(expected)), "class")) {
        a <- attr(object, s, exact = TRUE)
        b <- attr(expected, s, exact = TRUE)
        if (s == "filter_frag_params" && length(b)) {
            # Historical objects have no level-selection settings.
            for (name in intersect(c("msLevel", "effective_msLevel"),
                    names(a))) {
                expect_identical(a[[name]], 2L, label = name)
                a[[name]] <- NULL
            }
        }
        if (s %in% c("av_intra_params", "av_inter_params", "av_all_params") &&
            "blockSize" %in% names(a) && !"blockSize" %in% names(b)) {
            # Historical averages predate the block-size execution control.
            expect_identical(a$blockSize, 5000L, label = paste(s, "blockSize"))
            a$blockSize <- NULL
        }
        if (s == "grped_ms2") {
            a <- unrow(a)
            b <- unrow(b)
        }
        if (s == "grped_df" && "linkId" %in% names(a) && !"linkId" %in% names(b)) {
            expect_type(a$linkId, "integer")
            expect_false(anyNA(a$linkId) || anyDuplicated(a$linkId) > 0L)
            a$linkId <- NULL
        }
        expect_equal(a, b, label = paste("slot", s))
    }
    if (!isTRUE(getOption("msPurity.legacySlots", TRUE)))
        return(invisible(object))
    expect_identical(purityTable(object), object@puritydf)
    expect_equal(.parity_strip(groupedSpectra(object)),
                 .parity_strip(.label_ms2(object@grped_ms2)))
    afs <- allFragSpectra(object)
    if (nrow(object@all_frag_scans))
        expect_equal(.parity_strip(afs), .parity_strip(object@all_frag_scans))
    else expect_identical(nrow(afs), 0L)
    expect_equal(.parity_strip(averagedSpectra(object)),
                 .parity_strip(object@av_spectra))
    invisible(object)
}
