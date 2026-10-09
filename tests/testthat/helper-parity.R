# Parity references for the move to Spectra containers.
#
# Each workflow below runs msPurity from the raw msPurityData files and
# returns the results in their legacy shapes (puritydf, grped_df, grped_ms2,
# all_frag_scans, av_spectra, MSP text and the DIMS peak lists). The
# references in fixtures/parity were written by the code before the move,
# and every later change must reproduce them exactly.

.parity_dir <- function() {
    d <- testthat::test_path("fixtures", "parity")
    if (dir.exists(d)) d else file.path("tests", "testthat", "fixtures", "parity")
}

.parity_ref <- function(name) readRDS(file.path(.parity_dir(), paste0(name, ".rds")))

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
         grped_df = pa@grped_df,
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

# The legacy all_frag_scans names scans by a count that skips scans without
# peaks, so after such a scan its pid is wrong (the scan column is right).
# allFragSpectra() returns the table with each row's pid taken from its file
# and scan, and the purity and overall flags of that scan.
.parity_fix_allfrag <- function(view) {
    afs <- view$all_frag_scans
    if (!nrow(afs))
        return(view)
    p <- view$puritydf
    afs$pid <- p$pid[match(paste(afs$fileid, as.character(afs$scan)),
                           paste(as.character(p$fileid), p$seqNum))]
    afs$purity_pass_flag <- p$inPurity[match(afs$pid, p$pid)] > 0.7
    afs$pass_flag <- afs$purity_pass_flag & afs$intensity_pass_flag &
        afs$ra_pass_flag & afs$snr_pass_flag
    view$all_frag_scans <- afs
    view
}

.msp_text <- function(pth) {
    s <- readChar(pth, file.info(pth)$size)
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
    pa@grped_df <- g
    pa@grped_ms2 <- ms2
    pa
}

.parity_quiet <- function(expr) {
    utils::capture.output(res <- suppressWarnings(suppressMessages(expr)))
    res
}

# Run every purityA workflow and return the legacy views by stage name.
.parity_purityA <- function(via = "slots") {
    q <- .parity_quiet
    p <- unname(.lcmsms_paths())
    out <- list()

    pa1 <- q(purityA(p))
    out$purityA <- .parity_view(pa1, via)

    pa2 <- q(frag4feature(pa1, .golden_xcms("msms_only_xcmsnexp.rds")))
    out$frag4feature_xcmsnexp <- .parity_view(pa2, via)
    out$frag4feature_xset <- .parity_view(
        q(frag4feature(pa1, .golden_xcms("msms_only_xset.rds"))), via)
    out$frag4feature_group <- .parity_view(
        q(frag4feature(pa1, .golden_xcms("msms_only_xcmsnexp.rds"),
                       useGroup = TRUE)), via)

    pa3 <- q(filterFragSpectra(pa2, plim = 0.7, snr = 3, allfrag = TRUE))
    out$filter_allfrag <- .parity_view(pa3, via)
    pa3r <- q(filterFragSpectra(pa2, plim = 0.7, snr = 3, rmp = TRUE))
    out$filter_rmp <- .parity_view(pa3r, via)

    pa4 <- q(averageIntraFragSpectra(pa3))
    pa5 <- q(averageInterFragSpectra(pa4))
    pa6 <- q(averageAllFragSpectra(pa5))
    out$average_filtered <- .parity_view(pa6, via)
    # Averaging after filterFragSpectra(rmp = TRUE) failed in the reference
    # code when every spectrum of a feature lost all its peaks, so it has no
    # reference; test.2_purityA_fixes.R covers it.
    out$average_unfiltered <- .parity_view(q(averageAllFragSpectra(pa2)), via)
    out$msp <- q(.parity_msp(pa6))

    pa2d <- .parity_two_features(pa2)
    pa6d <- q(averageAllFragSpectra(q(averageInterFragSpectra(
        q(averageIntraFragSpectra(
            q(filterFragSpectra(pa2d, plim = 0.7, snr = 3))))))))
    out$two_features <- .parity_view(pa6d, via)
    out$two_features_msp <- q(.parity_msp(pa6d))
    out
}

.dims_paths <- function() {
    system.file("extdata", "dims", "mzML", package = "msPurityData")
}

# The DIMS workflow from the vignette, returning the peak lists by stage.
.parity_purityD <- function() {
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
}

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
        if (s == "grped_ms2") {
            a <- unrow(a)
            b <- unrow(b)
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
