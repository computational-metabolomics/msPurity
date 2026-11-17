context("frag4feature: retention time mapping and linking")

.f4f <- function(name) get(name, envir = asNamespace("msPurity"))
.q <- function(expr) suppressWarnings(suppressMessages(expr))
.links <- function(pa) {
    g <- pa@grped_df
    sort(paste(as.character(g$grpid), g$pid,
               if ("cid" %in% names(g)) g$cid else ""))
}

# Shift the adjusted times of the given samples of an xcmsSet, with their
# peaks and features, as if retention time correction had moved them.
.shift_xset <- function(x, samples, by) {
    for (s in samples) {
        x@rt$corrected[[s]] <- x@rt$corrected[[s]] + by
        r <- x@peaks[, "sample"] == s
        x@peaks[r, c("rt", "rtmin", "rtmax")] <-
            x@peaks[r, c("rt", "rtmin", "rtmax")] + by
    }
    if (length(samples) == length(x@rt$corrected))
        x@groups[, c("rtmed", "rtmin", "rtmax")] <-
            x@groups[, c("rtmed", "rtmin", "rtmax")] + by
    x
}

# Shift the MS/MS times of one file of a purityA object.
.shift_pa <- function(pa, fileid, by) {
    pa <- .f4f(".pa_update")(pa)
    sp <- pa@spectra
    f <- as.character(sp$fileid) == as.character(fileid)
    for (v in c("precursorRT", "retentionTime")) {
        x <- sp[[v]]
        x[f] <- x[f] + by
        sp[[v]] <- x
    }
    pa@spectra <- sp
    pa
}

test_that("times are interpolated between knots and shifted beyond them", {
    interp <- .f4f(".f4f_interp")
    raw <- c(10, 20, 30, 40)
    adj <- c(11, 23, 31, 39)
    expect_identical(interp(raw, raw, adj), adj)
    expect_equal(interp(c(15, 35), raw, adj), c(17, 35))
    # beyond the knots, the shift at the nearest end
    expect_equal(interp(c(5, 50), raw, adj), c(6, 49))
    expect_equal(interp(interp(c(12, 33), raw, adj), adj, raw), c(12, 33))
    expect_identical(interp(c(NA, 20), raw, adj), c(NA, 23))
})

test_that("files are matched by name, whatever the path separator", {
    idx <- .f4f(".f4f_file_index")
    expect_identical(.f4f(".f4f_basename")(c("C:\\data\\a.mzML", "/d/b.mzML",
                                              "c.mzML")),
                     c("a.mzML", "b.mzML", "c.mzML"))
    expect_identical(idx(c("/x/b.mzML", "/x/a.mzML", "/x/z.mzML"),
                         c("C:\\y\\a.mzML", "/y/b.mzML")),
                     c(2L, 1L, NA))
    expect_error(idx("a.mzML", c("/1/a.mzML", "/2/a.mzML")), "several files")
})

test_that("individual links are unchanged when a file's adjusted times shift", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xset.rds")
    ref <- .q(frag4feature(pa, x))
    shifted <- .shift_xset(x, 2, 20)
    expect_identical(.links(.q(frag4feature(pa, shifted))), .links(ref))
    # Without mapping, the shifted file loses its links
    raw <- .q(frag4feature(pa, shifted, convert2RawRT = FALSE))
    expect_lt(sum(raw@grped_df$fileid == 2), sum(ref@grped_df$fileid == 2))
})

test_that("grouped links are unchanged when every adjusted time shifts", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xset.rds")
    ref <- .q(frag4feature(pa, x, useGroup = TRUE))
    shifted <- .shift_xset(x, 1:2, 20)
    got <- .q(frag4feature(pa, shifted, useGroup = TRUE))
    expect_identical(.links(got), .links(ref))
    expect_equal(got@grped_df$precurMtchRTadj, ref@grped_df$precurMtchRTadj + 20)
    raw <- .q(frag4feature(pa, shifted, useGroup = TRUE, convert2RawRT = FALSE))
    expect_lt(nrow(raw@grped_df), nrow(ref@grped_df))
})

test_that("MS/MS times on XCMS scan times map to XCMS's adjusted times exactly", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    raw <- rtime(x, adjusted = FALSE, bySample = TRUE)
    adj <- rtime(x, adjusted = TRUE, bySample = TRUE)
    for (mode in c(FALSE, TRUE)) {
        g <- .q(frag4feature(pa, x, useGroup = mode))@grped_df
        expect_true(all(g$rtMapStatus == "aligned"))
        for (s in 1:2) {
            r <- g$fileid == s
            k <- match(g$precurMtchRT[r], raw[[s]])
            expect_true(any(!is.na(k)))
            expect_identical(g$precurMtchRTadj[r][!is.na(k)],
                             unname(adj[[s]][k[!is.na(k)]]))
        }
    }
})

test_that("filled peaks are left out unless includeFilled = TRUE", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    xf <- .q(fillChromPeaks(x, param = FillChromPeaksParam(expandMz = 0,
             expandRt = 0, ppm = 0), BPPARAM = BiocParallel::SerialParam()))
    for (conv in c(TRUE, FALSE)) {
        ref <- .q(frag4feature(pa, x, convert2RawRT = conv))
        got <- .q(frag4feature(pa, xf, convert2RawRT = conv))
        expect_identical(.links(got), .links(ref), label = paste("convert", conv))
        expect_false(any(got@grped_df$is_filled))
        all <- .q(frag4feature(pa, xf, convert2RawRT = conv,
                               includeFilled = TRUE))@grped_df
        expect_true(all(all$is_filled[!paste(all$grpid, all$pid, all$cid) %in%
                                       paste(got@grped_df$grpid, got@grped_df$pid,
                                             got@grped_df$cid)]))
    }
    # With conversion, filled peaks are reported with their raw bounds too
    all <- .q(frag4feature(pa, xf, includeFilled = TRUE))@grped_df
    expect_true(any(all$is_filled))
    expect_false(anyNA(all$rtmin[all$is_filled]))
})

test_that("files are linked by name, and files the XCMS object lacks are reported", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xset.rds")
    ref <- .q(frag4feature(pa, x))
    win <- x
    win@filepaths <- paste0("C:\\data\\", basename(x@filepaths))
    expect_identical(.links(.q(frag4feature(pa, win))), .links(ref))
    other <- x
    other@filepaths[2] <- file.path(dirname(x@filepaths[2]), "other.mzML")
    expect_warning(res <- suppressMessages(frag4feature(pa, other)),
                   "does not hold LCMSMS_2.mzML")
    expect_false(any(res@grped_df$fileid == 2))
    none <- x
    none@filepaths <- c("/a/one.mzML", "/a/two.mzML")
    expect_error(suppressMessages(frag4feature(pa, none)),
                 "None of the purityA files")
})

test_that("grouped mode honours intense and falls back when iMz is missing", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    tab <- purityTable(pa)
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    g <- .q(frag4feature(pa, x, useGroup = TRUE))@grped_df
    i <- match(g$pid, tab$pid)
    expect_true(all(g$precurMzSource[!is.na(tab$iMz[i])] == "iMz"))
    expect_identical(g$precurMtchMZ[!is.na(tab$iMz[i])],
                     tab$iMz[i][!is.na(tab$iMz[i])])
    g <- .q(frag4feature(pa, x, useGroup = TRUE, intense = FALSE))@grped_df
    i <- match(g$pid, tab$pid)
    expect_true(all(g$precurMzSource %in% c("aMz", "precursorMZ")))
    s <- .f4f(".f4f_scans")(data.frame(pid = 1:3, fileid = 1, seqNum = 1:3,
        inPurity = 1, precursorScanNum = 1:3, precursorRT = c(1, NA, 3),
        retentionTime = c(1.5, 2.5, 3.5), iMz = c(100, NA, NA),
        aMz = c(101, 201, NA), precursorMZ = c(102, 202, 302)),
        intense = TRUE, rtSource = "precursor")
    expect_identical(s$mz1, c(100, 201, 302))
    expect_identical(s$precurMzSource, c("iMz", "aMz", "precursorMZ"))
    expect_identical(s$rtRaw, c(1, 2.5, 3))
})

test_that("every scan has one best match, and rtTol and groupWindow change the windows", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    for (mode in c(FALSE, TRUE)) {
        g <- .q(frag4feature(pa, x, useGroup = mode))@grped_df
        expect_identical(as.vector(table(g$pid[g$matchRank == 1L])),
                         rep(1L, length(unique(g$pid))))
        wide <- .q(frag4feature(pa, x, useGroup = mode, rtTol = 5))
        expect_true(all(.links(.q(frag4feature(pa, x, useGroup = mode))) %in%
                        .links(wide)))
        expect_gt(nrow(wide@grped_df), nrow(g))
    }
    full <- .q(frag4feature(pa, x, useGroup = TRUE))@grped_df
    med <- .q(frag4feature(pa, x, useGroup = TRUE, groupWindow = "median"))@grped_df
    expect_true(all(med$rtminWin >= med$rtmin_full & med$rtmaxWin <= med$rtmax_full))
    expect_lte(nrow(med), nrow(full))
})

test_that("a Frag4FeatureParam with the RT settings links as the scalar call", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    param <- Frag4FeatureParam(useGroup = TRUE, rtTol = 2,
                               rtSource = "ms2", groupWindow = "median")
    scalar <- .q(frag4feature(pa, x, useGroup = TRUE, rtTol = 2,
                              rtSource = "ms2", groupWindow = "median"))
    typed <- .q(frag4feature(pa, x, param = param))
    expect_gt(nrow(typed@grped_df), 0L)
    expect_identical(typed@grped_df, scalar@grped_df)
    record <- tail(processingParameters(typed, "frag4feature"), 1L)[[1L]]
    expect_identical(record[names(as.list(param))], as.list(param))
    expect_error(frag4feature(pa, x, param = param, rtTol = 3),
                 "rtTol conflicts")
})

test_that("MS3 and MS4 rows keep their own precursor m/z and time", {
    skip_if_not_installed("xcms")
    prepared <- .msn_standard_feature("tryptophan", "MSn")
    pa <- .q(purityA(.msn_standard_file("tryptophan", "MSn"), cores = 1L,
                     plotP = FALSE, msLevel = 2:4,
                     parentReference = "acquisitionNum"))
    g <- .q(frag4feature(pa, prepared$xe, ppm = 10, msLevel = 2:4,
                         parentReference = "acquisitionNum"))@grped_df
    child <- g[g$msLevel > 2L, , drop = FALSE]
    expect_gt(nrow(child), 0L)
    expect_true(all(child$precurMzSource == "precursorMZ"))
    # No retention time correction, so the mapped time is the precursor time.
    expect_equal(child$precurMtchRTadj, child$precurMtchRT)
    expect_true(all(g$precurMzSource[g$msLevel == 2L] %in%
                    c("iMz", "aMz", "precursorMZ")))
})

test_that("MS/MS files outside the XCMS object are mapped by ms2RtMap or estimated", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    x1 <- filterFile(x, 1, keepAdjustedRtime = TRUE, keepFeatures = TRUE)
    f2 <- function(res) {
        g <- res@grped_df[res@grped_df$fileid == 2, ]
        sort(paste(g$grpid, g$pid))
    }
    # File 2 is not in x1: its times are compared unmapped, with a warning
    expect_warning(ref <- suppressMessages(frag4feature(pa, x1, useGroup = TRUE)),
                   "No retention time map for LCMSMS_2.mzML")
    expect_true(all(ref@grped_df$rtMapStatus[ref@grped_df$fileid == 2] == "missing"))
    expect_gt(length(f2(ref)), 0)

    # 30 s of drift loses most links of file 2
    late <- .shift_pa(pa, 2, 30)
    lost <- .q(frag4feature(late, x1, useGroup = TRUE))
    expect_lt(length(intersect(f2(lost), f2(ref))), length(f2(ref)) / 2)

    # A constant offset recovers them
    back <- .q(frag4feature(late, x1, useGroup = TRUE,
                            ms2RtMap = c(LCMSMS_2.mzML = -30)))
    expect_identical(f2(back), f2(ref))
    expect_true(all(back@grped_df$rtMapStatus[back@grped_df$fileid == 2] == "external"))
    fn <- .q(frag4feature(late, x1, useGroup = TRUE,
                          ms2RtMap = list(LCMSMS_2.mzML = function(t) t - 30)))
    expect_identical(f2(fn), f2(ref))

    # So does an estimate from the precursors
    est <- .q(frag4feature(late, x1, useGroup = TRUE, estimateRtMap = TRUE))
    d <- est@params$frag4feature$rt_diagnostics[[2]]
    expect_identical(d$status, "estimated")
    expect_gt(d$landmarks, 2L)
    expect_gt(length(intersect(f2(est), f2(ref))), 0.8 * length(f2(ref)))
    expect_identical(d$estimate_method, "median")

    # A drift that grows with time is followed as well
    drift <- .shift_pa(pa, 2, 0)
    sp <- drift@spectra
    f <- as.character(sp$fileid) == "2"
    for (v in c("precursorRT", "retentionTime")) {
        z <- sp[[v]]
        z[f] <- z[f] + 20 + 0.2 * z[f]
        sp[[v]] <- z
    }
    drift@spectra <- sp
    est <- .q(frag4feature(drift, x1, useGroup = TRUE, estimateRtMap = TRUE))
    expect_identical(est@params$frag4feature$rt_diagnostics[[2]]$estimate_method,
                     "loess")
    expect_gt(length(intersect(f2(est), f2(ref))), 0.8 * length(f2(ref)))

    # An XCMS object holding the file supplies its adjusted times
    obj <- .q(frag4feature(pa, x1, useGroup = TRUE, ms2RtMap = x))
    g <- obj@grped_df[obj@grped_df$fileid == 2, ]
    expect_true(all(g$rtMapStatus == "external"))
    raw <- rtime(x, adjusted = FALSE, bySample = TRUE)[[2]]
    adj <- rtime(x, adjusted = TRUE, bySample = TRUE)[[2]]
    k <- match(g$precurMtchRT, raw)
    expect_identical(g$precurMtchRTadj[!is.na(k)], unname(adj[k[!is.na(k)]]))
})

test_that("diagnostics are recorded for each file", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    res <- .q(frag4feature(pa, .golden_xcms("msms_only_xcmsnexp.rds")))
    d <- res@params$frag4feature$rt_diagnostics
    expect_length(d, 2L)
    expect_identical(vapply(d, `[[`, character(1), "status"), c("aligned", "aligned"))
    expect_identical(vapply(d, `[[`, integer(1), "linked"),
                     vapply(1:2, function(i)
                         length(unique(res@grped_df$pid[res@grped_df$fileid == i])),
                         integer(1)))
    expect_true(all(abs(vapply(d, `[[`, numeric(1), "median_offset")) < 10))
    expect_silent(jsonlite::toJSON(d, auto_unbox = TRUE))
})

test_that("plim ignores missing purity and deprecated arguments accept data frames", {
    skip_if_no_mspuritydata()
    suppressPackageStartupMessages(library(xcms))
    pa <- .golden_pa("1_purityA_pa.rds")
    x <- .golden_xcms("msms_only_xcmsnexp.rds")
    pa <- .f4f(".pa_update")(pa)
    sp <- pa@spectra
    p <- sp$inPurity
    p[seq(1, length(p), by = 3)] <- NA
    sp$inPurity <- p
    pa_na <- pa
    pa_na@spectra <- sp
    g <- .q(frag4feature(pa_na, x, plim = 0.5))@grped_df
    expect_false(anyNA(g$pid))
    expect_true(all(g$inPurity > 0.5))
    expect_error(res <- .q(frag4feature(pa, x,
                     grp_peaklist = data.frame(a = 1))), NA)
    expect_s4_class(res, "purityA")
})
