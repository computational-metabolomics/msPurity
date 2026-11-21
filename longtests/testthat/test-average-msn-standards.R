.average_standard <- local({
    cache <- list()
    function(compound) {
        if (is.null(cache[[compound]])) {
            prepared <- .msn_standard_feature(compound, "MSn")
            stopifnot(length(prepared$feature) == 1L)
            pa <- purityA(.msn_standard_file(compound, "MSn"), cores = 1L,
                plotP = FALSE, msLevel = 2:4,
                parentReference = "acquisitionNum")
            pa <- frag4feature(pa, prepared$xe, ppm = 10, msLevel = 2:4,
                parentReference = "acquisitionNum", convert2RawRT = FALSE)
            cache[[compound]] <<- list(pa = pa, feature = prepared$feature)
        }
        cache[[compound]]
    }
})

for (compound in c("tryptophan", "kynurenine")) local({
    name <- compound
    test_that(paste(name, "selects and pools the requested energy ranges"), {
        skip_if_not_installed("xcms")
        standard <- .average_standard(name)
        energies <- list(CID = seq(10, 60, 10), HCD = seq(10, 80, 10))
        sets <- list()
        for (activation in names(energies)) {
            selection <- list(msLevel = 2L, activationType = activation,
                collisionEnergy = energies[[activation]])
            sets[[activation]] <- AverageFragSpectraParam(select = selection,
                minfrac = 0, ppm = 10, av = "mean", sumi = FALSE)
            sets[[paste0(activation, "_pooled")]] <-
                AverageFragSpectraParam(select = selection,
                    combine = "collisionEnergy", minfrac = 0, ppm = 10,
                    av = "mean", sumi = FALSE)
        }
        out <- averageIntraFragSpectra(standard$pa, param = sets)
        s <- averagedSpectra(out, legacy = FALSE)
        s <- s[s$grpid == as.character(standard$feature)]

        for (activation in names(energies)) {
            separate <- s[s$msPurityParameterSet == activation]
            pooled <- s[s$msPurityParameterSet == paste0(activation, "_pooled")]
            expect_equal(sort(separate$collisionEnergy), energies[[activation]],
                info = activation)
            expect_equal(length(pooled), 1L, info = activation)
            expected_links <- do.call(rbind,
                .average_selected_links(separate))$featureLinkIndex
            expect_identical(sort(.average_selected_links(pooled)[[1L]]$
                featureLinkIndex), sort(expected_links), info = activation)
            expect_true(nrow(averagePeakMembers(pooled, pooled$spectrumId)) > 0L,
                info = activation)
        }
    })

    test_that(paste(name, "MS2-MS4 links average through the public route"), {
        skip_if_not_installed("xcms")
        standard <- .average_standard(name)

        out <- averageAllFragSpectra(standard$pa, rmp = TRUE,
            param = AverageFragSpectraParam(minfrac = 0, ppm = 10))
        averages <- averagedSpectra(out, legacy = FALSE)
        averages <- averages[averages$grpid == as.character(standard$feature)]

        # Every feature link of the standard is selected by exactly one
        # average of its own MS level (link counts as in
        # test-msn-integration.R).
        selected <- .average_selected_links(averages)
        levels <- rep(Spectra::msLevel(averages), vapply(selected, nrow,
            integer(1)))
        links <- do.call(rbind, selected)$featureLinkIndex
        expected <- if (name == "tryptophan") c(17L, 11L, 6L) else
            c(34L, 21L, 12L)
        expect_identical(tabulate(levels - 1L, nbins = 3L), expected)
        expect_false(anyDuplicated(links) > 0L)
        expect_identical(sort(links),
            which(out@grped_df$grpid == standard$feature))

        nonempty <- lengths(Spectra::mz(averages)) > 0L
        for (i in seq_along(averages)) {
            used <- averages$msPurityContributors[[i]]$featureLinkIndex
            expect_true(all(used %in% selected[[i]]$featureLinkIndex),
                info = i)
            level <- Spectra::msLevel(averages)[i]
            if (nonempty[i] && level > 2L)
                expect_identical(nrow(averages$msPurityPrecursorPath[[i]]),
                    level - 2L, info = i)
        }
        expect_false("msPurityAverageStatus" %in%
            Spectra::spectraVariables(averages))
        expect_true(any(Spectra::msLevel(averages)[nonempty] == 4L))
    })
})

test_that("ion-trap MS4 peaks combine with an absolute tolerance", {
    skip_if_not_installed("xcms")
    # Kynurenine MS4 is measured in the ion trap, where repeat fragment m/z
    # values differ by about 0.03-0.15 Da; MS2 and MS3 are FTMS.
    standard <- .average_standard("kynurenine")
    # Ion-trap MS4 averages: their selected links, and the number of averaged
    # peaks combining peaks of at least two different spectra.
    ms4 <- function(itms) {
        out <- averageAllFragSpectra(standard$pa, param = list(
            ftms = AverageFragSpectraParam(minfrac = 0, ppm = 10,
                select = list(scanFilterAnalyzer = "FTMS")),
            itms = itms))
        averages <- averagedSpectra(out, legacy = FALSE)
        averages <- averages[averages$grpid ==
            as.character(standard$feature) &
            averages$msPurityParameterSet == "itms"]
        selected <- vapply(.average_selected_links(averages), function(x)
            paste(sort(x$featureLinkIndex), collapse = ","), character(1))
        clusters <- unlist(averages$msPurityClusterContributors,
            recursive = FALSE)
        list(selected = sort(selected),
            combined = sum(lengths(lapply(clusters, unique)) >= 2L))
    }
    select <- list(scanFilterAnalyzer = "ITMS")
    relative <- ms4(AverageFragSpectraParam(select = select, minfrac = 0,
        ppm = 10))
    absolute <- ms4(AverageFragSpectraParam(select = select, minfrac = 0,
        ppm = 0, tolerance = 0.2))
    # Both settings form the same averages, so only the peak tolerance differs.
    expect_identical(absolute$selected, relative$selected)
    expect_true(any(grepl(",", relative$selected)))
    expect_identical(relative$combined, 0L)
    expect_gt(absolute$combined, 0L)
})

test_that("one-file standards form the same groups at every level", {
    skip_if_not_installed("xcms")
    standard <- .average_standard("tryptophan")
    param <- AverageFragSpectraParam(minfrac = 0, ppm = 10)
    out <- averageAllFragSpectra(standard$pa, param = param)
    out <- averageIntraFragSpectra(out, param = param)
    out <- averageInterFragSpectra(out, param = param)
    averages <- averagedSpectra(out, legacy = FALSE)
    averages <- averages[averages$grpid == as.character(standard$feature)]
    members <- function(level) {
        x <- averages[averages$av_level == level]
        sort(vapply(.average_selected_links(x), function(links)
            paste(sort(links$featureLinkIndex), collapse = ","), ""))
    }

    # With one file, intra-file grouping equals all-spectrum grouping.
    expect_identical(members("av_intra"), members("av_all"))
    inter <- averages[averages$av_level == "av_inter"]
    expect_true(all(lengths(.average_input_ids(inter)) == 1L))
    expect_identical(length(inter), sum(averages$av_level == "av_intra" &
        lengths(Spectra::mz(averages)) > 0L))
})

test_that("grouped feature associations support all three averaging levels", {
    skip_if_not_installed("xcms")
    standard <- .average_standard("tryptophan")
    prepared <- .msn_standard_feature("tryptophan", "MSn")
    linked <- frag4feature(standard$pa, prepared$xe, ppm = 10,
        msLevel = 2:4, parentReference = "acquisitionNum",
        convert2RawRT = FALSE, useGroup = TRUE)
    expect_gt(nrow(linked@grped_df), 0L)
    expect_false(any(c("sample", "cid") %in% names(linked@grped_df)))
    param <- AverageFragSpectraParam(minfrac = 0, ppm = 10, sumi = FALSE)

    out <- averageAllFragSpectra(linked, param = param)
    out <- averageIntraFragSpectra(out, param = param)
    out <- averageInterFragSpectra(out, param = param)
    expect_identical(out@grped_df, linked@grped_df)
    averages <- averagedSpectra(out, legacy = FALSE)
    averages <- averages[averages$grpid == as.character(prepared$feature)]

    for (level in c("av_all", "av_intra", "av_inter")) {
        selected <- averages[averages$av_level == level]
        expect_true(setequal(Spectra::msLevel(selected), 2:4), info = level)
        expect_true(any(lengths(Spectra::mz(selected)) > 0L),
            info = level)
        if (level == "av_intra")
            expect_true(all(selected$sample == "1"), info = level)
        contributors <- do.call(rbind, selected$msPurityContributors)
        expect_true(all(is.na(contributors$sample)), info = level)
        expect_true(all(contributors$fileid == 1L), info = level)
        expect_false("cid" %in% names(contributors), info = level)
        target_links <- which(out@grped_df$grpid ==
            as.character(prepared$feature))
        expect_true(all(contributors$featureLinkIndex %in% target_links),
            info = level)
        if (level != "av_inter") {
            selected_links <- do.call(rbind,
                .average_selected_links(selected))$featureLinkIndex
            expect_identical(sort(selected_links), target_links, info = level)
        }
    }

    # A copy reordered before its first average still addresses the same
    # association rows and their original positions within this feature.
    g <- linked@grped_df
    linked@grped_df <- g[rev(target_links), , drop = FALSE]
    reordered <- averageAllFragSpectra(linked, param = param)
    selected <- do.call(rbind, .average_selected_links(
        averagedSpectra(reordered, legacy = FALSE)))
    rows <- match(selected$linkId, g$linkId)
    expect_setequal(rows, target_links)
    expect_identical(selected$featureLinkIndex, rows)
    expect_identical(selected$contributorIndex, match(rows, target_links))
    expect_identical(unique(selected$associationId),
        names(tail(processingHistory(linked, "frag4feature"), 1L)))
})
