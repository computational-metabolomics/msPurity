.ancestry_fixture <- function() {
    data.frame(pid = 1:6, spectrum_id = paste0("source-", 1:6),
        source_spectrum_id = paste0("source-", 1:6),
        run_id = rep(c("one", "two"), each = 3L),
        ms_level = rep(2:4, 2L), acquisition_num = rep(101:103, 2L),
        scan_index = rep(1:3, 2L), rtime = rep(c(10, 11, 12), 2L),
        precursor_mz = rep(c(500, 150, 75), 2L), polarity = 1L,
        parent_reference_raw = rep(c(NA, "101", "102"), 2L),
        stringsAsFactors = FALSE)
}

.prepare_ancestry <- function(...) msPurity:::.pa_msn_ancestry(...)
.propagate_links <- function(...) msPurity:::.pa_msn_feature_links(...)

test_that("metadata ancestry resolves parents within their source run", {
    x <- .prepare_ancestry(.ancestry_fixture())
    expect_identical(x$parentReference, "acquisitionNum")
    expect_identical(x$scans$parent_pid, c(NA_integer_, 1L, 2L,
            NA_integer_, 4L, 5L))
    expect_identical(x$scans$root_pid, c(1L, 1L, 1L, 4L, 4L, 4L))
    expect_identical(x$scans$path_status,
        rep(c("root", "valid", "valid"), 2L))
    expect_identical(x$scans$path[[6L]], 4:6)
})

test_that("unknown metadata keys cannot resolve a literal NA reference", {
    metadata <- .ancestry_fixture()[1:2, ]
    metadata$acquisition_num <- NA_real_
    metadata$scan_index <- NA_real_
    metadata$parent_reference_raw[2L] <- "NA"
    out <- .prepare_ancestry(metadata, "acquisitionNum")$scans
    expect_identical(out$parent_status[2L], "unresolved_reference")
    expect_true(is.na(out$root_pid[2L]))
})

test_that("explicit source and position conventions resolve their own keys", {
    metadata <- .ancestry_fixture()
    metadata$parent_reference_raw <- c(NA, "source-1", "source-2",
        NA, "source-4", "source-5")
    source <- .prepare_ancestry(metadata, "sourceSpectrumId")
    metadata$parent_reference_raw <- rep(c(NA, "1", "2"), 2L)
    position <- .prepare_ancestry(metadata, "scanIndex")
    expect_identical(source$scans$parent_pid, position$scans$parent_pid)
    expect_identical(position$scans$root_pid, c(1L, 1L, 1L, 4L, 4L, 4L))
})

test_that("ambiguous conventions and duplicate source keys are rejected", {
    metadata <- .ancestry_fixture()[1:3, ]
    metadata$ms_level <- c(2L, 2L, 3L)
    metadata$acquisition_num <- c(2L, 1L, 3L)
    metadata$parent_reference_raw <- c(NA, NA, "1")
    expect_error(.prepare_ancestry(metadata), "Ambiguous parent")
    expect_identical(.prepare_ancestry(metadata,
            "acquisitionNum")$scans$parent_pid[3L], 2L)
    metadata$acquisition_num <- c(1L, 1L, 3L)
    expect_error(.prepare_ancestry(metadata, "acquisitionNum"),
        "Duplicate acquisitionNum")
    metadata$spectrum_id[2L] <- metadata$spectrum_id[1L]
    expect_error(.prepare_ancestry(metadata, "sourceSpectrumId"),
        "source spectrum identities must be unique")
})

test_that("automatic references compare only valid parent interpretations", {
    metadata <- .ancestry_fixture()[1:3, ]
    metadata$ms_level <- c(2L, 2L, 3L)
    metadata$acquisition_num <- 1:3
    metadata$source_spectrum_id <- as.character(1:3)
    metadata$spectrumId <- paste0("scan=", 1:3)
    metadata$spectrumRef <- c(NA, NA, "scan=1")
    metadata$parent_reference_raw <- c(NA, NA, "1")
    agreed <- .prepare_ancestry(metadata)
    expect_identical(agreed$scans$parent_pid, c(NA_integer_, NA_integer_, 1L))
    expect_identical(agreed$parentReference, "spectrumId")

    metadata$scan_index <- c(2L, 1L, 3L)
    metadata$rtime[2L] <- 14
    invalid_competitor <- .prepare_ancestry(metadata)
    expect_identical(invalid_competitor$scans$parent_pid,
        agreed$scans$parent_pid)
    metadata$rtime[2L] <- 11
    expect_error(.prepare_ancestry(metadata), "Ambiguous parent")
})

test_that("invalid and incomplete paths remain explicit and unlinked", {
    metadata <- .ancestry_fixture()
    metadata$parent_reference_raw[2L] <- "999"
    metadata$parent_reference_raw[6L] <- "101"
    out <- .prepare_ancestry(metadata, "acquisitionNum")$scans
    expect_identical(out$path_status,
        c("root", "unresolved_reference", "unresolved_reference",
            "root", "valid", "missing_stage"))
    metadata <- .ancestry_fixture()
    metadata$polarity[2L] <- 0L
    expect_identical(.prepare_ancestry(metadata,
            "acquisitionNum")$scans$path_status[2:3],
        rep("invalid_reference", 2L))
    metadata <- .ancestry_fixture()
    metadata$rtime[2L] <- 9
    expect_identical(.prepare_ancestry(metadata,
            "acquisitionNum")$scans$path_status[2L], "invalid_reference")
    metadata <- .ancestry_fixture()
    metadata$scan_index[2L] <- 0L
    expect_identical(.prepare_ancestry(metadata,
            "acquisitionNum")$scans$path_status[2L], "invalid_reference")
    for (mass in c(NA_real_, 0, -1)) {
        metadata <- .ancestry_fixture()
        metadata$precursor_mz[2L] <- mass
        expect_identical(.prepare_ancestry(metadata,
                "acquisitionNum")$scans$path_status[2:3],
            rep("missing_precursor", 2L), info = paste("precursor:", mass))
    }
})

test_that("each root association occurrence reaches its eligible descendants", {
    ancestry <- .prepare_ancestry(.ancestry_fixture(), "acquisitionNum")
    links <- data.frame(pid = c(1L, 1L, 1L, 4L),
        grpid = c("first", "first", "second", "other"),
        inPurity = c(0.5, 0.5, 0.5, 0.95), mzmed = c(501, 501, 502, 510))
    out <- .propagate_links(links, ancestry, plim = 0.7)
    expect_identical(out$pid, c(4L, 2L, 2L, 2L, 3L, 3L, 3L, 5L, 6L))
    expect_identical(out$grpid[2:4], c("first", "first", "second"))
    expect_equal(out$featureMz[2:4], c(501, 501, 502))
    expect_equal(out$precurMtchMZ[2:4], rep(150, 3L))
    expect_equal(out$precurMtchRT[2:4], rep(10, 3L))
    expect_equal(out$retentionTime[2:4], rep(11, 3L))
    expect_true(all(is.na(out$inPurity[out$msLevel > 2L])))
    expect_identical(out$msPurityRootPid, c(4L, rep(1L, 6L), 4L, 4L))
    expect_false(1L %in% out$pid)
    expect_identical(.propagate_links(links[FALSE, ], ancestry)$pid, integer())
})

test_that("child acquisition fields come from the child, not its root link", {
    metadata <- .ancestry_fixture()[1:3, ]
    metadata$spectrumId <- paste0("scan=", metadata$acquisition_num)
    metadata$nativeId <- c("scan=101", "scan=102", NA_character_)
    metadata$filterString <- c("FTMS ms2", "ITMS ms3", NA_character_)
    metadata$instrumentConfigurationId <- c("IC1", "IC2", NA_character_)
    metadata$isolationWindowTargetMz <- c(500, 150, 75)
    metadata$analyzer <- c("Orbitrap", "ion trap", NA_character_)
    metadata$activationMethod <- c("HCD", "CID", NA_character_)
    metadata$collisionEnergy <- c(30, 35, NA_real_)
    metadata$collisionEnergyUnit <- c("eV", "eV", NA_character_)
    metadata$grpid <- "source-annotation"
    metadata$featureMz <- 999
    links <- data.frame(pid = 1L, grpid = "feature", mzmed = 501,
        sample = "sample", cid = 19L,
        inPurity = 0.5, aPurity = 0.4, iPurity = 0.6,
        purity_pass_flag = FALSE, precursorNearest = 99L,
        acquisitionNum = 101L, precursorMZ = 500, precursorRT = 9,
        precursorIntensity = 1000, activationType = "HCD",
        nativeId = "scan=101", filterString = "FTMS ms2",
        instrumentConfigurationId = "IC1", isolationWindowTargetMz = 500)
    ancestry <- .prepare_ancestry(metadata, "acquisitionNum")
    out <- .propagate_links(links, ancestry, plim = 1)
    expected <- list(pid = 2:3, nativeId = c("scan=102", NA_character_),
        filterString = c("ITMS ms3", NA_character_),
        instrumentConfigurationId = c("IC2", NA_character_),
        isolationWindowTargetMz = c(150, 75), acquisitionNum = 102:103,
        spectrumId = c("scan=102", "scan=103"),
        spectrum_id = c("source-2", "source-3"),
        precursorMZ = c(150, 75), precursorRT = c(10, 11),
        massAnalyzer = c("ion trap", NA_character_),
        activationType = c("CID", NA_character_),
        collisionEnergy = c(35, NA_real_), grpid = rep("feature", 2L),
        featureMz = rep(501, 2L), cid = rep(19L, 2L))
    expect_identical(as.list(out[names(expected)]), expected)
    expect_true(all(is.na(out[c("inPurity", "aPurity", "iPurity",
                        "purity_pass_flag", "precursorNearest",
                        "precursorIntensity")])))
})

test_that("Spectra preparation keeps dataset keys and reads no peak arrays", {
    object <- Spectra::addProcessing(.ancestry_spectra(), function(x) {
        stop("Peak access must not be needed for ancestry")
    })
    metadata <- msPurity:::.pa_msn_spectra_metadata(object)
    out <- .prepare_ancestry(metadata)
    expect_identical(metadata$pid, c(7L, 40L, 41L))
    expect_identical(metadata$spectrum_id, as.character(101:103))
    expect_identical(metadata$acquisition_num, 11:13)
    expect_identical(metadata$source_index, 1:3)
    expect_identical(metadata$massAnalyzer, rep("Orbitrap", 3L))
    expect_identical(out$parentReference, "sourceSpectrumId")
    expect_identical(out$scans$parent_pid, c(NA_integer_, 7L, 40L))
    object$precScanNum <- c(NA_integer_, 11L, 12L)
    out <- .prepare_ancestry(msPurity:::.pa_msn_spectra_metadata(object))
    expect_identical(out$parentReference, "acquisitionNum")
    expect_identical(out$scans$parent_pid, c(NA_integer_, 7L, 40L))
    empty <- Spectra::Spectra(Spectra::MsBackendMemory())
    prepared <- .prepare_ancestry(msPurity:::.pa_msn_spectra_metadata(empty))
    expect_identical(nrow(prepared$scans), 0L)
    expect_identical(prepared$scans$pid, integer())
})

test_that("native Parquet keys do not imply acquisition-number references", {
    skip_if_not_installed("MsBackendParquet")
    skip_if_no_arrow()
    paths <- character()
    on.exit(unlink(paths, recursive = TRUE), add = TRUE)
    for (convention in c("sourceSpectrumId", "acquisitionNum")) {
        references <- if (convention == "sourceSpectrumId") c(NA, 1L, 2L) else
            c(NA, 11L, 12L)
        metadata <- S4Vectors::DataFrame(pid = c(7L, 40L, 41L),
            msLevel = 2:4, scanIndex = 1:3, acquisitionNum = 11:13,
            rtime = c(10, 11, 12), precursorMz = c(500, 150, 75),
            precScanNum = references, polarity = 1L,
            dataOrigin = "native-fixture.mzML")
        metadata$mz <- rep(list(100), 3L)
        metadata$intensity <- rep(list(10), 3L)
        path <- tempfile("ancestry-native-")
        paths <- c(paths, path)
        MsBackendParquet::createMsBackendParquetDataset(path = path,
            data = metadata)
        object <- Spectra::Spectra(path,
            source = MsBackendParquet::MsBackendParquet())
        source <- msPurity:::.pa_msn_spectra_metadata(object)
        out <- .prepare_ancestry(source)
        expect_identical(source$spectrum_id, as.character(1:3),
            info = convention)
        expect_equal(source$acquisition_num, 11:13, info = convention)
        expect_identical(out$parentReference, convention, info = convention)
        expect_equal(out$scans$parent_pid, c(NA, 7, 40), info = convention)
    }
})

test_that("native parent references stay distinct from backend keys", {
    object <- .ancestry_spectra()
    object$spectrumId <- paste0("scan=", object$acquisitionNum)
    object$spectrumRef <- c(NA_character_, "scan=11", "scan=12")
    object$precScanNum <- c(NA_integer_, 11L, 12L)
    metadata <- msPurity:::.pa_msn_spectra_metadata(object)
    native <- .prepare_ancestry(metadata, "spectrumId")
    automatic <- .prepare_ancestry(metadata)

    expect_identical(native$scans$parent_pid, c(NA_integer_, 7L, 40L))
    expect_identical(automatic$parentReference, "spectrumId")
    expect_identical(automatic$scans$parent_pid, native$scans$parent_pid)
    expect_identical(metadata$source_spectrum_id, as.character(101:103))
    expect_identical(metadata$spectrumId, paste0("scan=", 11:13))
    selected <- Spectra::filterPrecursorScan(object, acquisitionNum = 12L)
    expect_identical(selected$pid, native$scans$pid)

    reordered <- .prepare_ancestry(metadata[c(3, 1, 2), ])$scans
    expect_identical(reordered$parent_pid, c(40L, NA_integer_, 7L))
    object$precScanNum[3L] <- 11L
    metadata <- msPurity:::.pa_msn_spectra_metadata(object)
    expect_error(.prepare_ancestry(metadata), "Ambiguous parent")
    expect_error(.prepare_ancestry(metadata, "acquisitionNum"),
        "Native spectrumRef conflicts")
})

test_that("unresolved native references do not silently use numeric parents", {
    object <- .ancestry_spectra()
    object$spectrumId <- paste0("scan=", object$acquisitionNum)
    object$spectrumRef <- c(NA_character_, "absent", "scan=12")
    object$precScanNum <- c(NA_integer_, 11L, 12L)
    metadata <- msPurity:::.pa_msn_spectra_metadata(object)
    result <- .prepare_ancestry(metadata)$scans
    expect_identical(result$parent_status[2L], "unresolved_reference")
    expect_true(is.na(result$root_pid[3L]))
    metadata$spectrumId[2L] <- metadata$spectrumId[1L]
    expect_error(.prepare_ancestry(metadata), "Duplicate spectrumId")
})

# Two runs of one Parquet dataset: spectrum_id_ numbers both runs in one
# sequence, precScanNum names the precursor's acquisitionNum. In run b the
# MS3 scan's precScanNum (4) is the acquisitionNum of its MS2 parent and,
# by coincidence, also the spectrum_id_ of the run's MS1 scan.
.two_run_dataset <- function() {
    metadata <- S4Vectors::DataFrame(
        pid = 1:7, spectrum_id_ = 1:7,
        run_id = rep(c("a", "b"), c(3L, 4L)),
        msLevel = c(1L, 2L, 3L, 1L, 2L, 2L, 3L),
        scanIndex = c(1:3, 1:4),
        acquisitionNum = c(1L, 2L, 3L, 1L, 4L, 5L, 6L),
        rtime = c(1, 2, 3, 1, 2, 3, 4),
        precursorMz = c(NA, 200, 100, NA, 200, 210, 100),
        precScanNum = c(NA, 1L, 2L, NA, 1L, 1L, 4L), polarity = 1L,
        dataOrigin = rep(c("a.mzML", "b.mzML"), c(3L, 4L)))
    metadata$mz <- rep(list(100), 7L)
    metadata$intensity <- rep(list(10), 7L)
    Spectra::Spectra(metadata, source = Spectra::MsBackendMemory())
}

test_that("dataset keys are not read as parent scan numbers", {
    metadata <- msPurity:::.pa_msn_spectra_metadata(.two_run_dataset())
    expect_true(all(metadata$source_key_assigned))
    x <- .prepare_ancestry(metadata)
    expect_identical(x$parentReference, "acquisitionNum")
    expect_identical(x$scans$parent_pid, c(NA, 1L, 2L, NA, 4L, 4L, 5L))
    expect_identical(x$scans$root_pid[7L], 5L)
    ## The reading by key stays available when asked for.
    by_key <- .prepare_ancestry(metadata, "sourceSpectrumId")
    expect_identical(by_key$scans$parent_pid[7L], 4L)
    ## Keys a caller supplies may be what precScanNum names, so a
    ## coincidence with acquisitionNum remains ambiguous.
    metadata$source_key_assigned <- FALSE
    expect_error(.prepare_ancestry(metadata), "Ambiguous parent")
})
