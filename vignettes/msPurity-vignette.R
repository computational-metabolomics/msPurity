## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(msPurity)
library(xcms)


## ----study, message=FALSE-----------------------------------------------------
mzml_dir <- system.file("extdata", "lcms", "mzML", package = "msPurityData",
                        mustWork = TRUE)
files <- list.files(mzml_dir, pattern = "[.]mzML$", full.names = TRUE)
work <- tempfile("mspurity-")
dir.create(work)
study <- file.path(work, "study")
MsBackendParquet::mzMLToParquet(files, path = study, engine = "mzr",
                                verbose = FALSE)
raw <- Spectra::Spectra(Spectra::backendInitialize(
    MsBackendParquet::MsBackendParquet(), path = study))
is_dda <- startsWith(basename(Spectra::dataOrigin(raw)), "LCMSMS_")
dda <- raw[is_dda]


## ----acquired-purity, message=FALSE-------------------------------------------
rt <- Spectra::rtime(dda)
precursor <- Spectra::precursorMz(dda)
selected_ms2 <- rt >= 30 & rt <= 90 & precursor >= 100 & precursor <= 200
scans <- dda[which(Spectra::msLevel(dda) == 1L | selected_ms2)]
pa <- purityA(scans, cores = 1L)
purity <- purityTable(pa)
head(purity[c("filename", "seqNum", "precursorMZ", "precursorRT", "inPurity")])


## ----isolation-efficiency, message=FALSE--------------------------------------
origin <- Spectra::dataOrigin(scans)
weighted <- purityA(scans[origin == origin[1]], cores = 1L, iwNorm = TRUE,
                    iwNormFun = iwNormGauss(sdlim = 3, minOff = -0.5, maxOff = 0.5))
head(purityTable(weighted)[c("seqNum", "precursorMZ", "inPurity")])


## ----detect-features, message=FALSE-------------------------------------------
detect_features <- function(ms1, noise = 100) {
    origins <- unique(Spectra::dataOrigin(ms1))
    samples <- S4Vectors::DataFrame(sample_name = basename(origins), origin = origins)
    experiment <- MsExperiment::MsExperiment(sampleData = samples, spectra = ms1)
    experiment <- MsExperiment::linkSampleData(
        experiment, with = "sampleData.origin = spectra.dataOrigin")

    peaks <- findChromPeaks(experiment, param = CentWaveParam(
        snthresh = 3, noise = noise, ppm = 10, peakwidth = c(3, 30)),
        BPPARAM = BiocParallel::SerialParam())
    mz <- chromPeaks(peaks)[, "mz"]
    peaks <- filterChromPeaks(peaks, keep = mz >= 100 & mz <= 200)
    groupChromPeaks(peaks, param = PeakDensityParam(
        sampleGroups = rep("sample", nrow(samples)), minFraction = 0,
        bw = 5, binSize = 0.017))
}
ms1_window <- function(sp) {
    Spectra::filterRt(Spectra::filterMsLevel(sp, 1L), c(30, 90))
}
features <- detect_features(Spectra::filterMzRange(ms1_window(dda), c(100, 200)))
pa <- frag4feature(pa, xcmsObj = features)
linked <- groupedSpectra(pa, legacy = FALSE)
head(Spectra::spectraData(linked, c("precursorMz", "rtime")))


## ----average-and-export, message=FALSE----------------------------------------
pa <- filterFragSpectra(pa)
pa <- averageAllFragSpectra(pa)
length(averagedSpectra(pa, legacy = FALSE))
createMSP(pa, method = "all", msp_file_pth = file.path(work, "features.msp"))


## ----per-run-averages, eval=FALSE---------------------------------------------
# pa <- averageIntraFragSpectra(pa)
# pa <- averageInterFragSpectra(pa)


## ----parquet-results, message=FALSE-------------------------------------------
results <- createDatabase(pa, features, format = "parquet", study = study,
                          outDir = work, dbName = "results")
consensus <- readProcessedSpectra(results, msLevel = 2L, method = "all")
length(consensus)
msp <- system.file("extdata", "tests", "library", "mini_mona.msp",
                   package = "msPurity", mustWork = TRUE)
reference <- convertLibraryToParquet(msp, file.path(work, "library"),
                                     format = "msp", dialect = "mona")
search <- spectralMatching(results, reference, format = "parquet",
                           q_spectraTypes = "av_all")
head(search$matchedResults[c("qpid", "library_accession", "dpc", "rank")])


## ----predicted-purity, message=FALSE------------------------------------------
ms1_features <- detect_features(ms1_window(raw[!is_dda]), noise = 5e5)
area <- rowMeans(featureValues(ms1_features, value = "into"))
selected <- head(order(area, decreasing = TRUE, na.last = NA), 2)
px <- purityX(ms1_features, xgroups = selected, purityType = "purityFWmedian",
              offsets = c(0.5, 0.5), ilim = 0, cores = 1L)
px@predictions

