## ----packages, message=FALSE, warning=FALSE-----------------------------------
library(msPurity)
library(xcms)


## ----study, message=FALSE-----------------------------------------------------
mzml_dir <- system.file("extdata", "lcms", "mzML", package = "msPurityData",
                        mustWork = TRUE)
files <- list.files(mzml_dir, pattern = "^LCMSMS_.*[.]mzML$", full.names = TRUE)
work <- tempfile("mspurity-")
dir.create(work)
study <- file.path(work, "study")
MsBackendParquet::mzMLToParquet(files, path = study, engine = "mzr",
                                verbose = FALSE)
raw <- Spectra::Spectra(Spectra::backendInitialize(
    MsBackendParquet::MsBackendParquet(), path = study))


## ----features, message=FALSE--------------------------------------------------
ms1 <- Spectra::filterMsLevel(raw, 1L)
ms1 <- Spectra::filterMzRange(Spectra::filterRt(ms1, c(30, 90)), c(100, 200))
origins <- unique(Spectra::dataOrigin(ms1))
samples <- S4Vectors::DataFrame(sample_name = basename(origins), origin = origins)
experiment <- MsExperiment::MsExperiment(sampleData = samples, spectra = ms1)
experiment <- MsExperiment::linkSampleData(
    experiment, with = "sampleData.origin = spectra.dataOrigin")

peaks <- findChromPeaks(experiment, param = CentWaveParam(
    snthresh = 3, noise = 100, ppm = 10, peakwidth = c(3, 30)),
    BPPARAM = BiocParallel::SerialParam())
features <- groupChromPeaks(peaks, param = PeakDensityParam(
    sampleGroups = rep("sample", nrow(samples)), minFraction = 0,
    bw = 5, binSize = 0.017))
head(featureDefinitions(features)[, c("mzmed", "rtmed", "npeaks")])


## ----query, message=FALSE-----------------------------------------------------
rt <- Spectra::rtime(raw)
precursor <- Spectra::precursorMz(raw)
selected_ms2 <- rt >= 30 & rt <= 90 & precursor >= 100 & precursor <= 200
scans <- raw[which(Spectra::msLevel(raw) == 1L | selected_ms2)]
pa <- purityA(scans, cores = 1L)
pa <- frag4feature(pa, xcmsObj = features)
pa <- filterFragSpectra(pa)
pa <- averageAllFragSpectra(pa)
results <- createDatabase(pa, features, format = "parquet", study = study,
                          outDir = work, dbName = "results")


## ----reference, message=FALSE-------------------------------------------------
msp <- system.file("extdata", "tests", "library", "mini_mona.msp",
                   package = "msPurity", mustWork = TRUE)
reference <- convertLibraryToParquet(msp, file.path(work, "library"),
                                     format = "msp", dialect = "mona")


## ----default-matching, message=FALSE------------------------------------------
aligned <- spectralMatching(results, reference, format = "parquet",
                            q_spectraTypes = "av_all")
aligned$matchedResults[c("qpid", "library_accession", "rank", "dpc", "mcount")]


## ----blink-matching, message=FALSE--------------------------------------------
blink <- spectralMatching(results, reference, format = "parquet",
                          q_spectraTypes = "av_all", param = BlinkMatchingParam(),
                          updateDb = TRUE)
blink$matchedResults[c("qpid", "library_accession", "rank", "blinkCosine",
                       "blinkMatchingIonCount")]


## ----saved-evidence-----------------------------------------------------------
manifest <- MsBackendParquet::readManifest(results)
tables <- manifest$results$tables
read_table <- function(name) {
    path <- file.path(results, tables[[name]]$path)
    as.data.frame(arrow::open_dataset(path))
}
evidence <- read_table("evidence")
scores <- read_table("evidence_score")
head(merge(evidence[c("evidence_id_", "database_identifier", "rank")],
           scores[c("evidence_id_", "score_term", "score_value")],
           by = "evidence_id_"))
activity <- tail(manifest$provenance, 1)[[1]]
activity$parameters$scoring$method

