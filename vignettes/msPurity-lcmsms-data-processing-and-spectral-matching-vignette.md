---
title: "LC-MS/MS data processing and spectral matching workflow using msPurity and XCMS"
author: "Thomas N. Lawson"
date: "2026-11-20"
output: 
  BiocStyle::html_document:
    toc: true
    toc_depth: 3  
    number_sections: true  
    toc_float: true
    
bibliography: mspurity.bib
vignette: >
  %\VignetteIndexEntry{msPurity spectral matching}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
link-citations: true
---

# LC-MS/MS data processing and spectral matching workflow

## Read the two DDA runs as a Parquet study

This workflow links acquired MS2 spectra to chromatographic features,
averages the fragmentation spectra for each feature, and searches the same
queries and references with two scoring methods. The runs come from `msPurityData`, and the reference library ships with `msPurity`.


``` r
library(msPurity)
library(xcms)
```


``` r
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
```

## Detect chromatographic features with xcms 4

Limit peak detection to MS1 scans between 30 and 90 seconds and ions from
100 to 200 *m/z*. Keep the original spectra for precursor purity assessment
and MS2 processing. Each sample links to its scans by their common source
file name.


``` r
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
```

```
##          mzmed    rtmed npeaks
## FT001 100.0760 51.44504      1
## FT002 101.0035 60.76616      1
## FT003 101.0714 40.33178      1
## FT004 102.0553 36.42353      2
## FT005 102.0634 30.95213      1
## FT006 104.1073 50.08379      1
```

## Build the query and reference datasets

`purityA()` measures precursor ion purity for each fragmentation event
[@Lawson2017]. `frag4feature()` links the events to the chromatographic peaks
in their corresponding runs. Assess only MS2 events in the feature detection region, retaining all MS1
scans so each event keeps its original isolation-window context. Filter the
fragments and average the scans for
each feature across runs. The results dataset retains links to the source
measurements. A warning flags scans whose recorded precursor lies outside the
isolation bounds. Their purity is missing.


``` r
rt <- Spectra::rtime(raw)
precursor <- Spectra::precursorMz(raw)
selected_ms2 <- rt >= 30 & rt <= 90 & precursor >= 100 & precursor <= 200
scans <- raw[which(Spectra::msLevel(raw) == 1L | selected_ms2)]
pa <- purityA(scans, cores = 1L)
```

```
## Warning: 9 MS2 precursors are outside their recorded isolation bounds; purity
## is unavailable. Supply offsets when known.
## Warning: 9 MS2 precursors are outside their recorded isolation bounds; purity
## is unavailable. Supply offsets when known.
```

``` r
pa <- frag4feature(pa, xcmsObj = features)
pa <- filterFragSpectra(pa)
pa <- averageAllFragSpectra(pa)
results <- createDatabase(pa, features, format = "parquet", study = study,
                          outDir = work, dbName = "results")
```

Convert the bundled MoNA MSP file once. With your own data, supply the path
to the appropriate reference library.


``` r
msp <- system.file("extdata", "tests", "library", "mini_mona.msp",
                   package = "msPurity", mustWork = TRUE)
reference <- convertLibraryToParquet(msp, file.path(work, "library"),
                                     format = "msp", dialect = "mona")
```

## Rank library hits with the default score

Select the pooled MS2 averages (`av_all`). The default search considers
precursor-compatible library spectra, aligns fragment peaks [@Zhou2014],
and ranks hits within each query by dot product cosine (`dpc`). With the
default weights each aligned peak contributes $w = I^{0.5}(m/z)^2$
[@Horai2008]. The results also include a reverse dot product that discounts
unmatched query peaks and a composite dot product that incorporates
neighbouring peak intensity ratios [@stein1994optimization].


``` r
aligned <- spectralMatching(results, reference, format = "parquet",
                            q_spectraTypes = "av_all")
aligned$matchedResults[c("qpid", "library_accession", "rank", "dpc", "mcount")]
```

```
##   qpid  library_accession rank       dpc mcount
## 1    1 CCMSLIB00000577898    1 0.8739497      1
## 2    2           ML005101    1 0.9341948      2
## 3    3 CCMSLIB00003740024    1 0.5071041      2
## 4   19           PR100407    1 0.9693564      2
```

`qpid` identifies the query spectrum. The library accession identifies the
candidate, and `mcount` counts matched fragment peaks. Rank 1 is the best
score for that query. A library hit provides a candidate annotation.

## Search the same spectra with BLINK

BLINK computes an approximate cosine from binned fragment spectra
[@Harwood2023]. The defaults square-root transform and normalise intensities,
bin fragment masses at 0.001 Da, and apply a blur tolerance of 0.01 Da.
Multiple ions can contribute to a bin neighbourhood, so the cosine can
exceed one and the matching ion count is an estimate.


``` r
blink <- spectralMatching(results, reference, format = "parquet",
                          q_spectraTypes = "av_all", param = BlinkMatchingParam(),
                          updateDb = TRUE)
blink$matchedResults[c("qpid", "library_accession", "rank", "blinkCosine",
                       "blinkMatchingIonCount")]
```

```
##   qpid  library_accession rank blinkCosine blinkMatchingIonCount
## 1    1 CCMSLIB00000577898    1   0.5663570                     1
## 2    2           ML005101    1   0.8858762                     2
## 3    3 CCMSLIB00003740024    1   0.2584305                     2
## 4   19           PR100407    1   0.7607654                     2
```

## Read the saved evidence

`updateDb = TRUE` records the BLINK search in the query results dataset.
Use its manifest to locate the evidence and score tables. Each evidence row
identifies a ranked candidate and links to its scores by `evidence_id_`.


``` r
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
```

```
##   evidence_id_ database_identifier rank       score_term score_value
## 1            1  CCMSLIB00000577898    1 MSPURITY:0000034   0.5663570
## 2            1  CCMSLIB00000577898    1 MSPURITY:0000035   1.0000000
## 3            2            PR100407    1 MSPURITY:0000034   0.7607654
## 4            2            PR100407    1 MSPURITY:0000035   2.0000000
## 5            3            ML005101    1 MSPURITY:0000034   0.8858762
## 6            3            ML005101    1 MSPURITY:0000035   2.0000000
```

``` r
activity <- tail(manifest$provenance, 1)[[1]]
activity$parameters$scoring$method
```

```
## [1] "BLINK"
```

# References
