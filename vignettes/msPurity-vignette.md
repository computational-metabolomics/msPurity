---
title: "Using msPurity for Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics"
author: "Thomas N. Lawson"
date: "2018-11-12"
output: 
  BiocStyle::html_document:
    toc: true
bibliography: mspurity.bib
vignette: >
  %\VignetteIndexEntry{msPurity}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

# Introduction

The function of this R package is to assess the contribution of the targeted precursor in a fragmentation isolation window using a metric called “precursor purity”. 

What we call "Precursor purity" is a measure of the contribution of a selected precursor peak in an isolation window used for fragmentation. The simple calculation involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. When assessing MS/MS spectra this calculation is done before and after the MS/MS scan of interest and the purity is interpolated at the time of the MS/MS acquisition. The calculation is very similar to the "Precursor Ion Fraction" (PIF) metric described by [@Michalski2011] for proteomics with the exception that purity here is interpolated at the recorded point of MS2 acquisition using bordering full-scan spectra. Additionally, low abundance ions that are remove that are thought to have limited contribution to the resulting MS2 spectra and can optionally take into account the isolation efficiency of the mass spectrometer

There are two main use cases for the package

1.	Assessing precursor purity of previously acquired MS2 spectra: A user has acquired either LC-MS2 or DIMS2 spectra and an assessment is made of the precursor purity for each MS2 scan. `purityA`
2.	Assessing precursor purity of anticipated isolation windows for MS2 spectra: A user has acquired either LC-MS (`purityX`) or DIMS (`purityD`) full scan (MS1) data and an assessment is to be made of the precursor purity of detected features using anticipated or theoretical isolation windows. This information can then be used to guide further targeted MS2 experiments. 

The package has been developed to be used with DI-MS or LC-MS data and has been checked to work with the following vendor files after conversion to mzML: Thermo, Agilent and AB Sciex.

# Assessing precursor purity of previously acquired MS2 spectra

## purityA 

Given a vector of LC-MS/MS or DI-MS/MS mzML file paths the function `purityA` will calculate the precursor purity of each MS/MS scan. The output is a S4 class object where a dataframe of the purity results can be accessed using the appropriate slot (`@puritydf`).

The isolation widths will be determined automatically from the mzML file. For some mzML files this is not recorded and in these cases the offsets can be given as a parameter.

In the case of Agilent only the "narrow" isolation is supported. This roughly equates to +/- 0.65 Da (depending on the instrument). If the file is detected as originating from an Agilent instrument the isolation widths will automatically be set as +/- 0.65 Da (this can be overwritten with the `offsets` argument)

The purity dataframe (`pa@puritydf`) consists of the following columns:

* __pid__: unique id for MS/MS scan
* __fileid__: unqiue id for file
* __seqNum__: scan number
* __precursorIntensity__: precursor intensity value as defined from mzML file
* __precursorMZ__: precursor m/z value as defined from mzML file
* __precursorRT__: precursor RT value as defined from mzML file
* __precursorScanNum__: precursor scan number value as defined from mzML file
* __id__: unique id (redundant)
* __filename__: mzML filename
* __precursorNearest__: MS1 scan nearest to this MS/MS scan
* __aMz__: The _m/z_ value in the precursorNearest scan which most closely matches the precursorMZ value provided from the mzML file
* __aPurity__: The purity score for __aMz__ 
* __apkNm__: The number of peaks in the isolation window for __aMz__ 
* __iMz__: The _m/z_ value in the precursorNearest scan that is the most intense within the isolation window.
* __iPurity__: The purity score for __iMz__ 
* __ipkNm__: The number of peaks in the isolation window for __iMz__ 
* __inPurity__: The interpolated purity score
* __inpkNm__: The interpolated number of peaks in the isolation window




```r
library(msPurity)
msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
msPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "LCMS_")
```



























