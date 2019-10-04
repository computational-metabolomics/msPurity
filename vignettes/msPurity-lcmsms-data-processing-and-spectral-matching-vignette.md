---
title: "LC-MS/MS data processing and spectral matching workflow using msPurity and XCMS"
author: "Thomas N. Lawson"
date: "2019-09-05"
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
---

# LC-MS/MS data processing and spectral matching workflow

## Overview
The msPurity package can be used with XCMS as part of a data processing and annotation workflow for LC-MS/MS data

* Purity assessments
    +  (mzML files) -> purityA -> (pa)
* XCMS processing
    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
* Fragmentation processing
    + (xset, pa) -> frag4feature -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> spectralMatching -> (sqlite spectral database)


## XCMS processing
We first need to run XCMS so that we can later link the spectral matching result back to XCMS feature.

(Please use the appropiate settings for your data)


```r
library(msPurity)
mzMLpths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE)
xset <- xcms::xcmsSet(mzMLpths)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)
```


## Purity assesments and linking fragmentation to XCMS features

The `purityA` function is then called to calculate the precursor purity of the fragmentation results and the `frag4feature` function will link the 
fragmentation data back to the XCMS feature.


```r
pa  <- purityA(mzMLpths)
pa <- frag4feature(pa, xset)
```

## Filtering and averaging

The fragmentation can be filtered prior to averaging using the “filterFragSpectra” function


```r
pa <- filterFragSpectra(pa)
```


Averaging of the fragmentation spectra can be done with either “averageAllFragSpectra” or with “averageIntraFragSpectra” and averageInterFragSpectra”. This will depend if the user wishes to treat the fragmentation spectra from within a file and between files. Another alternative is to ignore the averaging completely and just use the non-averaged fragmentation spectra for the spectral matching.

If the inter and intra fragmentation scans are to be treated differently the following should be followed:


```r
pa <- averageIntraFragSpectra(pa) # use parameters specific to intra spectra 
pa <- averageInterFragSpectra(pa) # use parameters specific to inter spectra
```

If the inter and intra fragmentation scans are to be treated the same the following workflow should be used.


```r
pa <- averageAllFragSpectra(pa) 
```

## Creating a spectral-database

An SQLite database is then created of the LC-MS/MS experiment. The SQLite schema of the spectral database can be detailed [here](https://bioconductor.org/packages/release/bioc/vignettes/msPurity/inst/doc/msPurity-spectral-datatabase-schema.html).


```r
td <- tempdir()
q_dbPth <- createDatabase(pa, xset, outDir = td, dbName = 'lcmsms-processing.sqlite')
```

## Spectral matching

The spectralMatching function allows users to perform spectral matching to be performed for **Query** SQLite spectral-database against a **Library** SQLite spectral-database.

The query spectral-database in most cases should contain be the "unknown" spectra database generated the msPurity function createDatabase as part of a msPurity-XCMS data processing workflow.

The library spectral-database in most cases should contain the "known" spectra from either public or user generated resources. The library SQLite database by default contains data from MoNA including Massbank, HMDB, LipidBlast and GNPS. A larger database can be downloaded from here. To create a user generated library SQLite database the following tool can be used to generate a SQLite database from a collection of MSP files: msp2db. It should be noted though, that as long as the schema of the spectral-database is as described here, then any database can be used for either the library or query - even allowing for the same database to be used.

The spectral matching functionality has four main components, **spectral filtering**, **spectral alignment**, **spectral matching** and finally summarising the results.

**Spectral filtering** is simply filtering both the library and query spectra to be search against (e.g. choosing the library source, instrument, retention time, precursor PPM tolerance, xcms features etc).

The **spectral alignment** stage involves aligning the query peaks to the library peaks. The approach used is similar to modified pMatch algorithm described in [@Zhou2014].

The **spectral matching** of the aligned spectra is performed against a combined intensity and m/z weighted vector - created for both the query and library spectra (wq and wl). See below:

$$ \vec{w}=\vec{intensity}^x \cdot \vec{mz}^y$$

Where x and y represent weight factors, defaults to $x=0.5$ and $y=2$ as per MassBank recommendations for ESI based data [@Horai2008]. These can be adjusted by the user though.

The aligned weighted vectors are then matched using dot product cosine, reverse dot product cosine and the composite dot product. See below for dot product cosine equation.

$$ F_{dpc} = \frac{\sum \vec{w_{Q} }\cdot \vec{w_{L}}}{\sqrt{\sum \vec{w_{Q}^{2}}} \cdot \sqrt{\sum \vec{w_{L}^{2}}}} $$

The reverse dot product cosine (rpdc) uses the same algorithm as dpc but all peaks that do not match in the query spectra (based on the alignment) are omitted from the calculation. This will improve scores when the query spectra is noisy but should be used with caution as it might lead to more false positives.

The composite dot product cosine (cdpc) approach is also calculated - this appraoch is used in the NIST MS search tool and incorporates relative intensity of neighbouring peaks (see function $M_{rel}$ ), where $N$=number of peaks, $Q$=query, $L$=library, $L\&Q$= matching library and query peaks, $w$ is the weighted value and $n$ is either 1 (if the abundance ratio of the library, i.e. $\frac{w_{L,i}}{w_{L,i-1}}$, is $<$ than the abundance ratio of the query i.e. $\frac{w_{Q,i}}{w_{Q,i-1}}$) or -1 (if the abundance ratio of the library is $>$ than the abundance ratio of the query). The approach was first described in [@stein1994optimization].

$$ F_{rel} = \Bigg( \frac{1}{N_{L\&Q}-1} \Bigg) \cdot \sum_{i=2}^{N_{L\&Q}} \Bigg( \frac{w_{L,i}}{w_{L,i-1}} \Bigg)_{}{^n} \cdot \Bigg( \frac{w_{Q,i}}{w_{Q,i-1}} \Bigg)_{}{^{-n}}$$

$$ F_{cpdc} = \frac{1000}{N_{Q} + N_{L\&Q}} \cdot (N_{Q} \cdot F_{dpc} + N_{L\&Q} \cdot F_{rel}) $$

The following example shows how to match two xcms groups against two of the library spectral filtered by their MoNA/MassBank accession ids.


```r
result <- spectralMatching(q_dbPth, q_xcmsGroups = c(17, 41), l_accessions=c('CCMSLIB00000577898','CE000616'))
```

```
## Running msPurity spectral matching function for LC-MS(/MS) data
```

```
## Filter query dataset
```

```
## Filter library dataset
```

```
## aligning and matching
```

```
## Summarising LC feature annotations
```

```r
print(result)
```

```
## $q_dbPth
## [1] "/tmp/RtmpodCsKs/lcmsms-processing.sqlite"
## 
## $matchedResults
##    lpid qpid mid       dpc rdpc      cdpc mcount allcount   mpercent
## 1  5325 1661   1 0.8739497    1 0.8359519      1       22 0.04545455
## 2  5325 1662   2 0.8739497    1 0.8359519      1       22 0.04545455
## 3 53807 1664   3 0.9408905    1 0.8960862      1       20 0.05000000
## 4 53807 1665   4 0.9408905    1 0.8960862      1       20 0.05000000
##   library_rt query_rt rtdiff library_precursor_mz query_precursor_mz
## 1       <NA> 44.45066     NA               116.07           116.0705
## 2       <NA> 44.45066     NA               116.07           116.0705
## 3       <NA> 70.39686     NA            132.10191           132.1018
## 4       <NA> 70.39686     NA            132.10191           132.1018
##   library_precursor_ion_purity query_precursor_ion_purity
## 1                         <NA>                   0.997344
## 2                         <NA>                   1.000000
## 3                         <NA>                   1.000000
## 4                         <NA>                   1.000000
##    library_accession library_precursor_type library_entry_name
## 1 CCMSLIB00000577898                    M+H          L-PROLINE
## 2 CCMSLIB00000577898                    M+H          L-PROLINE
## 3           CE000616                 [M+H]+         Isoleucine
## 4           CE000616                 [M+H]+         Isoleucine
##                      inchikey library_source_name library_compound_name
## 1 ONIBWKKTOPOVIA-UHFFFAOYSA-N                gnps             L-PROLINE
## 2 ONIBWKKTOPOVIA-UHFFFAOYSA-N                gnps             L-PROLINE
## 3 AGPKZVBTJJNPAG-UHFFFAOYSA-N            massbank          L-ISOLEUCINE
## 4 AGPKZVBTJJNPAG-UHFFFAOYSA-N            massbank          L-ISOLEUCINE
## 
## $xcmsMatchedResults
##    pid grpid       mz    mzmin    mzmax       rt    rtmin    rtmax npeaks
## 1 1661    17 116.0705 116.0703 116.0706 44.45066 43.95639 44.90363      4
## 2 1662    17 116.0705 116.0703 116.0706 44.45066 43.95639 44.90363      4
## 3 1664    41 132.1018 132.1017 132.1018 70.39686 69.81528 70.75930      4
## 4 1665    41 132.1018 132.1017 132.1018 70.39686 69.81528 70.75930      4
##   mzML   LCMSMS_1   LCMSMS_2     LCMS_1     LCMS_2 grp_name  lpid mid
## 1    4  130337063  124086404  230819937  234755462  M116T44  5325   1
## 2    4  130337063  124086404  230819937  234755462  M116T44  5325   2
## 3    4 2879676651 2794252073 3455044747 3379259900  M132T70 53807   3
## 4    4 2879676651 2794252073 3455044747 3379259900  M132T70 53807   4
##         dpc rdpc      cdpc mcount allcount   mpercent library_rt query_rt
## 1 0.8739497    1 0.8359519      1       22 0.04545455       <NA> 44.45066
## 2 0.8739497    1 0.8359519      1       22 0.04545455       <NA> 44.45066
## 3 0.9408905    1 0.8960862      1       20 0.05000000       <NA> 70.39686
## 4 0.9408905    1 0.8960862      1       20 0.05000000       <NA> 70.39686
##   rtdiff library_precursor_mz query_precursor_mz
## 1     NA               116.07           116.0705
## 2     NA               116.07           116.0705
## 3     NA            132.10191           132.1018
## 4     NA            132.10191           132.1018
##   library_precursor_ion_purity query_precursor_ion_purity
## 1                         <NA>                   0.997344
## 2                         <NA>                   1.000000
## 3                         <NA>                   1.000000
## 4                         <NA>                   1.000000
##    library_accession library_precursor_type library_entry_name
## 1 CCMSLIB00000577898                    M+H          L-PROLINE
## 2 CCMSLIB00000577898                    M+H          L-PROLINE
## 3           CE000616                 [M+H]+         Isoleucine
## 4           CE000616                 [M+H]+         Isoleucine
##                      inchikey library_source_name library_compound_name
## 1 ONIBWKKTOPOVIA-UHFFFAOYSA-N                gnps             L-PROLINE
## 2 ONIBWKKTOPOVIA-UHFFFAOYSA-N                gnps             L-PROLINE
## 3 AGPKZVBTJJNPAG-UHFFFAOYSA-N            massbank          L-ISOLEUCINE
## 4 AGPKZVBTJJNPAG-UHFFFAOYSA-N            massbank          L-ISOLEUCINE
```


The output of spectralMatching returns a list containing the following elements:

**q_dbPth**: Path of the query database (this will have been updated with the annotation results if updateDb argument used)

**xcmsMatchedResults**
  If the qeury spectra had XCMS based chromotographic peaks tables (e.g c_peak_groups, c_peaks) in the sqlite database - it will
 be possible to summarise the matches for each XCMS grouped feature. The dataframe contains the following columns

* lpid - id in database of library spectra
* qpid - id in database of query spectra
* dpc - dot product cosine of the match
* rdpc - reverse dot product cosine of the match
* cdpc - composite dot product cosine of the match
* mcount - number of matching peaks
* allcount - total number of peaks across both query and library spectra
* mpercent - percentage of matching peaks across both query and library spectra
* library_rt - retention time of library spectra
* query_rt - retention time of query spectra
* rtdiff - difference between library and query retention time
* library_precursor_mz - library precursor mz
* query_precursor_mz - query precursor mz
* library_precursor_ion_purity - library precursor ion purity
* query_precursor_ion_purity - query precursor ion purity
* library_accession -  library accession value (unique string or number given to eith MoNA or Massbank data entires)
* library_precursor_type - library precursor type (i.e. adduct)
* library_entry_name - Name given to the library spectra
* inchikey - inchikey of the matched library spectra
* library_source_name - source of the spectra (e.g. massbank, gnps)
* library_compound_name - name of compound spectra was obtained from

**matchedResults**

All matched results from the query spectra to the library spectra. Contains the same columns as above
but without the XCMS details. This table is useful to observe spectral matching results
for all MS/MS spectra irrespective of if they are linked to XCMS MS1 features.


It should be noted that in a typical Data Dependent Acquisition (DDA) experiment not all the fragmentation scans collected can be linked backed to an associated XCMS features and in some cases the percentage of XCMS features with fragmentation spectra can sometimes be quite small.


# References