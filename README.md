[![Build Status](https://travis-ci.org/Viant-Metabolomics/msPurity.png)](https://travis-ci.org/Viant-Metabolomics/) [![Build Status windows](https://ci.appveyor.com/api/projects/status/github/viant-metabolomics/mspurity?branch=master&svg=true)](https://ci.appveyor.com/project/tomnl/mspurity/)

Build tests for travis and AppVeyor 


# msPurity

## About
The importance of assessing the contribution of the precursor peak within an isolation window targeted for fragmentation has been previously detailed in proteomics but to date there has been little attention of this data-processing technique in metabolomics. Here we present msPurity, a vendor independent R package for Liquid Chromatography-Mass spectrometry (LC-MS) and Direct Infusion-Mass Spectrometry (DI-MS) that calculates a simple metric to describe the contribution of the selected precursor peak.What we call "Precursor purity" is a measure of the contribution of a selected precursor peak in an isolation window used for fragmentation. The simple calculation involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. When assessing MS/MS spectra this calculation is done before and after the MS/MS scan of interest and the purity is interpolated at the time of the MS/MS acquisition. The calculation is very similar to the "Precursor Ion Fraction" (PIF) metric described by  \[1\] for proteomics with the exception that purity here is interpolated at the recorded point of MS2 acquisition using bordering full-scan spectra. Additionally, low abundance ions that are remove that are thought to have limited contribution to the resulting MS2 spectra and can optionally take into account the isolation efficiency of the mass spectrometer

Use the following link for more details:

http://bioconductor.org/packages/msPurity/


## Install

### Bioconductor

```
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
useDevel() # currently only available in development branch
biocLite("msPurity")

```

### Github

```

library(devtools)
library(testthat)

t <- try(install_github('Viant-Metabolomics/msPurity'))

if("try-error" %in% class(t)){
  print("INSTALLATION FAILED!")
} else {
  print("INSTALLATION SUCCESS!")
  print("PERFORM UNIT TESTING")
  # Perform unit test to make sure the package is working as expected
  tp <- test_package('msPurity')
  print(tp)
}
```



## Ref
[1] Michalski, A., Cox, J., & Mann, M. (2011). More than 100,000 detectable peptide species elute in single shotgun proteomics runs but the majority is inaccessible to data-dependent LC-MS/MS. Journal of Proteome Research, 10(4), pp. 1785-1793.



