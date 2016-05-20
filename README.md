# msPurity

## About
The importance of assessing the contribution of the precursor peak within an isolation window targeted for fragmentation has been previously detailed in proteomics but to date there has been little attention of this data-processing technique in metabolomics. Here we present msPurity, a vendor independent R package for Liquid Chromatography-Mass spectrometry (LC-MS) and Direct Infusion-Mass Spectrometry (DI-MS) that calculates a simple metric to describe the contribution of the selected precursor peak. What we call here “precursor purity” is calculated as per the Michalski approach \[1\] (intensity of selected precursor divided by total intensity of the isolation window) with the exception that the metric is interpolated at the time of the MS/MS scan. 

The R package can also be used to predict the precursor purity of subsequent MS/MS run from a prior MS full scan dataset.

## Workflows

![alt tag]()


## How to install the package and run the unit tests

```

library(devtools)
library(testthat)

token <- '74d3aa11143449965fd7cc66b57585bc6fc03758'
t <- try(install_github('Viant-Metabolomics/msPurity', auth_token=token))

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


