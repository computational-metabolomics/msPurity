============
msPurity: Package to assess precursor ion purity and perform spectral matching
============




**General:**

|Git| |Bioconda| |Build Status (Travis)|  |License| |DOI| |Paper|


**Bioconductor (release):**

|Bioconductor release availability| |Bioconductor release downloads|  |Bioconductor release build results|

**Bioconductor (devel):**

|Bioconductor devel availability| |Bioconductor devel downloads|  |Bioconductor release build results|

------------
Which version to use?
------------

Recommendation for most uses cases is to install and use the `Bioconductor stable version <http://bioconductor.org/packages/msPurity/>`_ of msPurity.

The code available from both the `Bioconductor development branch <http://bioconductor.org/packages/devel/bioc/html/msPurity.html>`_ and `the master branch on github <https://github.com/computational-metabolomics/mspurity>`_ has the newest functionality.

------------
About
------------

The importance of assessing the contribution of the precursor peak within an isolation window targeted for fragmentation has been previously detailed in proteomics but to date there has been little attention of this data-processing technique in metabolomics. Here we present msPurity, a vendor independent R package for Liquid Chromatography-Mass spectrometry (LC-MS) and Direct Infusion-Mass Spectrometry (DI-MS) that calculates a simple metric to describe the contribution of the selected precursor peak. What we call "Precursor ion purity" is a measure of the contribution of a selected precursor peak in an isolation window used for fragmentation. The simple calculation involves dividing the intensity of the selected precursor peak by the total intensity of the isolation window. When assessing MS/MS spectra this calculation is done before and after the MS/MS scan of interest and the purity is interpolated at the recorded time of the MS/MS acquisition. The calculation is similar to the "Precursor Ion Fraction" (PIF) metric described by  \[1\] for proteomics with the exception that metric here is interpolated at the recorded point of MS2 acquisition using bordering full-scan spectra, isotopic peaks can be removed, low abundance peaks are removed that are thought to have limited contribution to the resulting MS2 spectra and the isolation efficiency of the mass spectrometer can be used to normalise the intensities used for the calculation.

The software can either be used on fragmentation spectra that has been collected to assess the precursor ion purity or it can be used to calculate the anticipated precursor ion purity of a future run based on MS1 spectra.

Since the paper has been published the package now supports spectral matching, SQLite interaction and some of the functionality is now implemented in Galaxy.

Associated paper  `msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics. Analytical Chemistry <http://pubs.acs.org/doi/abs/10.1021/acs.analchem.6b04358>`_ [2]

Use the following links for more details:

* Bioconductor: http://bioconductor.org/packages/msPurity/
* Vignette: https://bioconductor.org/packages/devel/bioc/vignettes/msPurity/inst/doc/msPurity-vignette.html
* Manual: http://bioconductor.org/packages/devel/bioc/manuals/msPurity/man/msPurity.pdf
* Galaxy implementation: https://github.com/computational-metabolomics/mspurity-galaxy
* Bioconda (stable): https://anaconda.org/bioconda/bioconductor-mspurity
* Conda (dev and testing): https://anaconda.org/tomnl/bioconductor-mspurity



------------
Install
------------

Bioconductor
------------

.. code-block:: r

  if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
  BiocManager::install("msPurity", version = "3.8")



Github
------------

.. code-block:: r

  library(devtools)
  library(testthat)
  install_github('computational-metabolomics/msPurity')
 


------------
Ref
------------
[1] Michalski, A., Cox, J., & Mann, M. (2011). More than 100,000 detectable peptide species elute in single shotgun proteomics runs but the majority is inaccessible to data-dependent LC-MS/MS. Journal of Proteome Research, 10(4), pp. 1785-1793.

[2] Lawson, T.N., Weber, R.J., Jones, M.R., Chetwynd, A.J., Rodriguez Blanco, G.A., Di Guida, R., Viant, M.R. and Dunn, W.B., 2017. msPurity: Automated Evaluation of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in Metabolomics. Analytical Chemistry.


.. |Bioconductor release availability| image:: https://bioconductor.org/shields/availability/3.8/msPurity.svg
   :target: https://bioconductor.org/packages/release/bioc/html/msPurity.html#archives


.. |Bioconductor devel availability| image:: https://bioconductor.org/shields/availability/3.9/msPurity.svg
   :target: https://bioconductor.org/packages/devel/bioc/html/msPurity.html#archives

.. |Bioconductor release downloads| image:: https://bioconductor.org/shields/downloads/release/msPurity.svg
   :target: http://bioconductor.org/packages/stats/bioc/msPurity/

.. |Bioconductor devel downloads| image:: https://bioconductor.org/shields/downloads/devel/msPurity.svg
   :target: http://bioconductor.org/packages/stats/bioc/msPurity/


.. |Bioconductor release build results| image:: https://bioconductor.org/shields/build/release/bioc/msPurity.svg
   :target: http://bioconductor.org/checkResults/release/bioc-LATEST/msPurity/

.. |Bioconductor devel build results| image:: https://bioconductor.org/shields/build/devel/bioc/msPurity.svg
   :target: http://bioconductor.org/checkResults/devel/bioc-LATEST/msPurity/



.. |Build Status (Travis)| image:: https://img.shields.io/travis/computational-metabolomics/msPurity/master.svg?label=Travis
   :target: https://travis-ci.org/computational-metabolomics/msPurity

.. |Build Status (AppVeyor)| image:: https://ci.appveyor.com/api/projects/status/github/computational-metabolomics/mspurity?branch=master&svg=true
   :target: https://ci.appveyor.com/project/Tomnl/msPurity

.. |Git| image:: https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600
   :target: https://github.com/computational-metabolomics/msPurity

.. |Bioconda| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600
   :target: https://bioconda.github.io/recipes/bioconductor-mspurity/README.html

.. |License| image:: https://img.shields.io/badge/licence-GNU_v3-teal.svg?style=flat&maxAge=3600
   :target: https://www.gnu.org/licenses/gpl-3.0.html

.. |DOI| image:: https://img.shields.io/badge/DOI-10.18129/B9.bioc.msPurity-teal.svg?style=flat&maxAge=3600
   :target: https://doi.org/doi:10.18129/B9.bioc.msPurity

.. |Paper| image:: https://img.shields.io/badge/paper-Analytical_Chemistry-teal.svg?style=flat&maxAge=3600
   :target: http://doi.org/10.1021/acs.analchem.6b04358
