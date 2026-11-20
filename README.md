# msPurity

**An R/Bioconductor package for LC-MSn-based metabolomics, providing tools for
precursor ion purity assessment, processing fragmentation data and spectral
trees, and spectral matching.**

[![Git](https://img.shields.io/badge/repository-GitHub-blue.svg?style=flat&maxAge=3600)](https://github.com/computational-metabolomics/msPurity)
[![Checks](https://github.com/computational-metabolomics/msPurity/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/computational-metabolomics/msPurity/actions/workflows/check-bioc.yml)
[![Bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat&maxAge=3600)](https://bioconda.github.io/recipes/bioconductor-mspurity/README.html)
[![License](https://img.shields.io/badge/licence-GNU_v3-teal.svg?style=flat&maxAge=3600)](https://www.gnu.org/licenses/gpl-3.0.html)
[![DOI](https://img.shields.io/badge/DOI-10.18129/B9.bioc.msPurity-teal.svg?style=flat&maxAge=3600)](https://doi.org/10.18129/B9.bioc.msPurity)
[![Paper](https://img.shields.io/badge/paper-Analytical_Chemistry-teal.svg?style=flat&maxAge=3600)](https://doi.org/10.1021/acs.analchem.6b04358)

**Bioconductor release:** [![Bioconductor release
availability](https://bioconductor.org/shields/availability/release/msPurity.svg)](https://bioconductor.org/packages/release/bioc/html/msPurity.html#archives)
[![Bioconductor release
downloads](https://bioconductor.org/shields/downloads/release/msPurity.svg)](https://bioconductor.org/packages/stats/bioc/msPurity/)
[![Bioconductor release build
results](https://bioconductor.org/shields/build/release/bioc/msPurity.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/msPurity/)

**Bioconductor devel:** [![Bioconductor devel
availability](https://bioconductor.org/shields/availability/devel/msPurity.svg)](https://bioconductor.org/packages/devel/bioc/html/msPurity.html#archives)
[![Bioconductor devel
downloads](https://bioconductor.org/shields/downloads/devel/msPurity.svg)](https://bioconductor.org/packages/stats/bioc/msPurity/)
[![Bioconductor devel build
results](https://bioconductor.org/shields/build/devel/bioc/msPurity.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/msPurity/)

See the [documentation](https://bioconductor.org/packages/msPurity/) for worked
examples and the reference manual, and [NEWS](NEWS.md) for changes.

## Version 2

Version 2 extends the established precursor ion purity assessment and LC-MS/MS
processing workflow (version 1) with the following additional and improved functionalities:

- **Standardised data handling.** Use Spectra-based access to measured and
  processed spectra, with supported inputs from Spectra, MsExperiment and
  XcmsExperiment. MsBackendParquet connects Spectra to spectra stored in Parquet
  datasets.
- **MSn data processing.** Associate MSn spectra with chromatographic
  features through their acquisition relationships. Filter product peaks and
  average spectra using explicit acquisition conditions, retaining separate
  groups where precursor paths or conditions are incompatible. Assess precursor
  ion purity using neighbouring MS1 spectra for MS2 and the linked parent
  fragmentation spectrum for higher MSn levels.
- **Reproducible settings and results.** Reuse parameter objects for purity
  assessment, feature association, filtering and averaging. Processing history
  records the settings applied and the measurements used. Reopen selected
  processed spectra and their provenance without repeating the full workflow.
- **Spectral trees.** Explore Mass Spectral Trees linking successive spectra and
  Precursor-Fragment Trees linking individual precursor and product peaks.
  Follow parents, children and fragmentation paths, inspect supporting
  measurements, and exchange tree representations for further visualisation.
- **MSn spectral matching.** Compare measured or consensus spectra with other
  experimental spectra and reference libraries, accounting for MS level,
  precursor relationships and acquisition conditions. For approximate cosine
  similarity, `BlinkMatchingParam()` selects the BLINK method of
  [Harwood et al. (2023)](https://doi.org/10.1038/s41598-023-40496-9).
- **Spectrum selection and visualisation.** Select spectra and filter their
  peaks through Spectra. Linked views using lcmsPlot connect chromatographic
  peaks, fragmentation events, spectra and matching results.
- **Library exchange and reporting.** Import reference libraries through
  supported formats and Spectra backends, and convert them to Parquet for
  matching. Exchange spectra in formats such as MSP, and export feature
  information, spectral matches and spectrum references for downstream analysis.

The DIMS workflow, including `purityD()` and `dimsPredictPuritySingle()`, is
deprecated in version 2 and will receive no further bug fixes or feature
development. There is no direct replacement in msPurity. Existing calls remain
functional during deprecation; see
`help("msPurity-deprecated", package = "msPurity")` for the retirement path.

### Parquet connects raw data, processed results and libraries

Version 2 uses Parquet datasets to retain spectra alongside the information
needed to interpret and reuse them. The workflow connects three kinds of data:

- **Raw measurements:** convert MSn mzML data with
  `MsBackendParquet::mzMLToParquet()` and read the resulting spectra through
  Spectra. Spectrum identities and parent references link fragmentation scans
  back to their acquisition context.
- **Processed results:** export with `createDatabase(..., format = "parquet")`
  to store feature associations, filtering records, averaged spectra, processing
  settings and measurement provenance. Results reference the raw study spectra;
  keep the linked study dataset available when reopening raw measurements.
  `readProcessedSpectra()` and `readProcessedProvenance()` reopen selected
  results and their processing records without repeating the analysis.
- **Reference libraries:** use `convertLibraryToParquet()` to convert supported
  library inputs, including Spectra objects, into Parquet datasets for spectral
  matching.

Parquet is the result-storage route for MSn data processing and acquisition-aware
averages. The older SQLite route is limited to MS2 results; supported existing
SQLite databases can be converted with `convertSqliteToParquet()`.

## Version 1

Version 1 established the core workflow for assessing precursor ion purity and
processing fragmentation spectra:

- Assess precursor ion purity for acquired MS/MS spectra.
- Predict precursor ion purity for xcms LC-MS features before fragmentation data
  are acquired.
- Associate fragmentation spectra with xcms features.
- Filter and average MS/MS spectra within and across samples.
- Create spectral databases with feature information and annotations.
- Match query MS/MS spectra against library spectra.
- Export fragmentation spectra in MSP format.

These capabilities provide the foundation for version 2.

## What is precursor ion purity?

Precursor ion purity measures the contribution of a selected precursor peak to
the signal within the isolation window used for fragmentation. In its simplest
form, it is the intensity of the selected precursor divided by the total
intensity in that window.

For an MS/MS acquisition, msPurity uses the surrounding MS1 measurements to
estimate purity at the time of fragmentation. Options account for isotopic
peaks, low-intensity signals and the isolation efficiency of the instrument. The
available isolation metadata and surrounding measurements determine whether an
estimate can be made.

A high value indicates that the selected precursor dominates the measured
isolation-window signal. It does not guarantee an uncontaminated fragmentation
spectrum or establish a compound identification.

## Installation

### Bioconductor release

For routine analysis, install the stable release through Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("msPurity")
```

The version installed depends on your R and Bioconductor release. Check it with:

``` r
packageVersion("msPurity")
```

### Development version

The [Bioconductor development
version](https://bioconductor.org/packages/devel/bioc/html/msPurity.html) and
the [GitHub master
branch](https://github.com/computational-metabolomics/msPurity/tree/master)
provide access to ongoing development. GitHub and Bioconductor devel can differ;
use an R/Bioconductor environment compatible with the dependencies declared by
the source version you install.

To install from GitHub:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("computational-metabolomics/msPurity")
```

## Documentation and related tools

- [Bioconductor release
  documentation](https://bioconductor.org/packages/release/bioc/html/msPurity.html):
  installation, vignettes and reference manual for the stable release.
- [Bioconductor development
  documentation](https://bioconductor.org/packages/devel/bioc/html/msPurity.html):
  vignettes and manual for the Bioconductor development version.
- Installed documentation: `browseVignettes("msPurity")` and
  `help(package = "msPurity")`.
- [Release notes](NEWS.md): user-visible changes and compatibility notes.
- [Galaxy tools](https://github.com/computational-metabolomics/mspurity-galaxy):
  associated workflow tools; their supported functionality may differ from the R
  interfaces.
- [Bioconda](https://anaconda.org/bioconda/bioconductor-mspurity): an
  alternative distribution; check the available version before installation.

## Citation

If you use msPurity in your research, please cite:

Lawson, T.N., Weber, R.J., Jones, M.R., Chetwynd, A.J., Rodriguez Blanco, G.A.,
Di Guida, R., Viant, M.R. and Dunn, W.B. (2017). *msPurity: Automated Evaluation
of Precursor Ion Purity for Mass Spectrometry Based Fragmentation in
Metabolomics.* Analytical Chemistry, 89, 2432-2439.
[doi:10.1021/acs.analchem.6b04358](https://doi.org/10.1021/acs.analchem.6b04358).

```r
citation("msPurity")
```

Subsequent development of the version reported in the publication above,
including fragmentation spectrum processing, spectral database creation, spectral
matching and annotation handling, is described alongside the application of
Galaxy workflows in the following paper:

Lawson, T.N., Jones, M.R., Chetwynd, A.J., Sostare, E., Weidt, S., Mistrik, R.,
Dunn, W.B., Weber, R.J.M. and Viant, M.R. (2026). *Experimental and computational
approaches for deep metabolome annotation with application to the ecotoxicological
model organism Daphnia magna.* GigaScience, 15, giag055.
[doi:10.1093/gigascience/giag055](https://doi.org/10.1093/gigascience/giag055).

If you use BLINK for spectral matching, please also cite its original publication.

Harwood, T.V., Treen, D.G.C., Wang, M., de Jong, W., Northen, T.R. and
Bowen, B.P. (2023). *BLINK enables ultrafast tandem mass spectrometry cosine
similarity scoring.* Scientific Reports, 13, 13462.
[doi:10.1038/s41598-023-40496-9](https://doi.org/10.1038/s41598-023-40496-9).

## Authors and contributions

See [DESCRIPTION](DESCRIPTION) for the full list of authors and contributors and
their specific roles.

### Version 1

**Main contributor, author and previous maintainer:** Thomas N. Lawson.

Contributors: Ralf Weber, Martin Jones, Julien Saint-Vanne, Andris Jankevics and
Ossama Edbali.

Mark Viant and Warwick Dunn supervised Thomas N. Lawson's PhD research at the
University of Birmingham, of which this work formed part.

### Version 2

**Main contributors and authors:**

- Ossama Edbali (<o.edbali@bham.ac.uk>) — University of Birmingham (UK).
- Ralf Weber (<r.j.weber@bham.ac.uk>) — University of Birmingham (UK).

**Maintainer:** Ossama Edbali.

Ralf Weber supervised Ossama Edbali's PhD research at the University of
Birmingham, of which this work formed part.

### Community acknowledgements

msPurity builds on software developed by the
[R for Mass Spectrometry](https://www.rformassspectrometry.org/) and
[Bioconductor](https://bioconductor.org/) communities. We thank the developers
and contributors whose shared software, documentation and standards make this
work possible.

### Contributing

Bug reports, questions, suggestions and contributions are welcome through
[GitHub issues](https://github.com/computational-metabolomics/msPurity/issues).
For a bug report, include a small reproducible example, the error message and
`sessionInfo()`. For a substantial change, open an issue to discuss the scope
before submitting a pull request.
