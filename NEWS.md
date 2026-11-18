# msPurity 1.99.1

- **Raw MS1-MS4 data in Parquet studies** (#127): a `Spectra` on a dataset
  written by MsBackendParquet's `mzMLToParquet()` (0.99.12 or later) runs
  through `purityA()`, `frag4feature()`, `filterFragSpectra()` and averaging
  with the same results as the mzML files, for MS2 to MS4. Peaks are read
  from the dataset, never from the mzML files, so these need not exist; a
  dataset that has moved or changed since `purityA()` read it is reported as
  such. A selected or processed input numbers its raw peak rows as the
  matching mzML input does. An activation method the dataset records as a
  PSI-MS accession is read as its method name, so selectors match it as they
  match a filter string, and `activationTypeBasis` records which applied.
- With automatic `parentReference`, a Parquet dataset's own spectrum ids are
  no longer read as parent scan numbers once acquisition numbers resolve,
  which made parent references ambiguous in studies of several runs.
- **Parquet results written directly, for MS2 to MS4** (#127, #128):
  `createDatabase(format = "parquet")` no longer stages an SQLite database;
  the MS2 output is unchanged. Results version 1.1.0 adds, without changing
  what 1.0.0 datasets hold:
  - MS3 and MS4 scans in `x_mspurity_scan`, with `ms_level`, MSn stage
    purity and status, the activation method used and its basis, and each
    scan's parent and MS2 root, as references and as rows of the table;
  - MS3/MS4 feature links in the new `x_mspurity_msn_feature`, and the
    retention-time mapping of MS2 links in `spectrum_feature`;
  - filter flags of every filtered spectrum in `x_mspurity_scan_peak`,
    including those of feature-linked scans when `allfrag = FALSE`, which
    were lost before;
  - each step's processing history in `tool_provenance`.
- **Study links without the mzML files** (#127): with `study`, scans and the
  MS1 scans their MS2 precursors name are matched from the study's own
  metadata by native id, then scan number at the same position and level;
  files are no longer read or checksummed for this. When `purityA()` read
  a Parquet dataset it is the default `study`; `study = NA` opts out.
  Conditions the study does not record are reported as unavailable in the
  loss ledger.
- The SQLite writer refuses MS3/MS4 feature links as it refuses MS3/MS4
  spectra.

- **Retention-time correction in `frag4feature()`** (#137): precursor times
  of each file are mapped to the adjusted XCMS scale by interpolation and
  compared with peak and feature windows on that scale, in both modes.
  Before, individual mode converted peak bounds back to raw time by exact
  match, which silently dropped every filled peak, and grouped mode compared
  raw precursor times with adjusted feature windows, so `convert2RawRT` had
  no effect there.
- purityA files are matched to XCMS samples by file name instead of by
  position. Each file gets a mapping status (`aligned`, `unaligned`,
  `external`, `estimated` or `missing`) in the new `rtMapStatus` column, and
  files without a mapping are named in a warning.
- Filled peaks are excluded explicitly; `includeFilled = TRUE` links them,
  flagged by `is_filled`. `intense` applies in both modes, falling back from
  `iMz` to `aMz` to `precursorMZ`, and `precurMzSource` records the value
  used.
- New `frag4feature()` arguments: `ms2RtMap` maps DDA files outside the XCMS
  object (an XCMS object, named functions or named offsets),
  `estimateRtMap = TRUE` estimates the shift for them in grouped mode,
  `rtTol` widens the windows, `rtSource = "ms2"` uses the MS/MS scan time
  and `groupWindow = "median"` uses median feature bounds.
  `Frag4FeatureParam` holds the same settings, apart from `ms2RtMap`.
- New `grped_df` columns `matchRank` and `rtApexDist` order the matches of
  each scan; per-file diagnostics are stored in `pa@params$frag4feature`.
  Legacy columns keep their order and types. MS3 and MS4 rows keep their own
  precursor m/z and time.
- **Faster linking**: the nested `plyr::ddply` matching and the SNOW cluster
  are replaced by one `dplyr` interval join, so `dplyr (>= 1.1.0)` is
  required. On 32 files linking takes about 5 s instead of 870 s. Processing
  history records one effective worker, and grouped mode now records the RT
  conversion as used for association.
- Data frames passed to the deprecated `grp_peaklist`/`grpPeaklist`
  arguments no longer error, the `plim` filter handles missing purity, and
  the SQLite and Parquet writers and `createMSP()` prefer `fileid` over
  `sample` when mapping links to files.

# msPurity 1.99.0

- **MS3 and MS4 spectra** (#138): `purityA()`, `frag4feature()` and
  `filterFragSpectra()` take `msLevel` to process MS2 to MS4. `purityA()`
  calculates MS3/MS4 precursor purity from the linked parent spectrum;
  `frag4feature(msLevel = 2:4)` associates MS3/MS4 spectra with features
  through their validated MS2 ancestors, with `parentReference` choosing how
  parent spectra are identified. `filterFragSpectra()` applies purity
  thresholds to each spectrum's own assessment, so a failed parent does not
  exclude its descendants. MS3/MS4 averaging and processed export are not
  supported yet.
- **Parameter objects**: `PurityAParam`, `Frag4FeatureParam` and
  `FilterFragSpectraParam` extend `ProtGenerics::Param` and are passed
  through `param`. Purity and filtering settings can be assigned by MS level,
  activation method, collision energy or mass analyzer. Explicit arguments
  that conflict with `param` are an error.
- **Processing history** records the resolved settings and their
  assignments, the MS1 measurements used for each purity estimate, the
  inputs used for feature association, the RT conversion and the number of
  workers, and which filtering call produced each purity flag on extracted
  Spectra.
- **Bug fixes** (#133, #103, #107): each MS2 event uses its recorded
  isolation window, with explicit offsets taking precedence; an Agilent
  isolation width is no longer assumed. Isolation-efficiency weighting is
  corrected for asymmetric windows. Purity is `NA` with a distinct reason
  outside valid isolation bounds or when no neighbouring MS1 measurement can
  be used, and a later spectrum no longer replaces a missing preceding one.
- `purityTable()` returns MS2 rows only. Filtering reads unchanged peak
  arrays once, and equivalent parameter settings are resolved once.

# msPurity 1.39.4

- **Faster on Windows** (#136): raw spectra are read serially through
  Spectra. The default `bpparam()` on Windows is a SnowParam that started a
  new cluster for every read, which made `purityA()`, `frag4feature()` and
  the fragment filters many times slower there. msPurity still runs files in
  parallel through its `cores` argument. BiocParallel moves from Suggests to
  Imports.
- Examples load stored objects instead of rerunning the xcms and purityA
  pipelines, `spectralMatching()` examples use the bundled mini library
  instead of downloading the full one, and `createDatabase()`,
  `create_database()` and `createMSP()` examples write to `tempdir()`.

# msPurity 1.39.3

- **Parquet built on MsBackendParquet** (#126): manifests, mzStack
  conditions, and writing and reading spectra now use the MsBackendParquet
  dataset API, so msPurity and MsBackendParquet apply one set of mzStack
  rules. MsBackendParquet (>= 0.99.10), arrow and jsonlite move from Suggests
  to Imports.
- Parquet errors carry both `mzstack_<class>` and `msPurity_parquet_<class>`
  classes; `parquetCapabilities()` reports the `mzstack_*` names. A manifest
  commit refused by another writer keeps `msPurity_parquet_conflict`.
- Parquet datasets store run paths relative to the dataset, so a dataset can
  be moved or renamed. `spectrum_index` is int64 in every run, part files are
  named by MsBackendParquet, core variables with no values in a run are not
  written, and library sources without spectra are not written as empty
  runs. The Parquet format is unreleased, so no migration is needed.

# msPurity 1.39.2

- **Bug fix**: identifiers for processing records and Parquet datasets could
  repeat when created in quick succession, most often on macOS, where the
  clock that seeded them has microsecond resolution. A repeated processing
  identifier replaced an earlier record and made the processing history
  invalid. Identifiers now come from a random stream seeded once per R
  process, and the caller's random number stream is still left unchanged.

# msPurity 1.39.1

- **Parquet output**: new Parquet output,
  alongside SQLite. `createDatabase()`, `spectralMatching()` and
  `combineAnnotations()` take `format = "parquet"`; SQLite remains the
  default. Needs the suggested packages arrow and jsonlite.
- **Parquet conversion and validation**: new functions
  `convertSqliteToParquet()` converts an msPurity SQLite database to a Parquet
  dataset, `convertLibraryToParquet()` converts an msp2db or MSP spectral
  library or a Spectra object (such as an MSP file read with MsBackendMsp),
  `validateParquet()` checks a Parquet dataset against the specification, and
  `parquetCapabilities()` reports what the Parquet output supports.
- **Faster spectral matching on Parquet**: `spectralMatching()` on Parquet
  datasets gives the same scores as on SQLite and runs much faster against
  large libraries. New argument `topn` keeps the best-ranked matches per
  query. Their peaks are read through MsBackendParquet when it is installed.
- **SQLite output deprecated**: every function that writes SQLite raises one
  deprecation warning per call. Its output is unchanged.
- **Spectra containers**: msPurity reads raw data through Spectra and stores
  spectra in Spectra objects. All results are unchanged. Spectra and
  MsExperiment move from Suggests to Imports.
- **New slots**: purityA objects gain `spectra`, `fragSpectra` and `avSpectra`
  slots: the MS/MS scans with their purity as spectra variables, the scans
  linked to features (each once) with filter flags as peak variables, and the
  averaged spectra. They also gain a `params` slot recording the arguments of
  `purityA()` and `frag4feature()`. purityD objects gain `avSpectra` (the
  averaged peak lists) and `experiment` (the files as an MsExperiment) slots.
- **New accessors**: `purityTable()`, `groupedSpectra()`, `allFragSpectra()`,
  `averagedSpectra()` and `averagedPeaks()` return either the legacy data
  frames and lists or Spectra objects. Use them instead of reading slots
  directly.
- **Legacy slots**: the legacy slots (`puritydf`, `grped_ms2`,
  `all_frag_scans`, `av_spectra` and `avPeaks`) are still filled by default.
  `options(msPurity.legacySlots = FALSE)` leaves them empty; outputs are the
  same either way. The default is planned to change in the next release.
- **Saved objects**: `updateObject()` converts purityA and purityD objects
  saved by earlier versions; every method does this on its input.
- **More input types**: `purityA()` accepts Spectra, MsExperiment and
  XcmsExperiment objects as well as file paths, and reads each file once.
  `purityD()` accepts an MsExperiment. `frag4feature()` and
  `createDatabase()` accept XcmsExperiment objects.
- The `mzRback` argument is deprecated and ignored.
- **Processing history**: new functions `processingHistory()` and
  `processingParameters()` for purityA objects and extracted Spectra.
  `purityA()`, `frag4feature()`, `filterFragSpectra()` and the averaging
  functions each record their effective settings (defaults included), the
  software version, execution times and references to the spectra they read
  and wrote, not copies of peaks. Repeated calls keep earlier records. A
  subset of extracted Spectra returns only the history of its own spectra.
  Objects saved by earlier versions get their known settings imported, marked
  incomplete, with unknown dates left unknown.
- **Resolved purity settings**: `purityA()` records, for each file, the
  isolation offsets it used and where they came from (the caller, the
  recorded metadata, or the Agilent approximation), and the `iwNormFun` it
  applied when `iwNorm = TRUE`.
- **Filtering by acquisition condition**: new `FilterFragSpectraParam()`, a
  `ProtGenerics::Param` holding complete filter settings, and a `param`
  argument to `filterFragSpectra()` taking one such object or a named list of
  them. Its `select` field assigns settings to spectra by MS level, mass
  analyzer, activation type, collision energy and its unit, and polarity,
  without discarding spectra. Spectra matched by no set use the call's own
  arguments, and overlapping sets with different settings stop with an
  error. CID, HCD and ETD are recognised by abbreviation, CV name or
  accession. Where structured metadata are missing, ordinary scan filter
  strings supply CID/HCD and FTMS/ITMS. `plim = NA` skips purity filtering
  for the spectra a set matches. Filtering with parameter sets needs
  `cores = 1`. ProtGenerics is added to Imports.
- **Peak provenance**: filtered fragmentation spectra keep each peak's
  position in its source spectrum (`original_peak_index`), and source peak
  ids (`sid`) count the rows of the source data, including those of
  intervening MSn spectra. Averaged spectra record the feature links they
  were built from (`msPurityContributors`, `msPurityClusterContributors`) and
  whether their inputs share one acquisition configuration
  (`msPurityAcquisitionConfiguration`: uniform, incomplete, mixed or
  unknown). Averaging itself is unchanged, so MS2 averages can still combine
  acquisition conditions.
- **Stale results cleared**: running `frag4feature()` or
  `filterFragSpectra()` again clears the filter flags and averages that
  depended on the previous run, and recomputing intra-file averages clears
  the inter-file averages. Their history records are kept.
- `allFragSpectra()` names the scan of every peak correctly. The legacy
  `all_frag_scans` table, and the scan peaks in SQLite output, number scans
  with a count that skips scans without peaks, so after such a scan they name
  the wrong scan; they are kept as they were because the SQLite output is
  frozen.
- **Bug fix**: raw peaks and their filter flags could be assigned to the
  wrong spectrum when a file had empty MS2 spectra or MS3/MS4 spectra between
  them. Peaks are now matched by file and spectrum position. When converting
  SQLite to Parquet, raw peaks sharing a `pid` with an averaged spectrum were
  added to that average; an average now contains only its own peaks.
  Ambiguous assignments are refused and inconsistent identities in older
  databases raise a warning (#124).
- **Bug fix**: duplicate `seqNum` columns stopped new and RDS-reopened
  objects from being exported to SQLite or Parquet (#124).
- **Bug fix**: spectral matching could pair a fragment peak outside the m/z
  tolerance when it tied with an eligible peak on relative intensity
  difference or mass error. Only eligible, unused library peaks are now
  considered, in both the current and the deprecated aligners, for SQLite
  and Parquet (#125).
- **Bug fix**: averaging after `filterFragSpectra(rmp = TRUE)` failed when
  every spectrum of a feature lost all its peaks. Such features now get no
  averaged spectrum, and the other averages equal those of `rmp = FALSE`.
- **Bug fix**: `createMSP(method = "max")` failed when one peak of the chosen
  spectrum passed the filters, and could write the wrong scan when a scan was
  linked to a feature more than once. It now writes the scan with the highest
  precursor intensity and skips spectra with no passing peaks.
- **Bug fix**: `createDatabase()` called `featureValues()` without the
  `xcms::` prefix and failed unless xcms was attached.
- **Test data**: four reduced mzML files with MS1-only and MS1 to MS4
  acquisitions of L-tryptophan and L-kynurenine (#123). BiocParallel is added
  to Suggests.

# msPurity 1.37.3

- Update within the tests only (not within the R code base): library database path for lvl spectral matching testing (required after removal of the local SQLite database from msPurityData).

# msPurity 1.37.2

- Bug fix for subsetting the grouped data frame in the `average_xcms_grouped_msms_indiv` function

# msPurity 1.37.1

- Removed library_spectra.db direct reference from msPurityData package (too large for Bioconductor)
- spectralMatching() now downloads default library database from Zenodo (https://zenodo.org/records/18700802) and uses BiocFileCache with MD5 verification
- spectral_matching() (deprecated function) now requires explicit library_db_pth parameter; default library no longer available
- Added BiocFileCache to Imports for automatic caching of remote databases

# msPurity 1.30.1 / 1.31.1

- Update example for purityX causing build error

# msPurity 1.27.1

- createMSP fix - now uses the median precursor MZ and precursor RT in the MSP file

# msPurity 1.21.2

- Fixes due to mz and intensity being named in columns for mzR and XCMS
- Fixes for connection{base} file() opening changes (no longer accept w+a in function in R v4.2)
- Update tests for the above
- Update a reference of grpid in averageXFragSpectra functions more explicit (no change in functionality)
- Typo fix in vignette

# msPurity 1.21.1

- Bugfix for frag4feature for XCMS 3 compatability https://github.com/computational-metabolomics/msPurity/pull/93
- Remove imports that are no longer used

# msPurity 1.19.2

- XCMS 3 compatability update (M-R-JONES) https://github.com/computational-metabolomics/msPurity/pull/91

# msPurity 1.19.1

- Bug fix for flagRemove (full width was not calculated as expected)

# msPurity 1.17.2

- Update dev to match bug fixes in master

# msPurity 1.16.2

- Author list updated

# msPurity 1.17.1

- Update dev to match bug fixes in master

# msPurity 1.16.1

- Update of createDatabase to record all intra average spectra in database
- Add license and copyright info to code
- Add github workflow CI (and subsequent formatting updates to pass tests)

# msPurity 1.15.1

- Update dev to match bug fixes in master

# msPurity 1.14.1

- Update of rdpc algorithm (see https://github.com/computational-metabolomics/msPurity/issues/78)
- Update of align algorithm (see https://github.com/computational-metabolomics/msPurity/issues/79)
- Fix for spectralMatching of type 'scan' previously incorrectly outputing no matches

# msPurity 1.13.2

- Update dev to match bug fixes in master

# msPurity 1.12.2

- Fix to allow large files to be processed with combineAnnotation (in some cases the summary report should be excluded)

# msPurity 1.13.1

- Update dev to match bug fixes in master

# msPurity 1.12.1

- Bug fix for combineAnnotations - MS1 lookup not handling CAMERA adducts correctly
- Fillpeaks bug fix see github #68
- Documentation fixes

# msPurity 1.13.0

- Bioconductor dev (automatic version bump)

# msPurity 1.12.0

- Bioconductor release (automatic version bump)

# msPurity 1.11.5

- frag4feature fileid fix for conversion from factor to character
- Add missing plyr:: reference (thanks jsaintvanne)

# msPurity 1.11.3/4

- Overhaul of combineAnnotation function. Uses local database now as previously
  API calls would take too much time to finish and was not usable
- Various updates of createMSP to make compatible with Galaxy workflows
- Parameter added to purityA to allow user to change the PPM tolerance for MZ
  values between scans when calculated the interpolated precursor ion purity
- Update of spectralMatching results columns to include additional details
  (e.g. retention time)
- Update of spectralMatching so that either PostgreSQL or MySQL database can be
  used as input to either query or library

# msPurity 1.11.2

- Bug fix for EIC with MSMS data

# msPurity 1.11.1

- Bug fix for duplicate MSP spectra when not using metadata table
- Added xcms3 to xcmsSet conversion for "create database" code
- Fix for sirius combine annotations (incorrect column format)

# msPurity 1.11.0

- Bioconductor dev (automatic version bump)

# msPurity 1.10.0

- Bioconductor release (automatic version bump)

# msPurity 1.9.12

- Troublshoot mac "[MSData::Spectrum::getMZIntensityPairs()] Sizes do not match" error

# msPurity 1.9.11

- Bug fix for createMSP - now handles metadata with duplicate grpids

# msPurity 1.9.10

- Documentation updates

# msPurity 1.9.9

- NEW FUNCTION: createDatabase
  - Updated version of create_database that uses an updated schema
- NEW FUNCTION: spectralMatching
  - Updated version spectral_matching that allows for more flexibility
- Vignettes and documentation update

# msPurity 1.9.8

- NEW FUNCTION: filterFragSpectra (for purityA objects)
  - LC-MS/MS filtering of spectra (prior to averaging)
- NEW FUNCTIONS: averageIntraFragSpectra, averageInterFragSpectra, averageAllFragSpectra (for purityA objects)
  - LC-MS/MS averaging and filtering functionality
- NEW FUNCTION: createMSP
  - Create msp files from purityA objects where XCMS features have been linked to fragmentation spectra
- Updated create_database and spectral matching to have the option to use averaged fragmentation spectra

# msPurity 1.9.2

- Bug fix for groupPeaks and groupPeaksEx (the ppm argument was not working and there was inconsistent behaviour with larger datasets). Thanks to Elliot for spotting (https://github.com/litepalmer)
- Updated documentation for spectral matching

# msPurity 1.7.6

- create_sql bug fix for msms spectra table. Was crashing due to incorrect ordering.
- Add peak processing (flag remove peaks) for LC-MS XCMS data
- Add option to for frag4feature to just link to the grouped peak rtmin and rtmax (e.g. do not need to link to individual files)

# msPurity 1.6.1

- Bug fix. For pos/neg switching acquisition two files are can be generated when converting from RAW to mzML (1 for pos, 1 for neg).   The resulting files retention time scans were not being tracked properly in msPurity in these cases. This is now fixed. Thanks
  to Julien (https://github.com/jsaintvanne) for spotting the bug.

# msPurity 1.5.1

- Updates for database creation (can use CAMERA objects now)
- averageSpectra parameter 'MSFileReader' deprecated MSFileReader. Should use csvFile instead, MSFileReader option will still work but a warning will be given

# msPurity 1.4.1/2/3

- Updates for Galaxy for Spectral Matching
- Spectral matching ra_thres_t bugfix
- Separation of sqlite database creation. Now can be called on it's own or with frag4feature
  (allows the Galaxy tool to be simplified)

# msPurity 1.3.9

- Added very basic SIMS stitch compatibility
- pcalc can handle NAs
- Update of purityX to handle obiwarp RT correction (requires recording the RT RAW at an earlier step)
- bug fix for when library spectra is bigger than target spectra (thanks Martin)

# msPurity 1.3.1

- Add spectral matching functionality for LC-MS/MS

# msPurity 1.1.1

- Added pcalc functions to be used by user
- Added option to remove isotopes from calculation

# msPurity 0.99.10/11/12

- Troubleshooting windows build failure

# msPurity 0.99.9

- Offset bug fixed (previous only using extracting lower offset from mzML file)
- Updated handling of RT corrected xcmsSet objects for frag4feature function
- Additional column added for tracking ms/ms spectra

# msPurity 0.99.8

- User option to change the mzR backend library

# msPurity 0.99.4/5/6/7

- Troubleshooting mac build failure

# msPurity 0.99.3

- Grouping multiple peaklist into one wide dataframe
- Peaklists can now be averaged across each class using the function groupPeaks() for the class purityD
- A list of dataframes can also be grouped together using the function groupPeakEx()

# msPurity 0.99.2

- Updated class names purityPD to purityD
- Updated class names purityLC to purityX
- Updated vignette to reflect slightly different terminology
- Added normalised TIC option for purityD

# msPurity 0.99.0

- Initial release (2016-04-08).
