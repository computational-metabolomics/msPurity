The daily package check uses the bundled SQLite examples and the eight-spectrum
reference library. The weekly suite compares the complete legacy workflows,
both Spectra storage routes, Parquet and SQLite writers, saved objects and
spectral matches against their historical results. It also exercises full MS2–MS4
readback, XCMS feature detection across both input classes, and averaging across
batch sizes and worker counts.

From the package root, install msPurity and msPurityData (the release containing
`extdata/regression-v1`, first added in 1.41.2), then run
`Rscript --vanilla longtests/longtests.R`. Full-library matching uses the
default BiocFileCache entry and downloads it on first use. Bioconductor runs
this suite through `.BBSoptions`. The Linux pull request job runs it for changes
to the package or its checks.

The companion data package must propagate before enabling the long checks
against the published packages. Its SHA256SUMS records the original files.
To rebuild the daily SQLite examples, run
`python3 tests/testthat/fixtures/make-daily-databases.py /path/to/regression-v1`.
To intentionally regenerate a set of expected results, point
`MSPURITYDATA_SOURCE` at the companion source checkout before running the
corresponding script in `tests/testthat/fixtures`, then review the changes and
update the companion checksums.
