# Parquet export memory

The benchmark prepares one XcmsExperiment and purityA object per fixture size.
It alternates the two packaged msPurityData LC–MS/MS files under distinct sample
names, detects and groups chromatographic peaks, retains annotated scan peaks
with `allfrag = TRUE`, and runs intra-file, inter-file and all-file scalar
averaging. Each exporter reads the same saved input in a fresh R process.

| Input runs | Chromatographic peaks | Derived spectra | Bridge peak RSS | Direct peak RSS | Bridge elapsed | Direct elapsed |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 16,618 | 586 | 1,492 MiB | 1,317 MiB | 75.7 s | 55.4 s |
| 4 | 33,236 | 860 | 1,590 MiB | 1,387 MiB | 118.7 s | 99.3 s |

Peak RSS is the operating system's high-water mark for the entire R process,
including package loading, the input objects, Arrow, DuckDB and validation.
Elapsed time covers `createDatabase()`, including validation. These are single
measurements on Intel macOS 26.5.1 with R 4.6.0, xcms 4.10.1, Spectra 1.22.0,
MsBackendParquet 0.99.12, Arrow 24.0.0 and DuckDB 1.5.2. Timing and RSS vary
with the allocator and machine load.

The four-run fixture has 3,316 assessed spectra, 48,954 annotated scan peaks,
12,816 averaging peak edges and 41,006 source identifiers. The direct writer
produces 176,140 rows across 15 tables. The comparison checks all runs, table
columns and values, the results index and the manifest. It sorts table rows
by their declared keys, masks run-specific timestamps and dataset UIDs,
and allows physical part counts to vary. It checks the 132 additional loss
ledger rows against the original scalar averaging results. A second four-run export
with 32 spectra and 1,024 table rows per batch passed the same comparison across
30 native run parts and 197 result table parts.

Export reads up to 128 spectra together and writes result parts of up to 65,536
rows. Peak arrays and the provenance of a single spectrum stay together.
Temporary disk space holds buffered table parts, peak ranks, and the joins used
for feature links and validation. The input objects, spectrum metadata and
lookup indices remain resident, so total process memory also depends on input
size.

## Reproduction

Run the commands from the candidate msPurity checkout. Supply an empty working
directory on a filesystem with space for the mzML copies and both datasets.
The historical checkout supplies the SQLite staging route.

```sh
work=$(mktemp -d)
mkdir "$work/baseline"
git archive 3b1604fd46436bced1d7df9d2f2c528f60444209 | tar -x -C "$work/baseline"
for n in 2 4; do
    Rscript inst/benchmarks/parquet-export.R prepare "$work" "$n"
    /usr/bin/time -l Rscript inst/benchmarks/parquet-export.R bridge "$work" "$n" "$work/baseline"
    /usr/bin/time -l Rscript inst/benchmarks/parquet-export.R direct "$work" "$n" .
    Rscript inst/benchmarks/parquet-export.R compare "$work" "$n"
done
```

On macOS, divide the maximum resident set size in bytes by 1,048,576 for MiB.
GNU time reports maximum resident set size in KiB with `/usr/bin/time -v`.
