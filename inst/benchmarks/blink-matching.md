# BLINK and default fragment scoring

Eight experimental MS2 queries searched a Parquet library of 128 spectra in a
single process per scorer. Each spectrum contained 64 peaks, giving 512 query
peaks and 8,192 library peaks. Four precursor masses gave 32 eligible
references per query. Both searches scored and retained all 256 eligible pairs.

| Scorer | Search elapsed (s) | Process peak RSS (MiB) |
| --- | ---: | ---: |
| msPurity weighted dot product | 11.594 | 865.5 |
| BLINK approximate cosine | 4.797 | 851.3 |

These are single runs on Intel macOS 26.5.1, with R 4.6.0, msPurity 1.99.4,
Matrix 1.7-5, Spectra 1.22.0, MsBackendParquet 0.99.12 and Arrow 24.0.0.
Elapsed time covers the complete `spectralMatching()` call, including opening
the library, selecting candidates, reading peaks, scoring and ranking. RSS
covers the lifetime of each fresh R process, including package loading. The
operating system's file cache was left enabled.

Both methods used 5 ppm precursor tolerance on each side, zero precursor Da
tolerance, no rank cut and batches of eight references per query. The default
scorer used 10 ppm fragment tolerances, abundance exponent 0.5 and m/z exponent
2. BLINK used 0.001 Da bins, 0.01 Da blur tolerance and intensity exponent 0.5.
The deterministic input generator and search calls are in
[blink-matching.R](blink-matching.R).

## Reproduce

Install this checkout of msPurity, then run from the package directory. The
preparation step creates the library once. Each timed command starts a fresh
R process. On Linux, use GNU `/usr/bin/time -v` and divide the maximum resident
set size in KiB by 1,024. On macOS, divide the byte count by 1,048,576.

```sh
bench_dir=$(mktemp -d)
Rscript inst/benchmarks/blink-matching.R prepare "$bench_dir"
/usr/bin/time -l Rscript inst/benchmarks/blink-matching.R default "$bench_dir"
/usr/bin/time -l Rscript inst/benchmarks/blink-matching.R blink "$bench_dir"
```
