# Spectral matching benchmark: SQLite vs mzStack

Results of [spectral-matching.R](spectral-matching.R), run on 2026-10-06. The
benchmark times `spectralMatching()` on the SQLite route (`format = "sqlite"`) and the
mzStack route (`format = "mzstack"`) with the same inputs. It also checks that both
routes return identical scores.

## Setup

| | |
|---|---|
| Machine | Intel Core i7-9750H (6 cores, 12 threads), 16 GB RAM, macOS 26.5.1 |
| Software | R 4.6.0, msPurity 1.37.4, arrow 24.0.0, RSQLite 3.53.1 |
| Query | `extdata/tests/db/createDatabase_example.sqlite`: 23 averaged (`av_all`) spectra and no inter-file averages |
| Mini library | `extdata/tests/mzstack/library/mini_library.sqlite`: 8 spectra |
| Full library | the default msp2db library (Zenodo record 18700802): 201,096 spectra and 16,524,391 peaks |
| Parallelism | `cores = 1` throughout |

The full library holds four sources:

| Source | Spectra |
|---|---|
| lipidblast | 135,456 |
| massbank | 46,334 |
| gnps | 14,686 |
| hmdb | 4,620 |

## Matching

Times are elapsed seconds; where a scenario was run more than once, the median is
shown.

| Scenario | SQLite | mzStack | Speed-up | Matches | Identical scores |
|---|---|---|---|---|---|
| Mini library, `av_all` + `inter` (5 runs each) | 3.87 | 0.37 | 10x | 4 | yes |
| Full library, `av_all` (SQLite 1 run, mzStack 3 runs) | 563.8 | 4.20 | 134x | 123 | yes |
| Full library, `av_all` + `inter` (SQLite 1 run, mzStack 3 runs) | 512.5 | 3.32 | 154x | 123 | yes |

The two full-library scenarios ran the same 23 queries, because the query database has
no inter-file averages. The current script keeps only the first of them.

The individual runs were as follows.

- **Mini library, SQLite:** 2.57, 4.95, 3.87, 4.18, 3.21
- **Mini library, mzStack:** 0.42, 0.33, 0.37, 0.36, 0.51
- **Full library `av_all`, mzStack:** 4.58, 3.94, 4.20
- **Full library `av_all` + `inter`, mzStack:** 3.28, 3.34, 3.32

Writing the matches back took 3.48 s. That run used `updateDb = TRUE` on the mzStack
route against the full library, which appends evidence, scores, compounds and coverage
to the query dataset. The SQLite write was not timed.

### Why the gap grows with library size

On the SQLite route each query spectrum costs about 22 to 25 s against the full library,
so the total grows linearly with the number of queries. The library database has no
index on `library_spectra.library_spectra_meta_id`, so every query scans all 16.5
million peak rows.

The mzStack route works differently:

1. It reads the library metadata once.
2. It takes candidates from a precursor window.
3. It reads peaks only for those candidates.

Its time is mostly that fixed start-up cost, so the difference widens as the number of
query spectra grows.

## One-off conversions and disk size

| Step | Elapsed (s) |
|---|---|
| Query database to mzStack (`convertSqliteToMzstack`) | 2.0 |
| Mini library to mzStack (`convertLibraryToMzstack`) | 0.15 |
| Full library to mzStack (`convertLibraryToMzstack`) | 33.4 (36.6 s CPU, 1.77 GB peak RSS) |
| Full library to mzStack, earlier converter | 81.5 (2.85 GB peak RSS) |

The earlier converter held one data frame per library spectrum and read every peak into R
during validation. The current one produces identical output.

| Dataset | SQLite | mzStack | Ratio |
|---|---|---|---|
| Query | 2.3 MB | 1.1 MB | 2.2x smaller |
| Full library | 502.1 MB | 193.7 MB | 2.6x smaller |

## Caveats

- **Run counts:** each SQLite full-library scenario was run once because of its length,
  so its figure is a single observation.
- **Query set:** there are only 23 queries, all `av_all`. A larger query set would widen
  the gap further, for the reason given above.
- **Cache state:** timings include the operating system's file cache as it was at the
  time. No caches were flushed between runs.
