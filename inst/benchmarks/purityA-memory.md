# purityA workflow: time and memory with and without the legacy slots

`purityA-memory.R` runs `purityA()`, `frag4feature()`,
`filterFragSpectra(allfrag = TRUE)` and the three averaging steps on the two
LC-MS/MS files of msPurityData, first with
`options(msPurity.legacySlots = TRUE)` (the default) and then with `FALSE`.
Results below are from msPurity 1.37.4 on R 4.6.0, Spectra 1.22.0, macOS
(Intel), one core.

## Memory of the final object

| Slot | Legacy slots on | Legacy slots off |
|---|---:|---:|
| `puritydf` | 229 KB | empty |
| `grped_ms2` | 145 KB | empty |
| `all_frag_scans` | 1,916 KB | empty |
| `av_spectra` | 315 KB | empty |
| `grped_df` | 52 KB | 52 KB |
| `spectra` (MsBackendMzR) | 722 KB | 722 KB |
| `fragSpectra` (all MS/MS scans, with flags) | 1,941 KB | 1,941 KB |
| `avSpectra` | 109 KB | 109 KB |
| Whole object | 5,442 KB | 2,838 KB |

`fragSpectra` holds the same data as `all_frag_scans` in about the same
memory. It uses MsBackendDataFrame, because MsBackendMemory keeps each extra
peak variable as a data frame per spectrum and needed 5.5 MB for the same
spectra. With the legacy slots on, both representations are kept, so the
object is about twice the size it is with them off.

## Time per step (seconds)

| Step | Legacy slots on | Legacy slots off |
|---|---:|---:|
| `purityA()` | 14.2 | 11.0 |
| `frag4feature()` | 2.7 | 2.7 |
| `filterFragSpectra(allfrag = TRUE)` | 4.2 | 1.2 |
| `averageIntraFragSpectra()` | 0.35 | 0.32 |
| `averageInterFragSpectra()` | 0.47 | 0.44 |
| `averageAllFragSpectra()` | 0.41 | 0.39 |

With the legacy slots on, `filterFragSpectra(allfrag = TRUE)` also computes
the legacy `all_frag_scans` table the way it always has, which reads the raw
files a second time. Timings vary by a second or so between runs.
