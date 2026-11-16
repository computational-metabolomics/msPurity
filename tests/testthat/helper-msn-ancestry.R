.ancestry_spectra <- function() {
    metadata <- S4Vectors::DataFrame(pid = c(7L, 40L, 41L),
        spectrum_id_ = 101:103, run_id = "run", msLevel = 2:4,
        scanIndex = 1:3, acquisitionNum = 11:13,
        rtime = c(10, 11, 12), precursorMz = c(500, 150, 75),
        precScanNum = c(NA_integer_, 101L, 102L), polarity = 1L,
        analyzer = "Orbitrap", activationMethod = "CID",
        dataOrigin = "example.mzML")
    metadata$mz <- rep(list(100), 3L)
    metadata$intensity <- rep(list(10), 3L)
    Spectra::Spectra(metadata, source = Spectra::MsBackendMemory())
}
