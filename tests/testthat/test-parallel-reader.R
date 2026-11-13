test_that("parallel reading preserves file order, headers and peaks", {
    skip_if_not_installed("BiocParallel")
    filenames <- c("msPurity_L_tryptophan_6305_MSn.mzML.gz",
                   "msPurity_L_kynurenine_161166_MSn.mzML.gz")
    files <- file.path(system.file("extdata", "msn-standards",
        package = "msPurity", mustWork = TRUE), filenames)
    reference <- lapply(files, function(path) {
        raw <- mzR::openMSfile(path)
        on.exit(mzR::close(raw))
        list(header = mzR::header(raw), peaks = mzR::peaks(raw))
    })

    previous <- BiocParallel::registered()
    param <- BiocParallel::SnowParam(workers = 2L, type = "SOCK",
                                    fallback = FALSE)
    previous_options <- options(SnowParam = param)
    on.exit({
        options(previous_options)
        for (backend in rev(previous))
            BiocParallel::register(backend, default = TRUE)
        BiocParallel::bpstop(param)
    }, add = TRUE)
    BiocParallel::register(param, default = TRUE)
    BiocParallel::bpstart(param)
    expect_identical(BiocParallel::bpparam(), param)
    worker_pids <- unlist(BiocParallel::bplapply(seq_len(2L),
        function(i) Sys.getpid(), BPPARAM = param), use.names = FALSE)
    expect_length(unique(worker_pids), 2L)
    expect_false(Sys.getpid() %in% worker_pids)

    spectra <- msPurity:::.msp_read(files)

    counts <- vapply(reference, function(x) nrow(x$header), integer(1L))
    expect_identical(Spectra::dataOrigin(spectra),
                     rep(normalizePath(files), counts))
    per_file <- msPurity:::.msp_split_files(spectra)
    for (i in seq_along(files)) {
        expect_identical(msPurity:::.msp_header(per_file[[i]]),
                         reference[[i]]$header, label = filenames[i])
        expect_identical(msPurity:::.msp_peaks(per_file[[i]]),
                         reference[[i]]$peaks, label = filenames[i])
    }
})
