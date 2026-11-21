test_that("Parquet export holds configured averages with their provenance", {
    skip_if_no_arrow()
    skip_if_no_mspuritydata()
    skip_if_not_installed("xcms")
    suppressPackageStartupMessages(library(xcms))
    pa <- averageIntraFragSpectra(.golden_pa("3_filterFragSpectra_pa.rds"))
    configured <- averageAllFragSpectra(pa,
        param = AverageFragSpectraParam(combine = "collisionEnergy"))
    xe <- .golden_xcms("msms_only_xcmsnexp.rds")
    out <- tempfile("configured-")
    dir.create(out)
    path <- suppressMessages(createDatabase(configured, xe, outDir = out,
        dbName = "out.parquet", format = "parquet"))
    expect_false(any(validateParquet(path)$level == "MUST"))
    all <- .read_run(path, "av_all")
    expect_identical(sort(all$x_mspurity_spectrum_key),
        sort(configured@avSpectra$spectrumId[
            msPurity:::.pa_av_configured(configured@avSpectra)]))
    expect_true(all(grepl("^msPurity:av_all:[0-9]+:", all$id)))
    ## The intra-file averages computed without param keep their record.
    expect_true(length(grep("^av_intra_", vapply(.manifest(path)$runs,
        `[[`, "", "run_id"))) > 0L)
    mm <- .read_table(path, "merge_member")
    id <- all$spectrum_id_[1L]
    got <- mm[mm$merged_spectrum_id_ == id, ]
    ref <- averageMembers(configured, all$x_mspurity_spectrum_key[1L])
    expect_identical(got$member_rank[order(got$member_rank)], ref$rank)
    expect_identical(got$x_mspurity_contributes[order(got$member_rank)],
        ref$contributes)

    missing <- configured
    missing@params$averageGraph <- NULL
    expect_error(createDatabase(missing, xe, outDir = out, dbName = "m.parquet",
        format = "parquet"), "recorded provenance")
    changed <- configured
    changed@avSpectra <- Spectra::filterIntensity(changed@avSpectra,
        intensity = c(1e6, Inf))
    expect_error(createDatabase(changed, xe, outDir = out, dbName = "c.parquet",
        format = "parquet"), "Peak processing of averaged spectra")
})
