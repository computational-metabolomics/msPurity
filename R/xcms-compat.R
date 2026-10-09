# Accessors over the xcms result classes: xcmsSet (xcms 1), XCMSnExp (xcms 3)
# and XcmsExperiment / XcmsExperimentHdf5 (xcms 4).

#' Does `x` support chromPeaks(), featureDefinitions() and featureValues()?
#'
#' @noRd
.xcms_is_modern <- function(x) {
    is(x, "XCMSnExp") || is(x, "XcmsExperiment")
}

#' The raw data files of an xcms result, in sample order.
#'
#' @noRd
.xcms_files <- function(x) {
    if (is(x, "XCMSnExp"))
        return(x@processingData@files)
    if (is(x, "XcmsExperiment"))
        return(xcms::fileNames(x))
    x@filepaths
}

#' The sample class of each sample, or `NULL`.
#'
#' @noRd
.xcms_classes <- function(x) {
    if (is(x, "XCMSnExp"))
        return(x@phenoData@data$class)
    if (is(x, "XcmsExperiment"))
        return(MsExperiment::sampleData(x)$class)
    x@phenoData$class
}

#' The feature names (FT0001, ...).
#'
#' @noRd
.xcms_feature_names <- function(x) {
    if (is(x, "XCMSnExp"))
        return(x@msFeatureData$featureDefinitions@rownames)
    rownames(xcms::featureDefinitions(x))
}

#' Raw or adjusted retention times, as a list by sample.
#'
#' @noRd
.xcms_rtime_by_sample <- function(x, adjusted) {
    if (is(x, "XCMSnExp"))
        return(xcms::rtime(x, adjusted = adjusted, bySample = TRUE))
    sps <- xcms::spectra(x)
    rt <- if (adjusted && "rtime_adjusted" %in% Spectra::spectraVariables(sps))
        sps$rtime_adjusted else Spectra::rtime(sps)
    unname(split(rt, factor(MsExperiment::spectraSampleIndex(x),
                            levels = seq_along(xcms::fileNames(x)))))
}
