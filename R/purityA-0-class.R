######################################################################
# Create the base Experiment class
######################################################################
# An S4 class to assess precursor purity
#
# Given a vector of LC-MS/MS or DI-MS/MS mzML file paths calculate the precursor purity of each MS/MS scan. See
# purityA constructor function for more information
#
setClass(
  # Set the name for the class
  "purityA",

  # Define the slots
  slots = c(
    fileList = "vector",
    cores = "numeric",
    puritydf = "data.frame",
    grped_df = "data.frame",
    grped_ms2 = "list",
    mzRback = "character",
    db_path = "character",
    f4f_link_type = "character",
    av_spectra = "list",
    av_intra_params = "list",
    av_inter_params = "list",
    av_all_params = "list",
    filter_frag_params = "list"
  )
)

######################################################################
# Show method
######################################################################
#' @title Show method for purityA class
#' @description
#'
#' print statement for purityA class
#' @param object object; purityA object
#' @return a print statement of regarding object
#' @export
setMethod("show", "purityA", function(object) {
  print("purityA object for assessing precursor purity for MS/MS spectra")
})