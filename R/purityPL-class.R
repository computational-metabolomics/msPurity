######################################################################
# Create the base Experiment class
######################################################################
# An S4 class to predict precursor purity for LC-MS data
#
# The class used to predict the purity of an LC-MS dataset
setClass(
  # Set the name for the class
  "purityPL",

  # Define the slots
  slots = c(
    purityType = "character",
    cores = "numeric",
    offsets = "vector",
    fileignore = "vector",
    predictions = "data.frame"
  )
)

######################################################################
# Show method
######################################################################
#' @title Show method for purityPL
#'
#' @description
#' Show method for purityPL object
#'
#'
#' @param object = purityPL object
#' @return a print statement of regarding object
#' @export
setMethod("show", "purityPL", function(object) {
  print("purityPL object for predicting precursor purity from LC-MS spectra")
})
