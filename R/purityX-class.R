######################################################################
# Create the base Experiment class
######################################################################
# An S4 class to assess the anticipated purity (predicted purity) of XCMS features from an LC-MS run
#
setClass(
  # Set the name for the class
  "purityX",

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
#' @title Show method for purityX
#'
#' @description
#' Show method for purityX object
#'
#'
#' @param object object; purityX object
#' @return a print statement of regarding object
#' @export
setMethod("show", "purityX", function(object) {
  print("purityX object for assessing anticipated precursor ion purity from an XCMS dataset")
})
