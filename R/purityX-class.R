# msPurity R package for processing MS/MS data - Copyright (C)
#
# This file is part of msPurity.
#
# msPurity is a free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# msPurity is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with msPurity.  If not, see <https://www.gnu.org/licenses/>.



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
