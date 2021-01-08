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
    filter_frag_params = "list",
    all_frag_scans = "data.frame"
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
