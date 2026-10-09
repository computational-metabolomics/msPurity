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
    all_frag_scans = "data.frame",
    # Arguments of purityA() and frag4feature(), for provenance.
    params = "list",
    # MS/MS scans with their purity results as spectra variables.
    spectra = "Spectra",
    # Scans used downstream, each once, with filter flags as peak variables.
    fragSpectra = "Spectra",
    # Averaged spectra, one per feature, averaging level and sample.
    avSpectra = "Spectra"
  ),
  prototype = prototype(
    spectra = Spectra::Spectra(Spectra::MsBackendMemory()),
    fragSpectra = Spectra::Spectra(Spectra::MsBackendMemory()),
    avSpectra = Spectra::Spectra(Spectra::MsBackendMemory())
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
  if (.pa_is_current(object)) {
    cat(length(object@spectra), "MS/MS scans,",
        length(unique(object@grped_df$grpid)), "features with scans,",
        length(object@avSpectra), "averaged spectra\n")
  }
})
