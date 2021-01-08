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



#' @title Validate precursor purity predictions using LC-MS and LC-MS/MS dataset
#'
#' @description
#' The method is used to validate the precursor purity predictions made
#' from an LC-MS dataset
#'
#' @aliases validate
#'
#' @param pa object; purityA object
#' @param ppLCMS object; purityX object
#' @return purityA object
#' @export
setMethod(f="validate", signature="purityA",  definition= function(pa, ppLCMS){
  check <- merge(pa@grped_df, ppLCMS@predictions, by = "grpid", sort = FALSE, all.x = TRUE)
  # count how many ms/ms spectra there are per feature
  check$grpid <- as.numeric(check$grpid)

  #msmsperfeature <- median(daply(check,.(grpid), function(x){nrow(x)}))

  check <- check[order(check$grpid),]

  pdiff <- abs(check$median-check$inPurity)

  predictP <- check$median

  head(pa@grped_df)
  pa@grped_df <- cbind(pa@grped_df, predictP, pdiff)

  return(pa)
})

# \dontshow{
# datapth <- system.file("extdata", "lcms", package="msPurityData")
# files <- c("LCMS_1.mzML","LCMS_2.mzML", "LCMSMS_1.mzML", "LCMSMS_2.mzML")
# mzdatafiles <- vector()
# for (i in 1:length(files)){
#    mzdatafiles[i] <- file.path(datapth, "mzML", files[i])
# }
# xset <- readRDS(file.path(datapth, "example_results", "xset.rds"))
# xset@filepaths <- mzdatafiles
# }
# pa <- purityA(mzdatafiles, cores=1, interpol = "linear")
# pa <- frag4feature(pa, xset)
# ppLCMS <- purityPL(xset, fileignore = c(3,4), cores = 1, xgroups = 6)
# pa <- validate(pa, ppLCMS)


