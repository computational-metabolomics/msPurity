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



#' \code{msPurity} package
#'
#' \href{https://bioconductor.org/packages/release/bioc/html/msPurity.html}{msPurity Bioconductor}
#'
#' @docType package
#' @name msPurity
#' @importFrom magrittr %>%
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline legend lines plot points text
#' @importFrom stats approxfun dnorm median na.omit sd
#' @importFrom utils fix head packageVersion read.csv write.csv
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
