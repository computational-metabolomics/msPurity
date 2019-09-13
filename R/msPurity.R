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
