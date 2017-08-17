#' @title Get additional mzML meta
#'
#' @description
#' Extract the filter strings 'accession MS:1000512' from an mzML file. Called header in thermo software.
#' Enables quick access to various information regarding each scan
#' @param mzml_pth
#' @examples
#' mzml_pth <- system.file("extdata", "dims", "mzML", 'B02_Daph_TEST_pos.mzML', package="msPurityData")
#' meta_df <- get_additional_mzml_meta(mzml_pth)
#' @return dataframe of meta info
get_additional_mzml_meta <- function(mzml_pth){

  filters <- parse_file('.*accession=\"MS:1000512\".*name=\"filter string\".*value=\"(.*)\"', mzml_pth)
  scanids <- as.numeric(parse_file('.*spectrum.*index=\".*\".*id=\".*scan=([0-9]+)\".*', mzml_pth))

  if(length(filters)!=length(scanids)){ stop(paste('Scan id and filters (headers) do not match for file', mzml_pth)) }

  sims <- grepl('.*sim.*', filters, ignore.case = TRUE)
  locks <- grepl('.*lock.*', filters, ignore.case = TRUE)
  mass_ranges <- stringr::str_match(filters, '.*\\[([0-9]+\\.[0-9]+)-([0-9]+\\.[0-9]+).*\\]')[,2:3]

  mass_ranges <- apply(mass_ranges, 2, as.numeric)
  colnames(mass_ranges) <- c('scan_window_lower_limit', 'scan_window_upper_limit')


  meta_df <- data.frame('scanid'=scanids, 'header'=filters, 'sim'=sims, mass_ranges, stringsAsFactors = FALSE)
  return(meta_df)
}

parse_file <- function(pattern, inputfile){
  gm <- grep(pattern, readLines(inputfile), value=TRUE)
  sm <- stringr::str_match(gm, pattern)
  return(sm[,2])
}

