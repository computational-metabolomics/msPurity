#' @import dplyr
#' @import plyr
NULL

#' Calculate m/z range for given m/z value and ppm range
#' 
#' @param ppm Mass error range in ppm
#' @param mass m/z value for which to calculate range
#' @export

mzRange <- function(ppm, mass)
{
  delta <- ppm*mass*10^-6
  out <- c(mass-delta,mass+delta)
  out
}


#' Match list of m/z and/or RT values to msPurity grouped MS/MS spectral data
#' 
#' @param pa msPurity object
#' @param ppm Mass error range in ppm
#' @param rtwin RT window in seconds to use for feature matching
#' @param massTable Data frame containing list of masses to be matched.
#' @param mzColumn Name of the column in massTable for m/z values
#' @param rtColumn Name of the column in massTable for RT values. If set to NULL full RT range os measurement will be used
#' @export

matchTargetIons <- function (pa, ppm=5, rtwin=10, massTable, mzColumn="mz", rtColumn="rt"){
  
  out <- vector("list", nrow(massTable))
  
  for (cmpd in 1: nrow(massTable)){
    
    peakTargetRange <- mzRange(ppm, massTable[cmpd, mzColumn])
    
    if (!is.null(rtColumn)){
      rtTargetRange <- c(massTable[cmpd, rtColumn] - rtwin, massTable[cmpd, rtColumn] + rtwin)
    } else {
      rtTargetRange <- c(min(pa@grped_df$rtmin), max(pa@grped_df$rtmax))
    }
    
    # Match by m/z and RT window, and only for spectra matching filtering criteria, usually msPurity score, snr etc
    peakhits <- which(pa@grped_df$mz >= peakTargetRange[1] 
                      & pa@grped_df$mz <= peakTargetRange[2]
                      & pa@grped_df$rt >= rtTargetRange[1]
                      & pa@grped_df$rt <= rtTargetRange[2]
                      & pa@grped_df$purity_pass_flag ==TRUE)
    
    
    if (length(peakhits)>0) {
      df_out <- as.data.frame(pa@grped_df[peakhits, ] %>%
                                dplyr::group_by(msPurity_grpid = grpid) %>%
                                dplyr::summarise(precursor_mz = mean(mz),
                                                 precursor_rt = mean(rt),
                                                 sample=paste(unique(sample), collapse =", "),
                                                 filename=paste(filename, collapse = ", "),
                                                 cid = paste(unique(cid), collapse =", "),
                                                 pid = paste(unique(pid), collapse =", "))
      )
      
      out[[cmpd]] <- cbind(df_out, massTable[cmpd, ], row.names = NULL)
    }
    
  }
  
  out <- plyr::compact(out)
  out <- do.call (rbind, out)
  out
}
