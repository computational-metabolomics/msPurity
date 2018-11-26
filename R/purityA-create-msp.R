#' @title Create MSP file from purityA object
#'
#' @description
#'
#' Create a MSP file for all the fragmentation spectra that has been linked to an XCMS feature via frag4feature

#' @aliases createMSP
#
#' @param pa object;  purityA object
#' @param msp_file Name of the output msp file
#' @param metadata Data frame with additional coumpound infomation to include in msp output
#' @param metadata_cols Column names of meta data to incorporate in specta annotation
#' @param xcms_groupid XCMS group id's to extract ms/ms data for
#' @param win_format If set to TRUE will use MS Windows line endings, otherwise UNIX format
#' @param method "all" will export all matching ms/ms spectra to xcms features, "max" will use spectra with the highest inensity,
#'               "average-intra" will use the intra file averaged file, "average-inter" will use the inter file (within file)
#'                averaged spectra, "average-all" will use the averaged spectra (ignoring inter and intra)
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths, interpol = "linear")
#' pa <- frag4feature(pa, xset)
#' pa <- averageFragmentation(pa)
#' createMSP(pa)
#' @export
setMethod(f="createMSP", signature="purityA",
          definition = function(pa, msp_file, metadata, metadata_cols=c("Name", "ion"),
                                xcms_groupid=NULL, win_format=FALSE, method="all"){

            mspurity_to_msp(pa, msp_file, metadata, metadata_cols,
                            xcms_groupid, win_format, method)

          }
)

mspurity_to_msp <- function (pa, msp_file_pth, metadata, metadata_cols=c("Name", "ion"),
                             xcms_groupid=NULL, win_format=FALSE, method="all"){


  grped_df <- pa@grped_df
  msms <- pa@grped_ms2
  puritydf <- pa@puritydf


  of <- file(description = msp_file_pth, open = "w+a")
  if (is.null(xcms_groupid)){
    xcms_groupid <- as.numeric(names(pa@grped_ms2))
  }
  for(i in xcms_groupid){

    group_id <- which(grped_df$grpid==i)

    spec <- msms[[as.character(i)]]

    if (length(group_id)>=1){

      grpd <- grped_df[group_id,]

      if (opt$mode=="all"){

        for(jj in 1:length(group_id)){

          grpdj <- grpd[jj,]
          if ('sample' %in% colnames(grpd)){
            fileid = grpdj$sample
          }else{
            fileid = grpdj$fileid
          }

          spectrum <- spec[[jj]]

          name = paste( paste(metadata[xcms_groupid==i, metadata_cols], collapse = ", ")
                        ,", XCMS_group:" ,i, ", file:" ,fileid, sep='')

          write.msp(name,grpdj$precurMtchMZ,"",spectrum,of)
        }

      }else if (opt$mode=="max"){

        prec_int <- sapply(grpd$precurMtchID, function(x) puritydf[ which(puritydf$seqNum==x & puritydf$fileid==grpd$fileid[1]), "precursorIntensity"] )
        idx <- which(prec_int==max(prec_int))
        grpd <- grpd[idx,]

        if ('sample' %in% colnames(grpd)){
          fileid = grpd$sample
        }else{
          fileid = grpd$fileid
        }

        name = paste(i, fileid, grpd$pid, sep='-')
        write.msp(name,grpd$precurMtchMZ,"",specj[[idx]], of)

      }
    }
  }
}




write.msp <- function(name,precmz,prectype,spectra,ofile, win_format=FALSE){

  if (win_format==FALSE){
    line_end = "\n"
  } else {
    line_end = "\r\n"
  }

  cat(paste("NAME: ", name, line_end, sep = ""), file = ofile)

  cat(paste("PRECURSORMZ: ", precmz , line_end, sep = ""), file = ofile)

  cat("Comment:", line_end,file = ofile)

  cat(paste("Num Peaks: ", nrow(spectra), line_end, sep = ""), file = ofile)

  cat(paste(paste(spectra[,1], spectra[,2], sep = "\t"), sep = line_end), sep = line_end, file = ofile)

  cat(line_end, file = ofile)
}

