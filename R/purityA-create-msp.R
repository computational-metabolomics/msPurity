#' @title Create MSP file from purityA object
#'
#' @description
#'
#' Create an MSP file for all the fragmentation spectra that has been linked to an XCMS feature via frag4feature
#'
#' @aliases createMSP
#
#' @param pa object;  purityA object
#' @param msp_file_pth character; Name of the output msp file, if NULL the file "frag_spectra_{time stamp}.msp" will be created in the current directory
#' @param metadata data.frame; Data frame with additional coumpound infomation to include in msp output
#' @param metadata_cols vector; Column names of meta data to incorporate into name
#' @param xcms_groupids vector; XCMS group id's to extract ms/ms data for
#' @param method character; "all" will export all matching ms/ms spectra to xcms features, "max" will use spectra with the highest inensity,
#'               "av_intra" will use the intra file averaged spectra (within file), "av_inter" will use the inter file (across file)
#'                averaged spectra, "av_all" will use the averaged spectra (ignoring inter and intra)
#' @param adduct_split boolean; If either "adduct" or  MS$FOCUSED_ION: PRECURSOR_TYPE column is in metadata then each adduct will have it's own MSP spectra.
#'                     (Useful, if the MSP file will be used for further annotation)
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- frag4feature(pa, xset)
#' pa <- averageAllFragSpectra(pa)
#' createMSP(pa)
#' @export
#setMethod(f="createMSP", signature="purityA",
 #         definition = 
          createMSP <- function(pa, msp_file_pth=NULL, metadata=NULL, metadata_cols=c("CH$NAME", "MS$FOCUSED_ION: PRECURSOR_TYPE"),
                                xcms_groupids=NULL, method="all", adduct_split=TRUE){

cat("Building your MSP file...\n")
            mspurity_to_msp(pa, msp_file_pth, metadata, metadata_cols,
                            xcms_groupids, method, adduct_split)

          }
#)


mspurity_to_msp <- function (pa, msp_file_pth=NULL, metadata=NULL, metadata_cols=c("CH$NAME", "MS$FOCUSED_ION: PRECURSOR_TYPE"),
                             xcms_groupids=NULL, method="all", adduct_split){
print(paste("Methode :",method))
  if (is.null(msp_file_pth)){
    msp_file_pth <- paste('frag_spectra', format(Sys.time(), "%Y-%m-%d-%I%M%S"), '.msp', sep="")
  }

  grped_df <- pa@grped_df
  msms <- pa@grped_ms2
  puritydf <- pa@puritydf

  of <- file(description = msp_file_pth, open = "w+a")
  if (is.null(xcms_groupids)){
    xcms_groupids <- as.numeric(names(pa@grped_ms2))
  }

  for(grpid in xcms_groupids){

    group_id <- which(grped_df$grpid==grpid)

    specnum <- which(names(msms)==grpid)

    spec <- NULL
    for(n in 1:length(specnum)){
      spec <- c(spec,msms[[specnum[n]]])
    }

    if (length(group_id)>=1){

      grpd <- grped_df[group_id,]

      if (method=="all"){
print("method all")
        for(j in 1:length(group_id)){

          grpdj <- grpd[j,]
          print(grpdj)
          if ('sample' %in% colnames(grpd)){
            fileid = grpdj$sample
          }else{
            fileid = grpdj$fileid
          }

          spectrum <- spec[[j]]
          spectrum <- add_mzi_cols(spectrum)

          if(pa@f4f_link_type=='individual'){
            write.msp(grpdj$precurMtchMZ, grpdj$rt, grpid, fileid, spectrum, metadata,metadata_cols, of, method, adduct_split)
          }else{
            write.msp(grpdj$precurMtchMZ, grpdj$rtmed, grpid, fileid, spectrum, metadata,metadata_cols, of, method, adduct_split)
          }


        }

      }else if (method=="max"){
print("method max")

        prec_int <- puritydf[puritydf$pid %in% grpd$pid,'precursorIntensity']
        idx <- which(prec_int==max(prec_int))[1]  # if joint place, take the first one (very unlikely to occur)

        grpdi <- grpd[idx,]

        if ('sample' %in% colnames(grpdi)){
          fileid = grpdi$sample
        }else{
          fileid = grpdi$fileid
        }

        spec_max<- add_mzi_cols(spec[[idx]])

        write.msp(grpdi$precurMtchMZ,grpdi$rt, grpid, fileid, spec_max, metadata, metadata_cols, of, method, adduct_split)

      }else if (method=="av_inter"){
print("method av_inter")

        av_inter <- pa@av_spectra[[as.character(grpid)]]$av_inter

        if (!is.null(av_inter) && nrow(av_inter)>0){
          write.msp(grpd$mz[1], grpd$rt[1], grpid, NA, av_inter, metadata, metadata_cols, of, method, adduct_split)
        }


      }else if (method=="av_intra"){
print("method av_intra")

        av_intra <- pa@av_spectra[[as.character(grpid)]]$av_intra
        if (!is.null(av_intra) && length(av_intra)==0){
          next
        }


        for (j in 1:length(av_intra)){
          av_intra_j <- av_intra[j]
          fileid <- names(av_intra_j)
          if (!is.null(av_intra_j) && nrow(av_intra_j[[1]])>0){
            write.msp(grpd$mz[1], grpd$rt[1], grpid, fileid, av_intra_j[[1]], metadata, metadata_cols, of, method, adduct_split)
          }

        }

      }else if (method=="av_all"){
print("method av_all")

        av_all <- pa@av_spectra[[as.character(grpid)]]$av_all

        if (!is.null(av_all) && nrow(av_all)>0){
          write.msp(grpd$mz[1], grpd$rt[1], grpid, NA, av_all, metadata, metadata_cols, of, method, adduct_split)
        }


      }
    }
  }
  close.connection(of)
}




write.msp <- function(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method, adduct_split){

  if (adduct_split & sum(metadata$grpid==grpid) > 0){
    # To keep the naming consisten we stick with the massbank naming convention
    if ('adduct' %in% names(metadata)){
      names(metadata)[names(metadata)=='adduct'] <- 'MS$FOCUSED_ION: PRECURSOR_TYPE'
    }

    # check if in the precursor_type is in the columns
    if ('MS$FOCUSED_ION: PRECURSOR_TYPE' %in% names(metadata)){

      # extract the text, expecting to be in CAMERA format, e.g. "[M-H]- 88.016 [M-H-NH3]- 105.042"
      adduct_text <- gsub('[[:space:]]+[[:digit:]]+\\.[[:digit:]]+[[:space:]]*', ' ', metadata[metadata$grpid==grpid,'MS$FOCUSED_ION: PRECURSOR_TYPE'])
      # get a vector of the adducts

      adducts <- strsplit(adduct_text, ' ')[[1]]
      # loop through the adducts creating an appropiate MSP for each
      for (i in 1:length(adducts)){
        adduct <- adducts[i]
        metadata[metadata$grpid==grpid,'MS$FOCUSED_ION: PRECURSOR_TYPE'] = adduct
        write_msp_single(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method)
      }

    }else{

      # adduct set to split but adduct column not available
      write_msp_single(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method)
    }

  }else{

    # Ignore adduct splitting
    write_msp_single(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method)
  }

}

write_msp_single <- function(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method, msp_schema='massbank'){

  name <- concat_name(precmz, rtmed, grpid, fileid, metadata, metadata_cols)

  if (!is.null(metadata) && nrow(metadata[metadata$grpid==grpid,])>0){
    metadata <-  metadata[metadata$grpid==grpid,]
  }

  if(.Platform$OS.type == "unix") {
    line_end = "\n"
  } else {
    line_end = "\r\n"
  }
  if (!msp_schema=='massbank'){
    cat(paste("NAME: ", name, line_end, sep = ""), file = ofile)
    cat(paste("PRECURSORMZ: ", precmz , line_end, sep = ""), file = ofile)
  }

  cat(paste("CH$NAME: ", name, line_end, sep = ""), file = ofile)
  cat(paste("MS$FOCUSED_ION: PRECURSOR_M/Z ", precmz , line_end, sep = ""), file = ofile)
  cat(paste("AC$CHROMATOGRAPHY: RETENTION_TIME ", rtmed, line_end, sep = ""), file = ofile)

  if (!is.null(metadata)){
    cat(paste(names(metadata), unlist(metadata), line_end), sep="", file = ofile)
  }



  cat(paste("COMMENT: Exported from msPurity purityA object using function createMSP, using method '",
            method, "' msPurity version:", packageVersion("msPurity"), sep=''),
      line_end,file = ofile)

  cat(paste("PK$NUM_PEAK: ", nrow(spectra), line_end, sep = ""), file = ofile)
  cat(paste("PK$PEAK: m/z int. rel.int.", line_end, sep = ""), file = ofile)

  cat(paste(paste(spectra$mz, spectra$i,  round(spectra$i/max(spectra$i)*100,2), sep = "\t"), sep = line_end),
      sep = line_end, file = ofile)

  cat(line_end, file = ofile)

}

concat_name <- function(mz, rtmed, grpid, fileid=NA, metadata, metadata_cols){

  name <- paste(" MZ:" ,round(mz, 4)," | RT:" ,round(rtmed,1),  " | XCMS_group:" ,grpid, " | file:" ,fileid, sep='')

  if (!is.null(metadata) && nrow(metadata[metadata$grpid==grpid,])>0){

    name <- paste(paste(metadata[metadata$grpid==grpid, metadata_cols], collapse = " | "), name, sep=' | ')

  }

  return(name)

}

add_mzi_cols <- function(x){
  x <- data.frame(x)
  colnames(x) <- c('mz', 'i')
  return(x)
}

