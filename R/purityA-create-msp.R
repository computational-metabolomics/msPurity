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
#' @param filter boolean; TRUE if filtered peaks are to be removed
#' @param msp_schema character; Either MassBank (Europe) or MoNA style of MSP file format to be used ('massbank' or 'mona')
#' @param intensity_ra character; Either 'intensity', 'ra' (relative abundance) or 'intensity_ra' (intensity and relative abundance) to be written
#'                             to the MSP file
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
setMethod(f="createMSP", signature="purityA",
          definition = function(pa, msp_file_pth=NULL, metadata=NULL, metadata_cols=NULL,
                                xcms_groupids=NULL, method="all", adduct_split=TRUE, filter=TRUE,
                                msp_schema='massbank', intensity_ra='intensity_ra'){

            mspurity_to_msp(pa, msp_file_pth, metadata, metadata_cols,
                            xcms_groupids, method, adduct_split, filter, msp_schema,
                            intensity_ra)

          }
)


mspurity_to_msp <- function (pa, msp_file_pth=NULL, metadata=NULL, metadata_cols=NULL,
                             xcms_groupids=NULL, method="all", adduct_split=TRUE, filter=TRUE, msp_schema='massbank', intensity_ra='intensity_ra'){

  if (!msp_schema %in% c('massbank', 'mona')){
    stop('msp_schema should either be "massbank" or "mona"')
  }

  if (!intensity_ra %in% c('intensity', 'ra', 'intensity_ra')){
    stop('intensity_ra should either be "intensity", "ra" or "intensity_ra"')
  }


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

    spec <- msms[[as.character(grpid)]]


    if (length(group_id)>=1){

      grpd <- grped_df[group_id,]

      if (method=="all"){

        for(j in 1:length(group_id)){

          grpdj <- grpd[j,]
          if ('sample' %in% colnames(grpd)){
            fileid = grpdj$sample
          }else{
            fileid = grpdj$fileid
          }

          spectrum <- spec[[j]]

          if ((filter)  & ('pass_flag' %in% colnames(spectrum))){
             spectrum <- spectrum[spectrum[,'pass_flag']==1,]
          }

          spectrum <- add_mzi_cols(spectrum)


          write.msp(grpdj$precurMtchMZ, grpdj$rt, grpid, fileid, spectrum, metadata,metadata_cols, of, method, adduct_split, msp_schema, intensity_ra)


        }

      }else if (method=="max"){

        prec_int <- puritydf[puritydf$pid %in% grpd$pid,'precursorIntensity']
        idx <- which(prec_int==max(prec_int))[1]  # if joint place, take the first one (very unlikely to occur)

        grpdi <- grpd[idx,]

        if ('sample' %in% colnames(grpdi)){
          fileid = grpdi$sample
        }else{
          fileid = grpdi$fileid
        }
        specmax <- spec[[idx]]

        if ((filter) & ('pass_flag' %in% colnames(specmax))){
          specmax <- specmax[specmax[,'pass_flag']==1,]
        }
        specmax<- add_mzi_cols(specmax)


        write.msp(grpdi$precurMtchMZ,grpdi$rt, grpid, fileid, specmax, metadata, metadata_cols, of, method, adduct_split, msp_schema, intensity_ra)

      }else if (method=="av_inter"){

        av_inter <- pa@av_spectra[[as.character(grpid)]]$av_inter

        if (!is.null(av_inter) && length(av_inter)==0){
          next
        }

        if (filter){
          av_inter  <- av_inter[av_inter[,'pass_flag']==1,]
        }


        if (!is.null(av_inter) && nrow(av_inter)>0){
          write.msp(grpd$mz[1], grpd$rt[1], grpid, NA, av_inter, metadata, metadata_cols, of, method, adduct_split, msp_schema, intensity_ra)
        }


      }else if (method=="av_intra"){
        av_intra <- pa@av_spectra[[as.character(grpid)]]$av_intra


        if (!is.null(av_intra) && length(av_intra)==0){
          next
        }


        for (j in 1:length(av_intra)){
          av_intra_j_l <- av_intra[j]
          fileid <- names(av_intra_j_l)
          av_intra_j <- av_intra_j_l[[1]]
          if (filter){
            av_intra_j  <- av_intra_j[av_intra_j[,'pass_flag']==1,]
          }


          if (!is.null(av_intra_j) && nrow(av_intra_j)>0){
            write.msp(grpd$mz[1], grpd$rt[1], grpid, fileid, av_intra_j, metadata, metadata_cols, of, method, adduct_split, msp_schema, intensity_ra)
          }

        }

      }else if (method=="av_all"){

        av_all <- pa@av_spectra[[as.character(grpid)]]$av_all

        if (filter){
           av_all  <- av_all[av_all[,'pass_flag']==1,]
        }

        if (!is.null(av_all) && nrow(av_all)>0){
          write.msp(grpd$mz[1], grpd$rt[1], grpid, NA, av_all, metadata, metadata_cols, of, method, adduct_split, msp_schema, intensity_ra)
        }


      }
    }
  }
  close.connection(of)
}




write.msp <- function(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method, adduct_split, msp_schema, intensity_ra){

  if (msp_schema=='mona'){
    precursor_name <- 'PRECURSOR_TYPE:'
  }else{
    precursor_name <- 'MS$FOCUSED_ION: PRECURSOR_TYPE'
  }

  if (is.null(metadata_cols)){
    if (msp_schema=='mona'){
      metadata_cols <- c("NAME:", "PRECURSOR_TYPE:")
    }else{
      metadata_cols <- c("RECORD_TITLE:", "MS$FOCUSED_ION: PRECURSOR_TYPE")
    }
  }

  if (adduct_split & sum(metadata$grpid==grpid) > 0){
    # To keep the naming consisten we stick with the massbank naming convention

    if ('adduct' %in% names(metadata)){
      names(metadata)[names(metadata)=='adduct'] <- precursor_name
    }

    # check if in the precursor_type is in the columns
    if (precursor_name %in% names(metadata)){

      # extract the text, expecting to be in CAMERA format, e.g. "[M-H]- 88.016 [M-H-NH3]- 105.042"
      adduct_text <- gsub('[[:space:]]+[[:digit:]]+\\.[[:digit:]]+[[:space:]]*', ' ', metadata[metadata$grpid==grpid,precursor_name])
      # get a vector of the adducts

      adducts <- strsplit(adduct_text, ' ')[[1]]
      # loop through the adducts creating an appropiate MSP for each
      for (i in 1:length(adducts)){
        adduct <- adducts[i]
        metadata[metadata$grpid==grpid,precursor_name] = adduct
        write_msp_single(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method, msp_schema, intensity_ra)
      }

    }else{
      # adduct set to split but adduct column not available
      write_msp_single(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method, msp_schema, intensity_ra)
    }

  }else{
    # Ignore adduct splitting
    write_msp_single(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method, msp_schema, intensity_ra)
  }

}

write_msp_single <- function(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method, msp_schema='massbank', intensity_ra='intensity_ra'){
  name <- concat_name(precmz, rtmed, grpid, fileid, metadata, metadata_cols)

  if (!is.null(metadata) && nrow(metadata[metadata$grpid==grpid,])>0){
    metadata <-  metadata[metadata$grpid==grpid,]
  }

  if(.Platform$OS.type == "unix") {
    line_end = "\n"
  } else {
    line_end = "\r\n"
  }

  if (msp_schema=='mona'){
    cat(paste0("NAME: ", name, line_end), file = ofile)
    cat(paste0("PRECURSORMZ: ", precmz , line_end), file = ofile)
    cat(paste0("RETENTIONTIME: ", rtmed, line_end), file = ofile)
  }else{
    cat(paste0("RECORD_TITLE: ", name, line_end), file = ofile)
    cat(paste0("MS$FOCUSED_ION: PRECURSOR_M/Z ", precmz , line_end), file = ofile)
    cat(paste0("AC$CHROMATOGRAPHY: RETENTION_TIME ", rtmed, line_end), file = ofile)
  }




  if (!is.null(metadata)){
    metadata_to_write <- metadata[ , !(names(metadata) %in%  c('NAME:','RECORD_TITLE:', 'grpid', 'PRECURSORMZ:', 'RETENTIONTIME:',
                                                               'MS$FOCUSED_ION: PRECURSOR_M/Z', 'AC$CHROMATOGRAPHY: RETENTION_TIME'))]

    cat(paste(as.character(names(metadata_to_write)), ' ', as.character(unlist(metadata_to_write)), line_end, sep='', collapse=''), file = ofile)
  }


  cat(paste0("XCMS groupid (grpid): ", grpid, line_end), file = ofile)
  cat(paste0("COMMENT: Exported from msPurity purityA object using function createMSP, using method '",
            method, "' msPurity version:", packageVersion("msPurity")),
      line_end,file = ofile)

  if (msp_schema=='mona'){
    cat(paste0("Num Peaks: ", nrow(spectra), line_end), file = ofile)
  }else{
    cat(paste0("PK$NUM_PEAK: ", nrow(spectra), line_end), file = ofile)

  }

  if (intensity_ra=='intensity_ra'){
    if(msp_schema=='massbank'){
      cat(paste0("PK$PEAK: m/z int. rel.int.", line_end), file = ofile)
    }
    cat(paste(paste(spectra$mz, spectra$i,  round(spectra$i/max(spectra$i)*100,2), sep = "\t"), sep = line_end),
        sep = line_end, file = ofile)
  }else if (intensity_ra == 'ra'){
    if(msp_schema=='massbank'){
      cat(paste0("PK$PEAK: m/z rel.int.", line_end), file = ofile)
    }

    cat(paste(paste(spectra$mz, round(spectra$i/max(spectra$i)*100,2), sep = "\t"), sep = line_end),
        sep = line_end, file = ofile)
  }else{
    if(msp_schema=='massbank'){
      cat(paste0("PK$PEAK: m/z int.", line_end), file = ofile)
    }

    cat(paste(paste(spectra$mz, spectra$i, sep = "\t"), sep = line_end),
        sep = line_end, file = ofile)
  }


  cat(line_end, file = ofile)

}

concat_name <- function(mz, rtmed, grpid, fileid=NA, metadata, metadata_cols){

  name <- paste(" MZ:" ,round(mz, 4)," | RT:" ,round(rtmed,1),  " | grpid:" ,grpid, " | file:" ,fileid, sep='')

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

