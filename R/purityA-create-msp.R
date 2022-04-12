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



#' @title Using a purityA object, create an MSP file of fragmentation spectra
#'
#' @description
#' **General**
#'
#' Create an MSP file for all the fragmentation spectra that has been linked to an XCMS feature via frag4feature.
#' Can export all the associated scans individually or the averaged fragmentation spectra can be exported.
#'
#' Additional metadata can be included in a dataframe (each column will be added to metadata of the MSP spectra).
#' The dataframe must contain the column "grpid" corresponding to the XCMS grouped feature.
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> **purityA** -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> (optionally) xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'    +  --- *Older versions of XCMS* --- (mzML files) -> xcms.xcmsSet -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> frag4feature -> filterFragSpectra -> averageIntraFragSpectra -> averageIntraFragSpectra -> **createMSP** -> (MSP file)
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
#' @param include_adducts character; Additional adducts to include as a string seperated by white a space (e.g. \[M+H\]+ \[M+Na\]+)
#' @return Returns a MSP file with the selected spectra and metadata
#' @examples
#'
#' #====== XCMS =================================
#' ## Read in MS data
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#' #           package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #ms_data = readMSData(msmsPths, mode = 'onDisk', msLevel. = 1)
#'
#' ## Find peaks in each file
#' #cwp <- CentWaveParam(snthresh = 5, noise = 100, ppm = 10, peakwidth = c(3, 30))
#' #xcmsObj  <- xcms::findChromPeaks(ms_data, param = cwp)
#'
#' ## Optionally adjust retention time
#' #xcmsObj  <- adjustRtime(xcmsObj , param = ObiwarpParam(binSize = 0.6))
#'
#' ## Group features across samples
#' #pdp <- PeakDensityParam(sampleGroups = c(1, 1), minFraction = 0, bw = 30)
#' #xcmsObj <- groupChromPeaks(xcmsObj , param = pdp)
#'
#' #====== msPurity ============================
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa = pa, xcmsObj = xcmsObj)
#' #pa <- filterFragSpectra(pa, allfrag=TRUE)
#' #pa <- averageAllFragSpectra(pa)
#' #createMSP(pa)
#'
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "9_averageAllFragSpectra_with_filter_pa.rds",
#'                           package="msPurity"))
#' createMSP(pa)
#' @md
#' @export
setMethod(f="createMSP", signature="purityA",
          definition = function(pa, msp_file_pth=NULL, metadata=NULL, metadata_cols=NULL,
                                xcms_groupids=NULL, method="all", adduct_split=TRUE, filter=TRUE,
                                msp_schema='massbank', intensity_ra='intensity_ra', include_adducts=''
          ){

            mspurity_to_msp(pa, msp_file_pth, metadata, metadata_cols,
                            xcms_groupids, method, adduct_split, filter, msp_schema,
                            intensity_ra, include_adducts)

          }
)


mspurity_to_msp <- function (pa, msp_file_pth=NULL, metadata=NULL, metadata_cols=NULL,
                             xcms_groupids=NULL, method="all", adduct_split=TRUE, filter=TRUE,
                             msp_schema='massbank', intensity_ra='intensity_ra', include_adducts=''){

  if (!is.null(metadata) && !'grpid' %in% colnames(metadata)){
    metadata$grpid <- 1:nrow(metadata)
  }

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

  of <- file(description = msp_file_pth, open = "w+")
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
            spectrum <- spectrum[spectrum[,'pass_flag']==1,,drop=FALSE]
          }

          if (nrow(spectrum)>0){
            spectrum <- add_mzi_cols(spectrum)
            write.msp(grpdj$precurMtchMZ, grpdj$rt, grpid, fileid, spectrum, metadata,
                      metadata_cols, of, method, adduct_split, msp_schema, intensity_ra, include_adducts)
          }


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

        write.msp(grpdi$precurMtchMZ,grpdi$rt, grpid, fileid, specmax, metadata, metadata_cols, of,
                  method, adduct_split, msp_schema, intensity_ra, include_adducts)

      }else if (method=="av_inter"){

        av_inter <- pa@av_spectra[[as.character(grpid)]]$av_inter

        if (!is.null(av_inter) && length(av_inter)==0){
          next
        }

        if (filter){
          av_inter  <- av_inter[av_inter[,'pass_flag']==1,]
        }


        if (!is.null(av_inter) && nrow(av_inter)>0){
          write.msp(grpd$mz[1], grpd$rt[1], grpid, NA, av_inter, metadata, metadata_cols, of, method,
                    adduct_split, msp_schema, intensity_ra, include_adducts)
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
            write.msp(grpd$mz[1], grpd$rt[1], grpid, fileid, av_intra_j, metadata, metadata_cols,
                      of, method, adduct_split, msp_schema, intensity_ra, include_adducts)
          }

        }

      }else if (method=="av_all"){

        av_all <- pa@av_spectra[[as.character(grpid)]]$av_all

        if (filter){
          av_all  <- av_all[av_all[,'pass_flag']==1,]
        }



        if (!is.null(av_all) && nrow(av_all)>0){

          write.msp(grpd$mz[1], grpd$rt[1], grpid, NA, av_all, metadata, metadata_cols, of,
                    method, adduct_split, msp_schema, intensity_ra, include_adducts)
        }


      }
    }
  }
  close.connection(of)
}




write.msp <- function(precmz, rtmed, grpid, fileid, spectra, metadata, metadata_cols, ofile, method,
                      adduct_split, msp_schema, intensity_ra, include_adducts){

  # Kee the same name for adduct type
  if (msp_schema=='mona'){
    precursor_type <- 'PRECURSOR_TYPE:'
  }else{
    precursor_type <- 'MS$FOCUSED_ION: PRECURSOR_TYPE'
  }

  if ('adduct' %in% names(metadata)){
    # If adduct is used as metadata, we change the name to keep it consistent with the MoNA or MassBank format
    colnames(metadata)[colnames(metadata)=='adduct'] <- precursor_type
  }

  # loop through meta data (if metadata is null the sum will be zero)
  for(i in 1:sum(metadata$grpid==grpid)){
    if(i==0){
      # metadata is null and a single iteration has been performed (all that is required)
      next
    }
    metadatai <- metadata[metadata$grpid==grpid,][i,]


    if (precursor_type %in% names(metadatai) && sum(metadatai$grpid==grpid) > 0){
      # extract the text, expecting to be in CAMERA format, e.g. "[M-H]- 88.016 [M-H-NH3]- 105.042"
      adduct_text <- gsub('[[:space:]]+[[:digit:]]+\\.[[:digit:]]+[[:space:]]*', ' ', metadatai[metadatai$grpid==grpid,precursor_type])

      # add the user provided adduct text
      if(!include_adducts==''){
        adduct_text <- paste(adduct_text, include_adducts, sep = '')
      }


    }else{
      # otherwise just use the user provided text
      adduct_text <- include_adducts
    }

    # split the text into a vecotr
    adducts <- strsplit(adduct_text, ' ')[[1]]



    if (identical(adducts, character(0))){
      adducts <- ''
    }

    # Remove duplicate adducts
    adducts <- unique(adducts)

    if (adduct_split){
      # loop through the adducts creating an appropiate MSP for each
      for (i in 1:length(adducts)){

        adduct <- adducts[i]

        # Already created an output with an adduct if using default adduct
        if (!is.null(metadatai)){
          metadatai[metadatai$grpid==grpid,precursor_type] <- adduct
        }

        write_msp_single(precmz, rtmed, grpid, fileid, adduct, spectra, metadatai, metadata_cols, ofile,
                         method, msp_schema, intensity_ra)
      }

    }else{
      # Ignore adduct splitting
      write_msp_single(precmz, rtmed, grpid, fileid, adduct_text, spectra, metadatai, metadata_cols, ofile,
                       method, msp_schema, intensity_ra)
    }
  }


}



write_msp_single <- function(precmz, rtmed, grpid, fileid, adduct, spectra, metadata, metadata_cols,
                             ofile, method, msp_schema='massbank', intensity_ra='intensity_ra'){


  name <- concat_name(precmz, rtmed, grpid, fileid, adduct, metadata, metadata_cols)

  if (!is.null(metadata) && nrow(metadata[metadata$grpid==grpid,])>0){
    metadata <-  metadata[metadata$grpid==grpid,]
  }

  if(.Platform$OS.type == "unix") {
    line_end = "\n"
  } else {
    line_end = "\r\n"
  }

  if (msp_schema=='mona'){
    if ('NAME:' %in% colnames(metadata)){
      cat(paste0("NAME: ", metadata[,'NAME:'], line_end), file = ofile)
    }else{
      cat(paste0("NAME: ", name, line_end), file = ofile)
    }

    cat(paste0("PRECURSORMZ: ", precmz , line_end), file = ofile)
    cat(paste0("RETENTIONTIME: ", rtmed, line_end), file = ofile)
    if (is.null(metadata) && !adduct==''){
      cat(paste0("PRECURSOR_TYPE: ", adduct, line_end), file = ofile)
    }
  }else{
    if ('RECORD_TITLE:' %in% colnames(metadata)){
      cat(paste0("RECORD_TITLE: ", metadata[,'RECORD_TITLE:'], line_end), file = ofile)
    }else{
      cat(paste0("RECORD_TITLE: ", name, line_end), file = ofile)
    }
    cat(paste0("MS$FOCUSED_ION: PRECURSOR_M/Z ", precmz , line_end), file = ofile)
    cat(paste0("AC$CHROMATOGRAPHY: RETENTION_TIME ", rtmed, line_end), file = ofile)
    if (is.null(metadata) && !adduct==''){
      cat(paste0("MS$FOCUSED_ION: PRECURSOR_TYPE ", adduct, line_end), file = ofile)
    }
  }


  if (!is.null(metadata)){
    metadata_to_write <- metadata[ , !(names(metadata) %in%  c('NAME:','RECORD_TITLE:', 'grpid', 'PRECURSORMZ:', 'RETENTIONTIME:',
                                                               'MS$FOCUSED_ION: PRECURSOR_M/Z', 'AC$CHROMATOGRAPHY: RETENTION_TIME')),
                                   drop=FALSE]

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

concat_name <- function(mz, rtmed, grpid, fileid=NA, adduct, metadata, metadata_cols){
  if (adduct==''){
    adduct = NA
  }

  name <- paste("MZ:" ,round(mz, 4)," | RT:" ,round(rtmed,1),  " | grpid:" ,grpid, " | file:" ,fileid,
                " | adduct:", adduct, sep='')

  if (!is.null(metadata_cols) && !is.null(metadata) && nrow(metadata[metadata$grpid==grpid,])>0){

    name <- paste(paste(metadata_cols, ': ', metadata[metadata$grpid==grpid, metadata_cols], collapse = " | ", sep=''), name, sep=' | ')

  }

  return(name)

}

add_mzi_cols <- function(x){

  x <- data.frame(x)
  colnames(x) <- c('mz', 'i')
  return(x)
}

