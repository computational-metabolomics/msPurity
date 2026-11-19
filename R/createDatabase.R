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



#' @title Create database
#'
#' @description
#'
#' **General**
#'
#' Create an SQLite database of an LC-MS(/MS) experiment (replaces the create_database function).
#'
#' Schema details can be found [here](https://bioconductor.org/packages/release/bioc/vignettes/msPurity/inst/doc/msPurity-spectral-datatabase-schema.html).
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.findChromPeaks -> (optionally) xcms.adjustRtime -> xcms.groupChromPeaks -> (xcmsObj)
#'    +  --- *Older versions of XCMS* --- (mzML files) -> xcms.xcmsSet -> xcms.group -> xcms.retcor -> xcms.group -> (xcmsObj)
#'  * Fragmentation processing
#'    + (xcmsObj, pa) -> frag4feature -> filterFragSpectra -> averageAllFragSpectra -> **createDatabase** -> spectralMatching -> (sqlite spectral database)
#'
#' @param pa purityA object; Needs to be the same used for frag4feature function
#' @param xcmsObj xcms object of class XcmsExperiment, XCMSnExp or xcmsSet; Needs to be the same used for frag4feature function (this will be ignored when using xsa parameter)
#' @param xsa CAMERA object (optional); if CAMERA object is used, we ignore the xset parameter input and obtain all information
#'                          from the xset object nested with the CAMERA xsa object. Adduct and isotope information
#'                          will be included into the database when using this parameter. The underlying xset object must
#'                          be the one used for the frag4feature function
#' @param dbName character (optional); Name of the result database
#' @param grpPeaklist dataframe (optional); Can use any peak dataframe. Still needs to be derived from the xset object though
#' @param outDir character; Out directory for the SQLite result database
#' @param metadata list; A list of metadata to add to the s_peak_meta table
#' @param xset xcms object of class XCMSnExp or xcmsSet; (Deprecated - if provided, will replace variable 'obj')
#' @param format character; \code{"sqlite"} (the default) writes the SQLite database described above, and raises a
#'               deprecation warning: SQLite output is frozen and receives bug fixes only. \code{"parquet"} writes a
#'               Parquet dataset instead (experimental). See Details.
#' @param study character (optional, parquet only); path of the Parquet dataset the mzML files were converted into with
#'              MsBackendParquet::mzMLToParquet(). Scans and their parent scans are then referenced by their spectrum ids
#'              in that dataset. By default, when \code{purityA()} was given a Spectra on that dataset, it is the study;
#'              otherwise, and with \code{study = NA}, every source file is recorded as an external source.
#' @param studyKey character (parquet only); the key the study is declared under in the results dataset.
#' @param fileMap named character (optional, parquet only); source file path (as in \code{pa@fileList}) -> run_id.
#'                Required when two source files share a file name, and to match files whose recorded paths
#'                no longer exist.
#' @param overwrite logical (parquet only); replace an existing dataset at the destination.
#' @return For \code{format = "sqlite"}, the path to the SQLite database. For \code{format = "parquet"}, the path of
#'         the Parquet results dataset, invisibly.
#'
#' @details
#' **Parquet output (experimental)**
#'
#' With \code{format = "parquet"}, the results are written as a results dataset (a directory of Parquet
#' files with an \code{mzStack.json} manifest), following part 4 of the mzStack specification:
#'
#'  * xcms features, chromatographic peaks and their per-sample abundances become the \code{feature},
#'    \code{chromatographic_peak}, \code{chromatographic_peak_feature} and \code{abundance} tables, with
#'    \code{assay} and \code{sample} for the source files;
#'  * the links between MS2 scans and features become \code{spectrum_feature}, with the precursor m/z,
#'    error, retention time and purity each scan was matched at and how its retention time was mapped;
#'    MS3 and MS4 scans linked through their MS2 scan are in \code{x_mspurity_msn_feature};
#'  * averaged fragmentation spectra become spectra of the dataset, one run per source file for spectra
#'    averaged within a file (\code{av_intra_<run>}) and one run each for spectra averaged across files
#'    (\code{av_inter}, \code{av_all}), with their per-peak statistics;
#'  * the scans each averaged spectrum was made from are recorded in \code{merge_member};
#'  * every assessed MS2 to MS4 scan is listed in \code{x_mspurity_scan} with its MS2 precursor purity or
#'    MSn stage purity, the scan it fragments and the MS2 scan its path starts from, and the activation
#'    method msPurity used with where that came from; peaks flagged by \code{filterFragSpectra()} are in
#'    \code{x_mspurity_scan_peak}. Raw scans are referenced, not copied;
#'  * the processing parameters and history, what was converted and anything that could not be mapped are
#'    recorded in \code{tool_provenance}, \code{conversion}, \code{source_identifier} and \code{loss_ledger}.
#'
#' **Raw data in a Parquet study**
#'
#' Raw MS1-MS4 files converted with \code{MsBackendParquet::mzMLToParquet()} are read by passing
#' \code{purityA()} a \code{Spectra} on the dataset; msPurity then reads their peaks from the dataset, never
#' from the mzML files. A scan is linked to the study by the dataset's own key when \code{purityA()} read it
#' there, otherwise by native id, then by scan number at the same position, all from the study's own metadata.
#' An MS2 scan's parent is the study's MS1 scan its precursor scan number names. Acquisition conditions the
#' study does not record, such as the activation method of mzXML input, are reported in the loss ledger as
#' unavailable rather than inferred; where msPurity used a method named by a scan's filter string,
#' \code{activation_type_basis} says so.
#'
#' The dataset is written to \code{file.path(outDir, dbName)}; the default name ends in \code{.parquet}. It is built
#' beside its destination and moved into place only once complete, so an interrupted write leaves nothing behind.
#'
#' CAMERA annotation (\code{xsa}) and a custom \code{grpPeaklist} are not mapped to the Parquet output and are refused.
#'
#' @examples
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
#'
#' ## Run from previously generated data
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#'            package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'               "9_averageAllFragSpectra_with_filter_pa.rds", package="msPurity"))
#' xcmsObj <- readRDS(system.file("extdata", "tests", "xcms",
#'                                "msms_only_xcmsnexp.rds", package="msPurity"))
#' ## The saved objects record the file paths of the machine that made them
#' pa@fileList[1:2] <- msmsPths
#' xcmsObj@processingData@files <- msmsPths
#' db_pth <- createDatabase(pa, xcmsObj, outDir = tempdir(),
#'     metadata = list(polarity = "positive", instrument = "Q-Exactive"))
#'
#' @md
#' @details The SQLite database holds MS2 results only and refuses MS3/MS4
#'   feature links and processed spectra; the Parquet dataset holds MS2-MS4
#'   scans, links and filtered spectra. Averaged MS3/MS4 spectra cannot be
#'   exported yet.
#' @export
createDatabase <-  function(pa, xcmsObj, xsa=NULL, outDir='.', grpPeaklist=NA, dbName=NA, metadata=NA, xset = NA,
                            format = c("sqlite", "parquet"), study = NULL, studyKey = "study",
                            fileMap = NULL, overwrite = FALSE){
  pa <- .pa_update(pa)
  format <- match.arg(format)
  if (format == "parquet"){
    .pa_require_msn_export(pa)
    return(.createDatabase_parquet(pa = pa, xcmsObj = xcmsObj, xsa = xsa, outDir = outDir,
                                   grpPeaklist = grpPeaklist, dbName = dbName,
                                   metadata = metadata, xset = xset, study = study,
                                   studyKey = studyKey, fileMap = fileMap,
                                   overwrite = overwrite))
  }
  .pa_require_ms2_export(pa)
  db <- .createDatabase_sqlite(pa = pa, xcmsObj = xcmsObj, xsa = xsa, outDir = outDir,
                               grpPeaklist = grpPeaklist, dbName = dbName,
                               metadata = metadata, xset = xset)
  .msp_deprecate_sqlite("createDatabase()")
  db
}

# The SQLite writer behind createDatabase(), without the deprecation
# warning, for callers that warn once themselves or not at all.
.createDatabase_sqlite <- function(pa, xcmsObj, xsa=NULL, outDir='.', grpPeaklist=NA, dbName=NA, metadata=NA, xset = NA){
  ########################################################
  # Export the target data into sqlite database
  ########################################################

  if(!is.na(xset)){
    xcmsObj = xset
  }

  #if dbName is not defined, automatically generate a name
  if (is.na(dbName)){
    dbName <- paste('lcmsms_data', format(Sys.time(), "%Y-%m-%d-%I%M%S"), '.sqlite', sep="-")
  }

  #if a peaklist was not supplied to function, extract from xsa or obj.
  if (!is.data.frame(grpPeaklist)){
    grpPeaklist <- .db_peaklist(xcmsObj, xsa)
  }

  message("Creating a database of fragmentation spectra and LC features")
  targetDBpth <- export2sqlite(pa = pa, grpPeaklist = grpPeaklist, xcmsObj = xcmsObj,
                               xsa=xsa, outDir=outDir, dbName=dbName, metadata=metadata)

  return(targetDBpth)

}


# The feature table of `xcmsObj` (or of the CAMERA object `xsa`), with a
# leading grpid column.
.db_peaklist <- function(xcmsObj, xsa = NULL) {
  if (is.null(xsa)){
    if(.xcms_is_modern(xcmsObj)){
      grpPeaklist <- cbind(xcms::featureDefinitions(xcmsObj), xcms::featureValues(xcmsObj))
    }else if(is(xcmsObj, 'xcmsSet')){
      grpPeaklist <- xcms::peakTable(xcmsObj)
    }else{
      stop('createDatabase stopped as the "xcmsObj" (or "xset" if specified) argument is not of class "XcmsExperiment", "XCMSnExp" or "xcmsSet"')
    }
  }else{
    grpPeaklist <- CAMERA::getPeaklist(xsa)
  }
  data.frame(cbind('grpid'=1:nrow(grpPeaklist), grpPeaklist))
}

export2sqlite <- function(pa, grpPeaklist, xcmsObj, xsa, outDir, dbName, metadata){
  prep <- .db_prepare(pa = pa, grpPeaklist = grpPeaklist, xcmsObj = xcmsObj,
                      xsa = xsa, metadata = metadata)
  if (is.null(prep))
    return(NULL)
  dbPth <- file.path(outDir, dbName)
  con <- DBI::dbConnect(RSQLite::SQLite(),dbPth)
  for (t in prep$tables)
    custom_dbWriteTable(name_pk = t$pk, fks = t$fks, table_name = t$name,
                        df = t$df, con = con, pk_type = t$pk_type)
  DBI::dbDisconnect(con)
  return(dbPth)
}

# The tables createDatabase() writes, in write order, each with its primary
# key, foreign keys and data. NULL, with a message, when the files of `pa`
# and `xcmsObj` do not match.
#
# scanPeaks = "count" leaves out the raw scan peaks that are not annotated
# (allfrag = FALSE), which the SQLite writer re-reads from the mzML files,
# and records only how many there are, in `raw_scan_peaks`.
.db_prepare <- function(pa, grpPeaklist, xcmsObj, xsa, metadata,
                        scanPeaks = c("read", "count")){
  scanPeaks <- match.arg(scanPeaks)
  tables <- list()
  add <- function(name, pk, fks, df, pk_type = 'INTEGER')
    tables[[name]] <<- list(name = name, pk = pk, fks = fks, df = df,
                            pk_type = pk_type)
  raw_scan_peaks <- NULL

  if(!is.null(xsa)){
    # if user has supplied camera object we use the xset that the camera object
    # is derived from
    xcmsObj <- xsa@xcmsSet
    XCMSnExp_bool <- FALSE
  }else{
    #confirm whether xcmsObj is of class "XCMSnExp" (if not, it should be of class "xcmsSet")
    XCMSnExp_bool <- .xcms_is_modern(xcmsObj)
  }

  if (XCMSnExp_bool) {
    xfiles <- .xcms_files(xcmsObj)
    cond1 = (length(pa@fileList) > length(xfiles)) && (pa@f4f_link_type=='group')
    cond2 = !all(basename(pa@fileList)==basename(xfiles)) && (pa@f4f_link_type=='individual')
    cond3 = !all(names(pa@fileList)==basename(xfiles))
  }else{
    cond1 = (length(pa@fileList) > length(xcmsObj@filepaths)) && (pa@f4f_link_type=='group')
    cond2 = !all(basename(pa@fileList)==basename(xcmsObj@filepaths)) && (pa@f4f_link_type=='individual')
    cond3 = !all(names(pa@fileList)==basename(xcmsObj@filepaths))
  }

  if (cond1){
    # if more files in pa@filelist (can happen if some files were not processed with xcms because no MS1)
    # in this case we need to make sure any reference to a fileid is correct
    unevenFilelists = TRUE
  }else{
    unevenFilelists = FALSE
  }


  # if they are the same length, we check to make sure they are in the same order (only matters when
  # the f4f linking was for individual peaks)
  if(cond2){
    if(cond3){
      message('FILELISTS DO NOT MATCH')
      return(NULL)
    }else{
      if(XCMSnExp_bool){
        # XcmsExperiment keeps its file paths in its spectra.
        if (is(xcmsObj, 'XCMSnExp')) xcmsObj@processingData@files = unname(pa@fileList)
        xfiles <- unname(pa@fileList)
      }else{
        xcmsObj@filepaths <- unname(pa@fileList)
      }
    }
  }

  ###############################################
  # Add source
  ###############################################
  source <- data.frame(id=1,
                       name=paste('msPurity-database',  format(Sys.time(), "%Y-%m-%d-%I%M%S"), sep='-'),
                       parsing_software=paste('msPurity::createDatabase', packageVersion("msPurity")))

  add('source', 'id', NA, source)


  ###############################################
  # Add metab_compound (blank for time being)
  ###############################################
  metab_compound <- data.frame(inchikey_id=character(),
                               name=character(),
                               pubchem_id=character(),
                               chemspider_id=character(),
                               other_names=character(),
                               exact_mass=character(),
                               molecular_formula=character(),
                               molecular_weight=character(),
                               compound_class=character(),
                               smiles=character(),
                               created_at=character(),
                               updated_at=character())
  add('metab_compound', 'inchikey_id', NA, metab_compound, pk_type = 'TEXT')


  ###############################################
  # Add File info
  ###############################################
  nmsave <- names(pa@fileList) # this is for name tracking in Galaxy

  pa@fileList <- unname(pa@fileList)

  scaninfo <- purityTable(pa)
  fileList <- pa@fileList

  classInfo = .xcms_classes(xcmsObj)

  if(is.null(classInfo)){
    classInfo = rep(NA, length(fileList))
  }

  filedf <- data.frame(filename=basename(fileList),
                       filepth=fileList,
                       nm_save=nmsave,
                       fileid=seq(1, length(fileList)),
                       class=classInfo
  )

  add('fileinfo', 'fileid', NA, filedf)



  ###############################################
  # Add c_peaks (i.e. XCMS individual peaks)
  ###############################################

  if(XCMSnExp_bool){
    cPeaks <- xcms::chromPeaks(xcmsObj)
  }else{
    cPeaks <- xcmsObj@peaks
  }

  # Normally we expect the filelists to always be the same size, but there can be times when
  # MS/MS is collected without any full scan, or for some reasons it is not processed with xcms,
  # in these cases we need to ensure that fileids are correct
  if (unevenFilelists){
    if(XCMSnExp_bool){
      cPeaks[,'sample'] <- match(basename(xfiles[cPeaks[,'sample']]), filedf$filename)
    }else{
      cPeaks[,'sample'] <- match(basename(xcmsObj@filepaths[cPeaks[,'sample']]), filedf$filename)
    }
  }

  cPeaks <- data.frame(cbind('cid'=1:nrow(cPeaks), cPeaks))
  ccn <- colnames(cPeaks)

  colnames(cPeaks)[which(ccn=='sample')] <- 'fileid'
  colnames(cPeaks)[which(ccn=='into')] <- '_into'
  if ('i' %in% colnames(cPeaks)){
    cPeaks <- cPeaks[,-which(ccn=='i')]
  }
  fks_fileid <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'))
  add('c_peaks', 'cid', fks_fileid, cPeaks)


  ###############################################
  # Add c_peak_groups (i.e. XCMS grouped peaks) #chromPeaks(xset)
  ###############################################
  if (is.matrix(grpPeaklist)){
    grpPeaklist <- data.frame(grpPeaklist)
  }

  #check purpose
  colnames(grpPeaklist)[which(colnames(grpPeaklist)=='into')] <- '_into'

  if(XCMSnExp_bool){
    grpPeaklist$grp_name = .xcms_feature_names(xcmsObj)
    #convert list of peakIDX values to string - required for dbWriteTable, below
    grpPeaklist$peakidx = apply(grpPeaklist, 1, function(row){ paste(unlist(row['peakidx']), collapse = ', ')} )
  }else{
    grpPeaklist$grp_name <- xcms::groupnames(xcmsObj)
  }

  grpPeaklist <- grpPeaklist[order(grpPeaklist$grpid),]

  colnames(grpPeaklist)[colnames(grpPeaklist)=='rtmed'] = 'rt'
  colnames(grpPeaklist)[colnames(grpPeaklist)=='mzmed'] = 'mz'

  add('c_peak_groups', 'grpid', NA, grpPeaklist)

  ###############################################
  # Add MANY-to-MANY links for c_peak to c_peak_group
  ###############################################
  c_peak_X_c_peak_group <- getGroupPeakLink(xcmsObj = xcmsObj, method = 'medret', XCMSnExp_bool = XCMSnExp_bool)

  fks_for_cxg <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'),
                      'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks')
  )

  add('c_peak_X_c_peak_group', 'cXg_id', fks_for_cxg, c_peak_X_c_peak_group)

  ###############################################
  # Add s_peak_meta (i.e. scan information)
  ###############################################
  dropc <- c('filename', 'id')
  scaninfo <- scaninfo[,!colnames(scaninfo) %in% dropc]
  xx = c("name", "collision_energy", "ms_level", "accession", "resolution", "polarity",
         "fragmentation_type", "precursor_type", "instrument_type", "instrument",
         "copyright", "column", "mass_accuracy", "mass_error", "origin", "splash",
         "retention_index", "retention_time", "inchikey_id")
  scaninfo[xx] <- NA

  scaninfo$sourceid <- 1

  scaninfo$retention_time <- scaninfo$retentionTime
  scaninfo$precursor_mz <- scaninfo$precursorMZ
  scaninfo$spectrum_type <- 'scan'
  scaninfo <- update_cn_order(name_pk = 'pid',names_fk= 'fileid', df = scaninfo)



  ###############################################
  # Add s_peaks (all fragmentation spectra either scans or averaged)
  ###############################################
  # get all the fragmentation from the scans
    if (isTRUE(pa@filter_frag_params$allfrag)) {
        # These peak-row ids are local to assessed spectra; the scan route
        # retains the original all-MS>1 peak positions, including gaps.
        speaks <- allFragSpectra(pa, legacy = TRUE)
        if (!nrow(speaks)) {
            speaks <- .pa_empty_scan_peaks(scaninfo$pid)
            speaks$snr <- numeric()
            speaks$ra <- numeric()
            for (flag in c("purity_pass_flag", "intensity_pass_flag",
                    "ra_pass_flag", "snr_pass_flag", "pass_flag")) {
                speaks[[flag]] <- logical(nrow(speaks))
            }
            speaks <- speaks[, c("sid", "pid", "fileid", "scan", "mz", "i",
                    "snr", "ra", "type", "purity_pass_flag",
                    "intensity_pass_flag", "ra_pass_flag", "snr_pass_flag",
                    "pass_flag")]
        }
    } else if (scanPeaks == "count") {
        raw_scan_peaks <- .pa_scan_peak_count(pa, scaninfo$pid)
        speaks <- .pa_empty_scan_peaks(scaninfo$pid)
    } else {
        speaks <- getScanPeaks(pa)
    }
    speaks$grpid <- rep(NA_real_, nrow(speaks))

  av_list <- averagedSpectra(pa)
  if (length(av_list)>0){
    av_spectra <- plyr::ldply(av_list, getAvSpectraForGrp)
    colnames(av_spectra)[1] <- 'grpid'
    av_spectra$grpid <- names(av_list)[av_spectra$grpid]
    colnames(av_spectra)[colnames(av_spectra)=='sample'] <- 'fileid'

    colnames(av_spectra)[colnames(av_spectra)=='method'] = 'type'
    # The scan peak-row ids come first, whether or not they are kept.
    last_sid <- if (is.null(raw_scan_peaks)) max(c(0, speaks$sid))
                else attr(raw_scan_peaks, "last_sid")
    av_spectra$sid <- last_sid + seq_len(nrow(av_spectra))

    prePid  <- paste(av_spectra$grpid, av_spectra$fileid, av_spectra$type)
    newPids <- max(c(0, scaninfo$pid)) +
        as.numeric(factor(prePid, levels = unique(prePid)))

    # add new scaninfo details (with just pid for now)
    av_spectra$pid <- newPids

    topnav <- plyr::ddply(av_spectra, ~pid, getXcmsGrpDetails, grpPeaklist)

    grpidx <- which(grpPeaklist$grpid %in% topnav$grpid)
    if (is.null(topnav$fileid)){
      topnvfileids <- NA
    }else{
      topnvfileids <- topnav$fileid
    }

    scaninfo <- plyr::rbind.fill(scaninfo, data.frame(pid=topnav$pid,
                                                      fileid=topnvfileids,
                                                      spectrum_type=topnav$type,
                                                      precursor_mz=topnav$precusor_mz,
                                                      retention_time=topnav$retention_time,
                                                      inPurity=topnav$inPurity,
                                                      grpid=topnav$grpid))

    speaks <- merge(speaks, av_spectra, all = TRUE)

    speaks[is.na(speaks)] <- NA

  }


  if (!anyNA(metadata)){
    if(!is.null(metadata$polarity)){
      scaninfo$polarity <- metadata$polarity
    }

    if(!is.null(metadata$instrument_type)){
      scaninfo$instrument_type <- metadata$instrument_type
    }

    if(!is.null(metadata$instrument)){
      scaninfo$instrument <- metadata$instrument
    }
  }


  fks4speaks <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'),
                     'pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'))

  fks4smeta <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'))
  add('s_peak_meta', 'pid', fks4smeta, scaninfo)
  add('s_peaks', 'sid', fks4speaks, speaks)



  if (pa@f4f_link_type=='individual'){
    ###############################################
    # Add MANY-to-MANY links for c_peak to s_peak_meta
    ###############################################
    grpdf <- pa@grped_df
    c_peak_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'cid')])
    c_peak_X_s_peak_meta <- cbind('cXp_id'=1:nrow(c_peak_X_s_peak_meta), c_peak_X_s_peak_meta)

    fks_for_cXs <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'),
                        'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks'))

    add('c_peak_X_s_peak_meta', 'cXp_id', fks_for_cXs, c_peak_X_s_peak_meta)
  }else{
    ###############################################
    # Add MANY-to-MANY links for c_peak_group to s_peak_meta
    ###############################################
    grpdf <- pa@grped_df
    c_peak_group_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'grpid')])
    c_peak_group_X_s_peak_meta <- cbind('gXp_id'=1:nrow(c_peak_group_X_s_peak_meta), c_peak_group_X_s_peak_meta)

    fks_for_cXs <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'),
                        'grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'))

    add('c_peak_group_X_s_peak_meta', 'gXp_id', fks_for_cXs, c_peak_group_X_s_peak_meta)


  }



  if (!is.null(xsa)){
    ###############################################
    # Add CAMERA ruleset
    ###############################################
    if(is.null(xsa@ruleset)){
      rules_pos <- utils::read.table(system.file(file.path('rules', 'extended_adducts_pos.csv') , package = "CAMERA"), header = TRUE)
      rules_neg <- utils::read.csv(system.file(file.path('rules', 'extended_adducts_neg.csv') , package = "CAMERA"))
      rules <- rbind(rules_pos, rules_neg)

    }else{
      rules <- xsa@ruleset
    }
    rules$rule_id <- 1:nrow(rules)
    add('adduct_rules', 'rule_id', NA, rules)

    ###############################################
    # Add neutral mass groups
    ###############################################
    annoGrp <- data.frame(xsa@annoGrp)
    colnames(annoGrp)[1] <- 'nm_id'
    add('neutral_masses', 'nm_id', NA, annoGrp)

    ###############################################
    # Add adduct annotations
    ###############################################
    annoID <- data.frame(xsa@annoID)
    colnames(annoID) <- c('grpid', 'nm_id', 'rule_id', 'parentID')
    annoID  <- cbind('add_id'=1:nrow(annoID), annoID)

    fks_adduct <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_group'),
                       'nm_id'=list('new_name'='nm_id', 'ref_name'='nm_id', 'ref_table'='neutral_masses'),
                       'rule_id'=list('new_name'='rule_id', 'ref_name'='rule_id', 'ref_table'='adduct_rules')
    )

    add('adduct_annotations', 'add_id', fks_adduct, annoID)

    ###############################################
    # Add isotope annotations
    ###############################################
    isoID <- data.frame(xsa@isoID)
    colnames(isoID) <- c('c_peak_group1_id', 'c_peak_group2_id', 'iso', 'charge')
    isoID <- cbind('iso_id'=1:nrow(isoID), isoID)

    fk_isotope <- list('c_peak_group1_id'=list('new_name'='c_peak_group1_id',
                                               'ref_name'='grpid',
                                               'ref_table'='c_peak_group'),

                       'c_peak_group2_id'=list('new_name'='c_peak_group2_id',
                                               'ref_name'='grpid',
                                               'ref_table'='c_peak_group')
    )

    add('isotope_annotations', 'iso_id', fk_isotope, isoID)
  }

  list(tables = tables, raw_scan_peaks = raw_scan_peaks)
}


# How many raw peaks the scans `pids` have, and the last of their peak-row
# ids (`last_sid` attribute), numbered as getScanPeaks() numbers them: a
# running count over every MS>1 peak of the files, in file order. Taken from
# the source peak positions recorded at purityA(); objects without them are
# read from their files.
.pa_scan_peak_count <- function(pa, pids) {
    sp <- pa@spectra
    vars <- Spectra::spectraVariables(sp)
    k <- match(pids, sp$pid)
    if (length(pids) && !anyNA(k) &&
            all(c("msPurityRawPeakOffset", "msPurityRawPeakCount",
                  "msPurityRawPeakIdBasis") %in% vars) &&
            all(sp$msPurityRawPeakIdBasis[k] %in% "source_reference")) {
        offset <- sp$msPurityRawPeakOffset[k]
        count <- sp$msPurityRawPeakCount[k]
        if (!anyNA(offset) && !anyNA(count)) {
            n <- sum(count)
            attr(n, "last_sid") <- if (n) max((offset + count)[count > 0])
                                   else 0
            return(n)
        }
    }
    peaks <- getScanPeaks(pa)
    n <- nrow(peaks)
    attr(n, "last_sid") <- max(c(0, peaks$sid))
    n
}

getAvSpectraForGrp <- function(x){

  if (length(x$av_intra)>0){
    av_intra_df <- plyr::ldply(x$av_intra, .id = 'fileid')

    if (nrow(av_intra_df)==0){
      av_intra_df <- NULL
    }else{
      av_intra_df$method <- 'intra'
    }

  }else{
    av_intra_df <- NULL
  }

  if ((is.null(x$av_inter)) || (nrow(x$av_inter)==0)){
    av_inter_df <- NULL
  }else{
    av_inter_df <- x$av_inter
    av_inter_df$method <- 'inter'
  }

  if ((is.null(x$av_all)) || (nrow(x$av_all)==0)){
    av_all_df <- NULL
  }else{
    av_all_df <- x$av_all
    av_all_df$method <- 'all'
  }

  combined <- plyr::rbind.fill(av_intra_df, av_inter_df, av_all_df)

  return(combined)

}

getXcmsGrpDetails <- function(x, grpPeaklist){
  x <- x[1,]

  grpPeaklisti <- grpPeaklist[grpPeaklist$grpid==x$grpid,]

  x$retention_time <- grpPeaklisti$rt
  x$precusor_mz <- grpPeaklisti$mz
  return(x)
}

real_or_rest <- function(x){
  if(is.numeric(x)){
    return('REAL')
  }else{
    return('TEXT')
  }
}

get_column_info <- function(x, data_type){return(paste(x, data_type[x], sep = ' '))}

get_create_query <- function(pk, fks=NA, table_name, df, pk_type='INTEGER'){

  cns <- colnames(df)

  if (anyNA(fks)){
    cns_sml <- cns[which(!cns %in% pk)]
  }else{
    cns_sml <- cns[which(!cns %in% c(pk, names(fks)))]
  }

  data_type <- lapply(df[1, cns_sml], real_or_rest)

  colmninfo <- sapply(cns_sml, get_column_info, data_type=data_type)

  columninfo <- paste(colmninfo, collapse = ', ')

  pkinfo <- paste(pk, sprintf(' %s NOT NULL PRIMARY KEY', pk_type), sep='')
  if (anyNA(fks)){

    if (columninfo==''){
      allcolinfo <- pkinfo
    }else{
      allcolinfo <- paste(c( pkinfo, columninfo), collapse=', ')
    }

  }else{
    fks_info1 <- sapply(fks, function(x){
      paste(x$new_name, 'INTEGER')
    })

    fks_info2 <- sapply(fks, function(x){
      paste('FOREIGN KEY (', x$new_name, ') REFERENCES', x$ref_table, '(', x$ref_name, ')', sep=' ')
    })

    fksinfo <- paste(c(fks_info1, fks_info2), collapse = ', ')

    if (columninfo==''){
      allcolinfo <- paste(c(pkinfo, fksinfo), collapse=', ')
    }else{
      allcolinfo <- paste(c(pkinfo, columninfo,  fksinfo), collapse=', ')

    }


  }

  return(paste('CREATE TABLE', table_name, '(', allcolinfo, ')', sep=' '))

}

update_cn_order <- function(name_pk, names_fk, df){
  # primary key needs to be at the start
  # foreign keys at the end
  cn <- colnames(df)

  if (anyNA(names_fk)){
    columnorder <- c(name_pk, cn[!cn %in% name_pk])
  }else{
    columnorder <- c(name_pk, cn[!cn %in% c(name_pk, names_fk)], names_fk)
  }
  return(df[,columnorder])
}



# A table as written: primary key first, foreign keys last, and column names
# without "." or "-".
.db_frame <- function(df, name_pk, fks) {
  names_fk <- if (anyNA(fks)) NA else names(fks)
  df <- update_cn_order(name_pk=name_pk, names_fk=names_fk, df = df)
  names(df) <- gsub( ".",  "_", names(df), fixed = TRUE)
  names(df) <- gsub( "-",  "_", names(df), fixed = TRUE)
  df
}

custom_dbWriteTable <- function(name_pk, fks, df, table_name, con, pk_type='INTEGER'){
  df <- .db_frame(df, name_pk, fks)

  query <- get_create_query(pk=name_pk, fks=fks, table_name=table_name, df=df, pk_type=pk_type)

  sqr <- DBI::dbSendQuery(con, query)
  DBI::dbClearResult(sqr)

  head(df)

  DBI::dbWriteTable(con, name=table_name, value=df, row.names=FALSE, append=TRUE)


}



getGroupPeakLink <- function(xcmsObj, method='medret', XCMSnExp_bool){

  if(XCMSnExp_bool){
    gidx <- xcms::featureDefinitions(xcmsObj)$peakidx
    bestpeaks <- xcms::featureValues(xcmsObj, method = method, value = 'index')
    sids = xcms::chromPeaks(xcmsObj)[,'sample']
    peaks_df = data.frame(xcms::chromPeaks(xcmsObj))
  }else{
    gidx <- xcmsObj@groupidx
    bestpeaks <- xcms::groupval(xcmsObj, method=method, value = 'index')
    sids = xcmsObj@peaks[,'sample']
    filenames = rownames(xcmsObj@phenoData)
    peaks_df = data.frame(xcmsObj@peaks)
  }

  idis <- unlist(plyr::dlply(peaks_df, ~sample, function(x){ 1:nrow(x)}))

  #for(i in 1:1000){
  for(i in 1:length(gidx)){
    gidxi <- gidx[[i]]
    bpi <- bestpeaks[i,]

    grpid <- rep(i, length(gidxi))
    #sid=sids[gidxi]

    im <- cbind('grpid'=grpid, 'cid'=gidxi, 'idi'=idis[gidxi], 'bestpeak'=((gidxi %in% bpi) * 1))
    #im <- data.frame(im)
    # multipling by 1 converts TRUE FALSE to 1 0

    if(i==1){
      allm <- im
    }else{
      allm <- rbind(allm, im)
    }
  }

  rownames(allm) <- 1:nrow(allm)
  allm <- data.frame(allm)
  allm$cXg_id <- 1:nrow(allm)
  return(allm)

}

# The SQLite schema holds MS2 scans, feature links and spectra only.
.pa_require_ms2_export <- function(pa) {
    if (any(.pa_av_configured(pa@avSpectra)))
        stop("Database export of averages computed with ",
            "AverageFragSpectraParam is not supported. Retain the processed ",
            "Spectra object, or average without param before export.")
    for (field in c("fragSpectra", "avSpectra")) {
        sp <- methods::slot(pa, field)
        if (length(sp) && any(Spectra::msLevel(sp) > 2L, na.rm = TRUE))
            stop("Processed MS3/MS4 database export is not supported by ",
                "the SQLite writer. Use createDatabase(format = \"parquet\") ",
                "or export MS2-only results.")
    }
    g <- pa@grped_df
    if (nrow(g) && "msLevel" %in% names(g) && any(g$msLevel > 2L, na.rm = TRUE))
        stop("MS3/MS4 feature links are not supported by the SQLite writer. ",
            "Use createDatabase(format = \"parquet\") or link MS2 only.")
    invisible(TRUE)
}

# The Parquet writer holds MS2-MS4 scans, links and filtered spectra; averaged
# spectra above MS2 are not produced yet.
.pa_require_msn_export <- function(pa) {
    sp <- pa@avSpectra
    if (length(sp) && any(Spectra::msLevel(sp) > 2L, na.rm = TRUE))
        stop("Averaged MS3/MS4 spectra cannot be exported yet.")
    invisible(TRUE)
}
