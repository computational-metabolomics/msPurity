#' @title Create database [deprecated]
#'
#' @description
#' Create and SQLite database of an LC-MS(/MS) experiment
#'
#' msPurity::create_database is deprecated. Please use msPurity::createDatabase for future use
#'
#' @param pa purityA object; Needs to be the same used for frag4feature function
#' @param xset xcms object; Needs to be the same used for frag4feature function (this will be ignored when using xsa parameter)
#' @param xsa CAMERA object \[optional\]; if CAMERA object is used, we ignore the xset parameter input and obtain all information
#'                          from the xset object nested with the CAMERA xsa object. Adduct and isotope information
#'                          will be included into the database when using this parameter. The underlying xset object must
#'                          be the one used for the frag4feature function
#' @param db_name character \[optional\]; Name of the result database
#' @param grp_peaklist dataframe \[optional\]; Can use any peak dataframe. Still needs to be derived from the xset object though
#' @param out_dir character; Out directory for the SQLite result database
#' @return path to SQLite database and database name
#'
#' @examples
#'
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#' #            package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #xset <- xcms::xcmsSet(msmsPths)
#' #xset <- xcms::group(xset)
#'
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xset)
#' #pa <- averageAllFragSpectra(pa)
#' #db_pth <- create_database(pa, xset)
#'
#' # Run from previously generated data
#' pa <- readRDS(system.file("extdata", "tests", "purityA",
#'                           "9_averageAllFragSpectra_with_filter_pa.rds",
#'                           package="msPurity"))
#' xset <- readRDS(system.file("extdata","tests", "xcms",
#'                             "msms_only_xset.rds", package="msPurity"))
#'
#' # Need to ensure the filelists are matching
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML",
#'                                    package="msPurityData"),
#'                                    full.names = TRUE, pattern = "MSMS")
#' pa@fileList[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
#' pa@fileList[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
#' xset@filepaths[1] <- msmsPths[basename(msmsPths)=="LCMSMS_1.mzML"]
#' xset@filepaths[2] <- msmsPths[basename(msmsPths)=="LCMSMS_2.mzML"]
#' db_pth <- create_database(pa, xset)
#' @export
create_database <-  function(pa, xset, xsa=NULL, out_dir='.', grp_peaklist=NA, db_name=NA){
  ########################################################
  # Export the target data into sqlite database
  ########################################################
  if(!is.null(xset)){
    xset <- getxcmsSetObject(xset)
  }

  if (is.na(db_name)){
    db_name <- paste('lcmsms_data', format(Sys.time(), "%Y-%m-%d-%I%M%S"), '.sqlite', sep="-")
  }

  if (!is.data.frame(grp_peaklist)){
    if (is.null(xsa)){
      grp_peaklist <- xcms::peakTable(xset)
    }else{
      grp_peaklist <- CAMERA::getPeaklist(xsa)
    }

    grp_peaklist <- data.frame(cbind('grpid'=1:nrow(grp_peaklist), grp_peaklist))
  }

  message("Creating a database of fragmentation spectra and LC features")
  target_db_pth <- export_2_sqlite(pa, grp_peaklist, xset, xsa, out_dir, db_name)

  return(target_db_pth)

}


export_2_sqlite <- function(pa, grp_peaklist, xset, xsa, out_dir, db_name){

  if(!is.null(xsa)){
    # if user has supplied camera object we use the xset that the camera object
    # is derived from
    xset <- xsa@xcmsSet

  }


  if ((length(pa@fileList) > length(xset@filepaths)) && (pa@f4f_link_type=='group')){
    # if more files in pa@filelist (can happen if some files were not processed with xcms because no MS1)
    # in this case we need to make sure any reference to a fileid is correct
    uneven_filelists = TRUE
  }else{
    uneven_filelists = FALSE
  }


  # if they are the same length, we check to make sure they are in the same order (only matters when
  # the f4f linking was for individual peaks)
  if(!all(basename(pa@fileList)==basename(xset@filepaths)) && (pa@f4f_link_type=='individual')){
      if(!all(names(pa@fileList)==basename(xset@filepaths))){
        message('FILELISTS DO NOT MATCH')
        return(NULL)
      }else{
        xset@filepaths <- unname(pa@fileList)
      }
    }



  db_pth <- file.path(out_dir, db_name)

  con <- DBI::dbConnect(RSQLite::SQLite(),db_pth)

  ###############################################
  # Add File info
  ###############################################
  nm_save <- names(pa@fileList) # this is for name tracking in Galaxy
  pa@fileList <- unname(pa@fileList)

  scan_info <- pa@puritydf
  fileList <- pa@fileList


  filedf <- data.frame(filename=basename(fileList),
                       filepth=fileList,
                       nm_save=nm_save,
                       fileid=seq(1, length(fileList)),
                       class=xset@phenoData$class
                       )

  custom_dbWriteTable(name_pk = 'fileid', fks=NA, table_name = 'fileinfo', df=filedf, con=con)



  ###############################################
  # Add c_peaks (i.e. XCMS individual peaks)
  ###############################################

  c_peaks <- xset@peaks

  # Normally we expect the filelists to always be the same size, but there can be times when
  # MS/MS is collected without any full scan, or for some reasons it is not processed with xcms,
  # in these cases we need to ensure that fileids are correct
  if (uneven_filelists){
    c_peaks[,'sample'] <- match(basename(xset@filepaths[c_peaks[,'sample']]), filedf$filename)
  }

  c_peaks <- data.frame(cbind('cid'=1:nrow(c_peaks), c_peaks))
  ccn <- colnames(c_peaks)

  colnames(c_peaks)[which(ccn=='sample')] <- 'fileid'
  colnames(c_peaks)[which(ccn=='into')] <- '_into'
  if ('i' %in% colnames(c_peaks)){
    c_peaks <- c_peaks[,-which(ccn=='i')]
  }
  fks_fileid <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'))
  custom_dbWriteTable(name_pk = 'cid', fks=fks_fileid, table_name = 'c_peaks', df=c_peaks, con=con)


  ###############################################
  # Add c_peak_groups (i.e. XCMS grouped peaks)
  ###############################################
  if (is.matrix(grp_peaklist)){
    grp_peaklist <- data.frame(grp_peaklist)
  }
  colnames(grp_peaklist)[which(colnames(grp_peaklist)=='into')] <- '_into'
  grp_peaklist$grp_name <- xcms::groupnames(xset)
  custom_dbWriteTable(name_pk = 'grpid', fks=NA, table_name = 'c_peak_groups', df=grp_peaklist, con=con)

  ###############################################
  # Add s_peak_meta (i.e. scan information)
  ###############################################
  dropc <- c('filename')
  scan_info <- scan_info[,!colnames(scan_info) %in% dropc]
  scan_info <- update_cn_order(name_pk = 'pid',names_fk= 'fileid', df = scan_info)
  custom_dbWriteTable(name_pk = 'pid', fks=fks_fileid, table_name = 's_peak_meta', df=scan_info, con=con)

  ###############################################
  # Add s_peaks (i.e. the mz, i from each scan)
  ###############################################
  # ensure the filedf is in the same order as the scan_info file
  # other wise the s_peak_meta ids might not match up (should be fixed from previous matching)
  # filedf <- filedf[as.numeric(unique(scan_info$fileid)),]

  scanpeaks_frag <- plyr::ddply(filedf, ~ fileid, scan_peaks_4_db)

  comb <- paste(scanpeaks_frag[,1], scanpeaks_frag[,2], sep=' ')
  scanpeaks_frag <- cbind(1:nrow(scanpeaks_frag), cumsum(!duplicated(comb)), scanpeaks_frag)
  colnames(scanpeaks_frag) <- c('sid','pid', 'fileid', 'scan', 'mz', 'i')
  fks_pid <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'))

  custom_dbWriteTable(name_pk = 'sid', fks=append(fks_fileid, fks_pid),
                      table_name = 's_peaks', df=scanpeaks_frag, con=con)


  ###############################################
  # Add MANY-to-MANY links for c_peak to c_peak_group
  ###############################################
  c_peak_X_c_peak_group <- get_group_peak_link(xset)

  fks_for_cxg <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'),
                      'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks')
                      )

  custom_dbWriteTable(name_pk = 'cXg_id', fks=fks_for_cxg,
                      table_name ='c_peak_X_c_peak_group', df=c_peak_X_c_peak_group, con=con)

  if (pa@f4f_link_type=='individual'){
    ###############################################
    # Add MANY-to-MANY links for c_peak to s_peak_meta
    ###############################################
    grpdf <- pa@grped_df
    c_peak_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'cid')])
    c_peak_X_s_peak_meta <- cbind('cXp_id'=1:nrow(c_peak_X_s_peak_meta), c_peak_X_s_peak_meta)

    fks_for_cXs <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'),
                        'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks'))

    custom_dbWriteTable(name_pk = 'cXp_id', fks=fks_for_cXs,
                        table_name ='c_peak_X_s_peak_meta', df=c_peak_X_s_peak_meta, con=con)
  }else{
    ###############################################
    # Add MANY-to-MANY links for c_peak_group to s_peak_meta
    ###############################################
    grpdf <- pa@grped_df
    c_peak_group_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'grpid')])
    c_peak_group_X_s_peak_meta <- cbind('gXp_id'=1:nrow(c_peak_group_X_s_peak_meta), c_peak_group_X_s_peak_meta)

    fks_for_cXs <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'),
                        'grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'))

    custom_dbWriteTable(name_pk = 'gXp_id', fks=fks_for_cXs,
                        table_name ='c_peak_group_X_s_peak_meta', df=c_peak_group_X_s_peak_meta, con=con)


  }

  if (length(pa@av_spectra)>0){

    av_spectra <- plyr::ldply(pa@av_spectra, get_av_spectra_for_db)

    if (nrow(av_spectra)==0){
      message('No average spectra to use for database')
    }else{
      # for some reason the names are not being saved for the list as a column, so we just get them back
      colnames(av_spectra)[1] <- 'grpid'
      av_spectra$grpid <- names(pa@av_spectra)[av_spectra$grpid]

      colnames(av_spectra)[2] <- 'fileid'
      av_spectra$avid <- 1:nrow(av_spectra)
      fks_for_av_spectra <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'),
                                 'grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups')
      )

      custom_dbWriteTable(name_pk = 'avid', fks=fks_for_av_spectra,
                          table_name ='av_peaks', df=av_spectra, con=con)
    }



  }


  if (!is.null(xsa)){
    ###############################################
    # Add CAMERA ruleset
    ###############################################
    if(is.null(xsa@ruleset)){
      rules_pos <- utils::read.table(system.file(file.path('rules', 'extended_adducts_pos.csv') , package = "CAMERA"), header = T)
      rules_neg <- utils::read.csv(system.file(file.path('rules', 'extended_adducts_neg.csv') , package = "CAMERA"))
      rules <- rbind(rules_pos, rules_neg)

    }else{
      rules <- xsa@ruleset
    }
    rules$rule_id <- 1:nrow(rules)
    custom_dbWriteTable(name_pk = 'rule_id', fks=NA,
                        table_name ='adduct_rules', df=rules, con=con)

    ###############################################
    # Add neutral mass groups
    ###############################################
    annoGrp <- data.frame(xsa@annoGrp)
    colnames(annoGrp)[1] <- 'nm_id'
    custom_dbWriteTable(name_pk = 'nm_id', fks=NA,
                        table_name ='neutral_masses', df=annoGrp, con=con)

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

    custom_dbWriteTable(name_pk = 'add_id', fks=fks_adduct,
                        table_name ='adduct_annotations', df=annoID, con=con)

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

    custom_dbWriteTable(name_pk = 'iso_id', fks=fk_isotope,
                        table_name ='isotope_annotations', df=isoID, con=con)
  }


  DBI::dbDisconnect(con)
  return(db_pth)

}


get_av_spectra_for_db <- function(x){

  if (length(x$av_intra)>0){
    av_intra_df <- plyr::ldply(x$av_intra)

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



scan_peaks_4_db <- function(x){

  mr <- mzR::openMSfile(as.character(x$filepth))
  scanpeaks <- mzR::peaks(mr)
  scans <- mzR::header(mr)
  names(scanpeaks) <- seq(1, length(scanpeaks))

  scanpeaks_df <- plyr::ldply(scanpeaks[scans$seqNum[scans$msLevel>1]], .id=TRUE)
}


custom_dbWriteTable <- function(name_pk, fks, df, table_name, con, pk_type='INTEGER'){
  if (anyNA(fks)){
    names_fk =  NA
  }else{
    names_fk =names(fks)
  }
  df <- update_cn_order(name_pk=name_pk, names_fk=names_fk, df = df)

  names(df) <- gsub( ".",  "_", names(df), fixed = TRUE)
  names(df) <- gsub( "-",  "_", names(df), fixed = TRUE)

  query <- get_create_query(pk=name_pk, fks=fks, table_name=table_name, df=df, pk_type=pk_type)

  sqr <- DBI::dbSendQuery(con, query)
  DBI::dbClearResult(sqr)

  head(df)

  DBI::dbWriteTable(con, name=table_name, value=df, row.names=FALSE, append=TRUE)


}



get_group_peak_link <- function(xset, method='medret'){

  gidx <- xset@groupidx
  bestpeaks <- xcms::groupval(xset, method=method)

  sids = xset@peaks[,'sample']
  filenames = rownames(xset@phenoData)

  idis <- unlist(plyr::dlply(data.frame(xset@peaks), ~sample, function(x){ 1:nrow(x)}))

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
