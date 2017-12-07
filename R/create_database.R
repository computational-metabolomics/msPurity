#' @title Create database
#'
#' @description
#' Create and SQLite database of an LC-MS(/MS) experiment
#'
#' @param pa purityA object; Needs to be the same used for frag4feature function
#' @param xset xcms object; Needs to be the same used for frag4feature function
#' @param db_name character [optional]; Name of the result database
#' @param grp_peaklist dataframe [optional]; Can use any peak dataframe. Still needs to be derived from the xset object though
#' @param out_dir character; Out directory for the SQLite result database
#' @return path to SQLite database and database name
#' @export
create_database <-  function(pa, xset, out_dir, grp_peaklist=NA, db_name=NA){
  ########################################################
  # Export the target data into sqlite database
  ########################################################
  if (is.na(db_name)){
    db_name <- paste('lcmsms_data', format(Sys.time(), "%Y-%m-%d-%I%M%S"), '.sqlite', sep="-")
  }

  if (!is.data.frame(grp_peaklist)){
    grp_peaklist <- xcms::peakTable(xset)
    grp_peaklist <- data.frame(cbind('grpid'=1:nrow(grp_peaklist), grp_peaklist))
  }

  message("Creating a database of fragmentation spectra and LC features")
  target_db_pth <- export_2_sqlite(pa, grp_peaklist, xset, out_dir, db_name)

  return(target_db_pth)

}


export_2_sqlite <- function(pa, grp_peaklist, xset, out_dir, db_name){

  db_pth <- file.path(out_dir, db_name)

  con <- DBI::dbConnect(RSQLite::SQLite(),db_pth)

  ###############################################
  # Add File info
  ###############################################
  scan_info <- pa@puritydf
  filepth_df <- data.frame(cbind('filename'=basename(pa@fileList), 'filepth'=pa@fileList))
  filedf <- unique(scan_info[ ,c('fileid', 'filename')])
  filedf <- merge(filedf, filepth_df)

  filedf <- filedf[,c('fileid', 'filename', 'filepth')]

  custom_dbWriteTable(name_pk = 'fileid', fks=NA, table_name = 'fileinfo', df=filedf, con=con)

  ###############################################
  # Add c_peaks (i.e. XCMS individual peaks)
  ###############################################

  c_peaks <- xset@peaks

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

  custom_dbWriteTable(name_pk = 'grpid', fk=NA, table_name = 'c_peak_groups', df=grp_peaklist, con=con)

  ###############################################
  # Add s_peak_meta (i.e. scan information)
  ###############################################
  dropc <- c('filename')
  scan_info <- scan_info[,!colnames(scan_info) %in% dropc]
  scan_info <- update_cn_order(name_pk = 'pid',names_fk= 'fileid', df = scan_info)
  custom_dbWriteTable(name_pk = 'pid', fk=fks_fileid, table_name = 's_peak_meta', df=scan_info, con=con)

  ###############################################
  # Add s_peaks (i.e. the mz, i from each scan)
  ###############################################
  # ensure the filedf is in the same order as the scan_info file
  # other wise the s_peak_meta ids might not match up
  filedf <- filedf[as.numeric(unique(scan_info$fileid)),]

  scanpeaks_frag <- plyr::ddply(filedf, ~ fileid, scan_peaks_4_db)

  comb <- paste(scanpeaks_frag[,1], scanpeaks_frag[,2], sep=' ')
  scanpeaks_frag <- cbind(1:nrow(scanpeaks_frag), cumsum(!duplicated(comb)), scanpeaks_frag)
  colnames(scanpeaks_frag) <- c('sid','pid', 'fileid', 'scan', 'mz', 'i')
  fks_pid <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'))

  custom_dbWriteTable(name_pk = 'sid', fk=append(fks_fileid, fks_pid),
                      table_name = 's_peaks', df=scanpeaks_frag, con=con)


  ###############################################
  # Add MANY-to-MANY links for c_peak to c_peak_group
  ###############################################
  grpdf <- pa@grped_df
  c_peak_X_c_peak_group <- unique(grpdf[ ,c('grpid', 'cid')])
  c_peak_X_c_peak_group <- cbind('cXg_id'=1:nrow(c_peak_X_c_peak_group), c_peak_X_c_peak_group)

  fks_grpid <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'))
  fks_cid <- list('cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks'))

  custom_dbWriteTable(name_pk = 'cXg_id', fk=append(fks_grpid, fks_cid),
                      table_name ='c_peak_X_c_peak_group', df=c_peak_X_c_peak_group, con=con)

  ###############################################
  # Add MANY-to-MANY links for c_peak to s_peak_meta
  ###############################################
  c_peak_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'cid')])
  c_peak_X_s_peak_meta <- cbind('cXp_id'=1:nrow(c_peak_X_s_peak_meta), c_peak_X_s_peak_meta)

  fks_pid <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'))

  custom_dbWriteTable(name_pk = 'cXp_id', fk=append(fks_pid, fks_cid),
                      table_name ='c_peak_X_s_peak_meta', df=c_peak_X_s_peak_meta, con=con)

  DBI::dbDisconnect(con)
  return(db_pth)

}

real_or_rest <- function(x){
  if(is.numeric(x)){
    return('REAL')
  }else{
    return('TEXT')
  }
}

get_column_info <- function(x, data_type){return(paste(x, data_type[x], sep = ' '))}

get_create_query <- function(pk, fks=NA, table_name, df){

  cns <- colnames(df)

  if (anyNA(fks)){
    cns_sml <- cns[which(!cns %in% pk)]
  }else{
    cns_sml <- cns[which(!cns %in% c(pk, names(fks)))]
  }

  data_type <- lapply(df[1, cns_sml], real_or_rest)

  colmninfo <- sapply(cns_sml, get_column_info, data_type=data_type)

  columninfo <- paste(colmninfo, collapse = ', ')

  pkinfo <- paste(pk, ' INTEGER NOT NULL PRIMARY KEY', sep='')
  if (anyNA(fks)){

    if (columninfo==''){
      allcolinfo <- pkinfo
    }else{
      allcolinfo <- paste(c( pkinfo, columninfo), collapse=', ')
    }

  }else{
    fks_info1 <- sapply(fks, function(x){
      paste(x$ref_name, 'INTEGER')
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
  scanpeaks_df <- plyr::ldply(scanpeaks[scans$seqNum[scans$msLevel>1]], .id=T)
}


custom_dbWriteTable <- function(name_pk, fks, df, table_name, con){
  if (anyNA(fks)){
    names_fk =  NA
  }else{
    names_fk =names(fks)
  }
  df <- update_cn_order(name_pk=name_pk, names_fk=names_fk, df = df)
  query <- get_create_query(pk=name_pk, fks=fks, table_name=table_name, df=df)
  sqr <- DBI::dbSendQuery(con, query)
  DBI::dbClearResult(sqr)
  DBI::dbWriteTable(con, name=table_name, value=df, row.names=F, append=T)


}
