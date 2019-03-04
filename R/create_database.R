#' @title Create database
#'
#' @description
#' Create and SQLite database of an LC-MS(/MS) experiment
#'
#' @param pa purityA object; Needs to be the same used for frag4feature function
#' @param xset xcms object; Needs to be the same used for frag4feature function (this will be ignored when using xsa parameter)
#' @param xsa CAMERA object [optional]; if CAMERA object is used, we ignore the xset parameter input and obtain all information
#'                          from the xset object nested with the CAMERA xsa object. Adduct and isotope information
#'                          will be included into the database when using this parameter. The underlying xset object must
#'                          be the one used for the frag4feature function
#' @param db_name character [optional]; Name of the result database
#' @param grp_peaklist dataframe [optional]; Can use any peak dataframe. Still needs to be derived from the xset object though
#' @param out_dir character; Out directory for the SQLite result database
#' @return path to SQLite database and database name
#'
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
#' db_pth <- create_database(pa, xset)
#' @export
create_database <-  function(pa, xset, xsa=NULL, out_dir='.', grp_peaklist=NA, db_name=NA, grped_df=NULL, fileMS1=NULL, fileMS2=NULL){
  ########################################################
  # Export the target data into sqlite database
  ########################################################
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
  target_db_pth <- export_2_sqlite(pa, grp_peaklist, xset, xsa, out_dir, db_name,grped_df, fileMS1, fileMS2)

  return(target_db_pth)

}


export_2_sqlite <- function(pa, grp_peaklist, xset, xsa, out_dir, db_name, grped_df, fileMS1, fileMS2){

  if(!is.null(xsa)){
    # if user has supplied camera object we use the xset that the camera object
    # is derived from
    xset <- xsa@xcmsSet

  }


  if(pa@f4f_link_type=='group'){  
    # if more files in pa@filelist (can happen if some files were not processed with xcms because no MS1)
    # in this case we need to make sure any reference to a fileid is correct
    uneven_filelists = TRUE
  }else{
    uneven_filelists = FALSE
  }



  db_pth <- file.path(out_dir, db_name)

  con <- DBI::dbConnect(RSQLite::SQLite(),db_pth)

  ###############################################
  # Add File info
  ###############################################
  nm_save <- names(pa@fileList[fileMS2]) # this is for name tracking in Galaxy

  pa@fileList[fileMS2] <- unname(pa@fileList[fileMS2])
  MS1name <- paste(gsub("(^.+)\\..*$","\\1",basename(xset@filepaths[fileMS1])))
  MS2name <- paste(gsub("(^.+)\\..*$","\\1",basename(nm_save)))

  scan_info <- pa@puritydf[which(pa@puritydf[,"fileid"] == fileMS2),]
  fileList <- pa@fileList[fileMS2]

  filedf <- data.frame(filename=basename(fileList),
                       filepth=fileList,
                       nm_save=nm_save,
                       fileid=seq(1, length(fileList)),
                       class=as.character(xset@phenoData[which(xset@phenoData[,"sample_name"]==MS2name),"class"])
                       )

  custom_dbWriteTable(name_pk = 'fileid', fks=NA, table_name = paste0('fileinfo_',MS1name,"_",MS2name), df=filedf, con=con)



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
  custom_dbWriteTable(name_pk = 'cid', fks=fks_fileid, table_name = paste0('c_peaks_',MS1name,"_",MS2name), df=c_peaks, con=con)


  ###############################################
  # Add c_peak_groups (i.e. XCMS grouped peaks)
  ###############################################
  if (is.matrix(grp_peaklist)){
    grp_peaklist <- data.frame(grp_peaklist)
  }
  colnames(grp_peaklist)[which(colnames(grp_peaklist)=='into')] <- '_into'
  custom_dbWriteTable(name_pk = 'grpid', fks=NA, table_name = paste0('c_peak_groups_',MS1name,"_",MS2name), df=grp_peaklist, con=con)

  ###############################################
  # Add s_peak_meta (i.e. scan information)
  ###############################################
  dropc <- c('filename')
  scan_info <- scan_info[,!colnames(scan_info) %in% dropc]
  scan_info <- update_cn_order(name_pk = 'pid',names_fk= 'fileid', df = scan_info)
  custom_dbWriteTable(name_pk = 'pid', fks=fks_fileid, table_name = paste0('s_peak_meta_',MS1name,"_",MS2name), df=scan_info, con=con)

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
                      table_name = paste0('s_peaks_',MS1name,"_",MS2name), df=scanpeaks_frag, con=con)


  ###############################################
  # Add MANY-to-MANY links for c_peak to c_peak_group
  ###############################################
  c_peak_X_c_peak_group <- get_group_peak_link(xset)

  fks_for_cxg <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'),
                      'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks')
                      )

  custom_dbWriteTable(name_pk = 'cXg_id', fks=fks_for_cxg,
                      table_name = paste0('c_peak_X_c_peak_group_',MS1name,"_",MS2name), df=c_peak_X_c_peak_group, con=con)

  if (pa@f4f_link_type=='individual'){
    ###############################################
    # Add MANY-to-MANY links for c_peak to s_peak_meta
    ###############################################
    grpdf <- grped_df
    c_peak_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'cid')])
    c_peak_X_s_peak_meta <- cbind('cXp_id'=1:nrow(c_peak_X_s_peak_meta), c_peak_X_s_peak_meta)

    fks_for_cXs <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'),
                        'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks'))

    custom_dbWriteTable(name_pk = 'cXp_id', fks=fks_for_cXs,
                        table_name = paste0('c_peak_X_s_peak_meta_',MS1name,"_",MS2name), df=c_peak_X_s_peak_meta, con=con)
  }else{
    ###############################################
    # Add MANY-to-MANY links for c_peak_group to s_peak_meta
    ###############################################
    grpdf <- grped_df
    c_peak_group_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'grpid')])
    c_peak_group_X_s_peak_meta <- cbind('gXp_id'=1:nrow(c_peak_group_X_s_peak_meta), c_peak_group_X_s_peak_meta)

    fks_for_cXs <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'),
                        'grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'))

    custom_dbWriteTable(name_pk = 'gXp_id', fks=fks_for_cXs,
                        table_name = paste0('c_peak_group_X_s_peak_meta_',MS1name,"_",MS2name), df=c_peak_group_X_s_peak_meta, con=con)


  }

  if (length(pa@av_spectra)>0){

    av_spectra <- plyr::ldply(pa@av_spectra, get_av_spectra_for_db)

    # for some reason the names are not being saved for the list as a column, so we just get them back
    colnames(av_spectra)[1] <- 'grpid'
    av_spectra$grpid <- names(pa@av_spectra)[av_spectra$grpid]

    colnames(av_spectra)[2] <- 'fileid'
    av_spectra$avid <- 1:nrow(av_spectra)
    fks_for_av_spectra <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'),
                               'grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups')
                               )

    custom_dbWriteTable(name_pk = 'avid', fks=fks_for_av_spectra,
                        table_name = paste0('av_peaks_',MS1name,"_",MS2name), df=av_spectra, con=con)

  }


  if (!is.null(xsa)){
    ###############################################
    # Add CAMERA ruleset
    ###############################################
    if(is.null(xsa@ruleset)){
      rules_pos <- read.table(system.file(file.path('rules', 'extended_adducts_pos.csv') , package = "CAMERA"), header = T)
      rules_neg <- read.csv(system.file(file.path('rules', 'extended_adducts_neg.csv') , package = "CAMERA"))
      rules <- rbind(rules_pos, rules_neg)

    }else{
      rules <- xsa@ruleset
    }
    rules$rule_id <- 1:nrow(rules)
    custom_dbWriteTable(name_pk = 'rule_id', fks=NA,
                        table_name = paste0('adduct_rules_',MS1name,"_",MS2name), df=rules, con=con)

    ###############################################
    # Add neutral mass groups
    ###############################################
    annoGrp <- data.frame(xsa@annoGrp)
    colnames(annoGrp)[1] <- 'nm_id'
    custom_dbWriteTable(name_pk = 'nm_id', fks=NA,
                        table_name = paste0('neutral_masses_',MS1name,"_",MS2name), df=annoGrp, con=con)

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
                        table_name = paste0('adduct_annotations_',MS1name,"_",MS2name), df=annoID, con=con)

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
                        table_name = paste0('isotope_annotations_',MS1name,"_",MS2name), df=isoID, con=con)
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

  if(DBI::dbExistsTable(con,table_name)){
    DBI::dbRemoveTable(con,table_name)
  }
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
