#' @title Spectral matching
#'
#' @description
#' Perform spectral matching to spectral libraries using dot product cosine on a LC-MS/MS dataset and link to XCMS features.
#'
#' @param pa purityA object; Needs to be the same used for frag4feature function
#' @param xset xcms object; Needs to be the same used for frag4feature function
#' @param ra_thres_t numeric; Relative abundance threshold for target spectra
#'                      (Peaks below this RA threshold will be excluded)
#' @param ra_thres_l numeric; Relative abundance threshold for library spectra
#' @param cores numeric; Number of cores to use
#' @param pol character; Polarity ['positive' or 'negative']
#' @param ppm_tol_prod numeric; PPM tolerance to match to product
#' @param ppm_tol_prec numeric; PPM tolerance to match to precursor
#' @param score_thres numeric; Dot product cosine score threshold
#' @param out_dir character; Out directory for the SQLite result database
#' @param topn numeric [optional]; Only use top n matches
#' @param db_name character [optional]; Name of the result database
#' @param grp_peaklist dataframe [optional]; Can use any peak dataframe. Still needs to be derived from the xset object though
#'                                   (e.g. can use CAMERA peaklist)
#' @param library_db_pth character [optional]; path to library spectral SQLite database. Defaults to msPurityData package data.
#' @param instrument_types vector [optional]; Vector of instrument types, defaults to all
#' @param library_sources vector [optional]; Vector of library sources. Default option is for massbank only but the 'lipidblast'
#'                                    library is also available
#' @param scan_ids vector [optional]; Vector of unique scan ids calculated from msPurity "pid". These scans will on
#'                        used for the spectral matching. All scans will be used if set to NA
#'
#'
#' @return list of database details and dataframe summarising the results for the xcms features
#' @examples
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- frag4feature(pa, xset)
#' #NOTE that scan_ids here are refer the unique scan id calculated by purityA (pids).
#' #Only required if you want to limit the spectral matching to certain scans
#' result <- spectral_matching(pa, xset, scan_ids = c(1120,  366, 1190, 601,  404,1281, 1323, 1289))
#' @export
spectral_matching <- function(pa, xset, ra_thres_l=0, ra_thres_t=2, cores=1, pol='positive', ppm_tol_prod=10, ppm_tol_prec=5,
                                     score_thres=0.6, out_dir='.', topn=NA, db_name=NA, grp_peaklist=NA, library_db_pth=NA,
                                     instrument_types=NA, library_sources='massbank', scan_ids=NA){
  message("Running msPurity spectral matching function for LC-MS(/MS) data")

  if (is.na(db_name)){
    db_name <- paste('lcmsms_data', format(Sys.time(), "%Y-%m-%d-%I%M%S"), '.sqlite', sep="-")
  }

  if (!is.data.frame(grp_peaklist)){
    grp_peaklist <- xcms::peakTable(xset)
    grp_peaklist <- data.frame(cbind('grpid'=1:nrow(grp_peaklist), grp_peaklist))
  }

  if (is.na(library_db_pth)){
    library_db_pth <- system.file("extdata", "library_spectra", "library_spectra.db", package="msPurityData")
  }


  ########################################################
  # Export the target data into sqlite database
  ########################################################
  message("Creating a database of fragmentation spectra and LC features")
  target_db_pth <- export_2_sqlite(pa, grp_peaklist, xset, out_dir, db_name)

  ########################################################
  # Perform the spectral matching
  ########################################################
  message("Performing spectral matching")
  matched <- match_2_library(target_db_pth,
                  library_db_pth,
                  ra_thres_t=ra_thres_t,
                  ra_thres_l=ra_thres_l,
                  cores=cores,
                  pol=pol,
                  ppm_tol_prod=ppm_tol_prod,
                  ppm_tol_prec=ppm_tol_prec,
                  topn=topn,
                  instrument_types = instrument_types,
                  library_sources = library_sources,
                  scan_ids = scan_ids)

  ########################################################
  # Create a summary table for xcms grouped objects
  ########################################################
  if (matched){
    message("Summarising LC features annotations")
    xcms_summary_df <- get_xcms_sm_summary(target_db_pth, topn=topn, score_f=score_thres)
  }else{
    xcms_summary_df <- NA
  }


  return(list('db_name' = db_name, 'out_dir' = out_dir, 'xcms_summary_df' = xcms_summary_df))
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


ra_calc <- function(x){
  x$ra <- (x$i/max(x$i))*100
  return(x)
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

get_target_spectra_list <- function(x, scan_info){
  list(x, scan_info[scan_info$pid==unique(x$pid),]$precursorMZ)
}

match_2_library <- function(target_db_pth, library_db_pth, instrument_types=NA, mslevel=NA, mslevel_match=TRUE,
                            ra_thres_t=2, ra_thres_l=0, cores=1, pol, ppm_tol_prod=100, ppm_tol_prec=50, topn=5,
                            library_sources=NA, scan_ids=NA){




  ########################################################
  # Get target spectra
  ########################################################
  # Get all of the target MS2 data
  conT <- DBI::dbConnect(RSQLite::SQLite(), target_db_pth)

  # Get scan peaks for target spectra
  target_spectra <- DBI::dbGetQuery(conT, 'SELECT * FROM s_peaks ' )

  if (!anyNA(scan_ids)){
    target_spectra <- target_spectra[target_spectra[,'pid'] %in% scan_ids,]
  }

  # Calculate relative abundance
  target_spectra <- plyr::ddply(target_spectra, ~ pid, ra_calc)

  # Relative abundnace threshold
  target_spectra <- target_spectra[target_spectra$ra>=ra_thres_t,]

  # Flag for spectra type
  target_spectra$type <- 1


  ########################################################
  # Get library spectra
  ########################################################
  # Get all of the library peaks based on parameters given
  conL <- DBI::dbConnect(RSQLite::SQLite(), library_db_pth)

  # Match each target spectra to the library spectra
  # Only keep matches with certain criteria certain criteria

  if (anyNA(library_sources)){
    library_meta_query <- sprintf("SELECT * FROM library_spectra_meta WHERE polarity = '%s'", pol)
  }else{
    l_source_str <- paste("'",paste(library_sources, collapse = "', '"), "'", sep='')
    library_meta_query <- sprintf("SELECT m.*, s.name AS source_name FROM library_spectra_meta as m
                                    LEFT JOIN library_spectra_source AS s ON s.id=m.library_spectra_source_id
                                    WHERE m.polarity = '%s' AND s.name IN (%s)", pol, l_source_str)
  }


  if (!anyNA(instrument_types)){
    # instrument_types <- c('CE-ESI-TOF', 'ESI-ITFT', 'ESI-ITTOF', 'ESI-QTOF', 'LC-ESI-IT',
    #                       'LC-ESI-ITFT', 'LC-ESI-ITTOF','LC-ESI-QFT', 'LC-ESI-QIT', 'LC-ESI-QQ', 'LC-ESI-QTOF', 'LC-ESI-TOF')
    ints_string <- paste("'",paste(instrument_types, collapse = "', '"), "'", sep='')
    library_meta_query <- paste(library_meta_query, sprintf(" AND instrument_type IN (%s)", ints_string))

  }


  library_meta <- DBI::dbGetQuery(conL, library_meta_query)

  if (!is.na(mslevel)){
    library_meta <- library_meta[library_meta$ms_level==mslevel,]
  }

  library_meta_ids <- paste(library_meta$id, collapse=", ")

  library_spectra <- DBI::dbGetQuery(conL, sprintf('SELECT * FROM library_spectra WHERE library_spectra_meta_id IN (%s)', library_meta_ids) )

  library_spectra <- plyr::ddply(library_spectra, ~ library_spectra_meta_id, ra_calc)

  library_spectra <- library_spectra[library_spectra$ra>ra_thres_l,] # mass bank default does not do this filter

  library_spectra$type <- 2


  ########################################################
  # Weight spectra
  ########################################################
  # weighted intensity value
  target_spectra$w <- (target_spectra$ra^0.5)*(target_spectra$mz^2)
  library_spectra$w <- (library_spectra$ra^0.5)*(library_spectra$mz^2)


  ########################################################
  # Get target spectra in format ready for matching
  ########################################################
  # Order by relative abundance for target spectra
  target_spectra <- target_spectra[order(target_spectra[,'pid'], -target_spectra[,'ra']), ]

  scan_info <- DBI::dbGetQuery(conT,'SELECT * FROM s_peak_meta' )

  # we loop through the target spectra as list
  target_spectra_list <- plyr::dlply(target_spectra, ~ pid, get_target_spectra_list,
                                     scan_info=scan_info)


  ########################################################
  # Perform spectral matching
  ########################################################
  # Convert library spectra to matrix for speed
  library_spectra <- as.matrix(library_spectra[,c('mz','ra','type','w','library_spectra_meta_id')])

  if (cores>1){
    cl<-parallel::makeCluster(cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    parallel = TRUE
  }else{
    parallel = FALSE
  }

  allmatches_l <- plyr::llply(target_spectra_list, .parallel = parallel, match_targets,
                        library_spectra=library_spectra,
                        library_precs=library_meta$precursor_mz,
                        ppm_tol_prod=ppm_tol_prod,
                        ppm_tol_prec=ppm_tol_prec)

  allmatches <- plyr::ldply(allmatches_l, .id = 'pid')


  ########################################################
  # Submit to database
  ########################################################
  if (nrow(allmatches)>0){
    if (!is.na(topn)){
      allmatches <- plyr::ddply(allmatches, ~ pid, get_topn)
    }

    library_meta_f <- library_meta[library_meta$id %in% allmatches$lid, ]

    colnames(library_meta_f)[which(colnames(library_meta_f)=='id')] = 'lid'

    custom_dbWriteTable(name_pk = 'lid', fks = NA,
                        df=library_meta_f, table_name = 'library_meta', con = conT)

    allmatches$mid <- 1:nrow(allmatches)

    fks_lid <- list('lid'=list('new_name'='lid', 'ref_name'='lid', 'ref_table'='library_spectra_meta'))
    fks_pid <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'))

    custom_dbWriteTable(name_pk = 'mid', fks = append(fks_lid, fks_pid),
                        df=allmatches, table_name = 'matches', con = conT)
    matched = TRUE
  }else{
    matched = FALSE
  }

  DBI::dbDisconnect(conT)
  DBI::dbDisconnect(conL)

  return(matched)

}

get_topn <- function(x){
  if(nrow(x)>topn){
    return(x[topn,])
  }else{
    return(x)
  }

}

get_xcms_sm_summary <- function(target_db_pth, topn=NA, score_f=0.3, frag_nm_f=1){

  # Get all of the target fragmentation data
  conT <- DBI::dbConnect(RSQLite::SQLite(), target_db_pth)

  XLI <- DBI::dbGetQuery(conT, 'SELECT * FROM c_peak_groups
   LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.grpid=c_peak_groups.grpid
   LEFT JOIN c_peaks on c_peaks.cid=cXg.cid
   LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.cid=c_peaks.cid
   LEFT JOIN s_peak_meta ON cXs.pid=s_peak_meta.pid
   LEFT JOIN matches ON matches.pid=s_peak_meta.pid
   LEFT JOIN library_meta ON matches.lid=library_meta.lid
   WHERE matches.score IS NOT NULL')

  XLI <- XLI[XLI$score>=score_f & XLI$match>=frag_nm_f,]

  XLI <- XLI[order(XLI$grpid, -XLI$score),]

  # Only use topn of hits per scan
  if (!is.na(topn)){
    XLI <- plyr::ddply(XLI, ~ pid, check_topn, topn=topn)
  }

  # Summaries the annotation hits
  xcms_mtch <- plyr::ddply(XLI, ~ grpid, get_ann_summary)

  if(nrow(xcms_mtch)==0){
    print('NO MATCHES FOR XCMS')
    DBI::dbDisconnect(conT)
    return(0)
  }

  DBI::dbWriteTable(conT, name='xcms_match', value=xcms_mtch, row.names=F, append=T)

  return(xcms_mtch)


}

check_topn <- function(x, topn){
  if(nrow(x)>topn){
      return(x[topn,])
    }else{
      return(x)
   }
}



median_match_results <- function(y){
  c('best_median_score'=median(y$score), 'best_median_perc_mtch'=median(y$perc_mtch), 'best_median_match'=median(y$match))
}

get_ann_summary <- function(x){
  # get the 'best' match based on the best scored compounds with
  # the same name
  med_results <- plyr::ddply(x, ~ name, median_match_results)

  med_results <- med_results[order(med_results$best_median_score, decreasing = TRUE),]
  colnames(med_results)[1] <- 'best_name'
  best_match <- med_results[1,]

  unique_names <- med_results$name

  # Get all the matches
  allnames <- paste(unlist(x$name),collapse=", ")
  allpids <- paste(unlist(x$pid),collapse=", ")
  alllids <- paste(unlist(x$lid),collapse=", ")
  allscores <- paste(unlist(round(x$score,3)),collapse=", ")

  snm <- length(unique(unlist(x$pid)))

  out_v <- c(unlist(best_match), 'all_names'=allnames, 'all_pids'=allpids,
             'all_lids'=alllids, 'all_scores'=allscores, 'numb_unique_scans'=snm)

}

match_targets <- function(target_peaks_list, library_spectra, ppm_tol_prod=100, ra_diff=10, library_precs, ppm_tol_prec=50){

  # Get target peaks and precursor of target
  target_peaks <- target_peaks_list[[1]]
  target_prec <- target_peaks_list[[2]]

  tcnm = c('mz','ra','type','w')

  # Convert target peaks to matrix for speed
  if(is.vector(target_peaks) | nrow(target_peaks)==1){
    target_peaks <- target_peaks[tcnm]
    target_peaks <- matrix(as.numeric(target_peaks), nrow=1, ncol=4)
    colnames(target_peaks) <- tcnm
  }else{
    pid <- unique(target_peaks[,'pid'])
    target_peaks <- as.matrix(target_peaks[,tcnm])
  }

  # get unique ids of library_spectra precursor 'meta' info
  library_meta_ids <- unique(library_spectra[,'library_spectra_meta_id'])

  # Get ppm difference between precursor from library and target
  out_prec <- outer(target_prec, library_precs, '-')
  colprec <- library_precs[col(out_prec)] # so we can divide by column if library
  ppmdiff_prec <- as.vector(abs(1e6*(out_prec/colprec)))

  # Filter the library spectra that does not meet the precursor tolerance check
  library_spectra_red <- library_spectra[library_spectra[,'library_spectra_meta_id'] %in% library_meta_ids[ppmdiff_prec<=ppm_tol_prec],]

  # Exit if no peaks left (if vector it means it is just 1 row)
  if(is.vector(library_spectra_red)){
    lnms <- names(library_spectra_red)
    library_spectra_red <- matrix(library_spectra_red, nrow=1, ncol=length(library_spectra_red))
    colnames(library_spectra_red) <- lnms
  }else if (nrow(library_spectra_red)==0){
    return(NULL)
  }

  # calculate ppm error of all peaks (I am presuming this is faster than doing it library spectra at a time...)
  out_peaks_mz <- outer(target_peaks[,'mz'], library_spectra_red[,'mz'], '-')
  colmz <- library_spectra_red[,'mz'][col(out_peaks_mz)] # so we can calculate column wide
  ppmdiff_all <- abs(1e6*(out_peaks_mz/colmz))

  # Calculate the percentage error of RA to library
  out_peaks_ra <- outer(target_peaks[,'ra'], library_spectra_red[,'ra'], '-')
  colra <- library_spectra_red[,'ra'][col(out_peaks_ra)] # so we can calculate column wide
  idiff_all <- abs(outer(target_peaks[,'ra'], library_spectra_red[,'ra'], '-')/colra)*100

  # Get index for difference matrices
  gidx <- c(0,cumsum(table(library_spectra_red[,'library_spectra_meta_id'])))

  lnm <- length(unique(library_spectra_red[,'library_spectra_meta_id']))


  l = list()

  library_meta_info_ids <- unique(library_spectra_red[,'library_spectra_meta_id'])

  for(j in 1:lnm){

    uid <- library_meta_info_ids[j]

    mtch<- matchi(ppmdiff = ppmdiff_all[,(gidx[j]+1):gidx[j+1]],
                  idiff = idiff_all[,(gidx[j]+1):gidx[j+1]],
                  library_peaks = library_spectra_red[library_spectra_red[,'library_spectra_meta_id']==uid,],
                  target_peaks = target_peaks,
                  ppm_tol_prod=ppm_tol_prod,
                  ra_diff=ra_diff)

    l[[j]] <- c(mtch, uid)
  }

  d <-plyr::ldply(l)


  if(nrow(d)==0){
    return(NULL)
  }else{

    colnames(d) <- c('score','perc_mtch','match','lid')
    d <- d[d[,'score']>0 & !is.na(d[,'score']),]
    return(d)
  }

}


matchi <-function(library_peaks, target_peaks, ppmdiff, idiff, ppm_tol_prod=50, ra_diff=10){

  if (!any(ppmdiff<ppm_tol_prod)){
    return(c(NA,NA,NA))
  }

  if(is.vector(library_peaks)){
    ppmdiff <- matrix(ppmdiff)
    idiff <- matrix(idiff)
    library_peaks <- library_peaks[c('mz','ra','type','w')]
    library_peaks <- matrix(library_peaks, nrow=1, ncol=4)
    colnames(library_peaks) <- c('mz','ra','type','w')

  }else{
    library_peaks <- library_peaks[,c('mz','ra','type','w')]
  }

  # Following the pMatch-hammer method for peak matching but with slight variation that we also check the percentage difference for
  # the relative intensity as well
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3928522/

  # we keep a vector detailing which library peaks we have already matched
  lp_remain <- seq(1, ncol(idiff))

  # Loop through the differences for the target to library
  allpeaks <-  list()
  mcount = NULL

  for (i in 1:nrow(idiff)){

    ppmD <- ppmdiff[i,]
    iD <- idiff[i,]

    if(all(is.na(ppmD))){
      allpeaks[[i]] <- rbind(target_peaks[i,],c(0, 0, 2,0))
      mcount <- append(mcount, 'miss')
      next
    }

    bm <- min(ppmD, na.rm = T)

    if (bm>ppm_tol_prod){
      allpeaks[[i]] <- rbind(target_peaks[i,],c(0, 0, 2,0))
      mcount <- append(mcount, 'miss')
      next
    }

    # First check to see if there is a matching intensity value within ra_diff (default 10%)
    intenc <- iD[ppmD<ppm_tol_prod & iD<ra_diff & !is.na(ppmD) & !is.na(ra_diff)]

    if (!identical(unname(intenc), numeric(0))){
      matchi <- match(min(intenc, na.rm = T), iD)
    }else{
      matchi <- match(min(ppmD[ppmD<ppm_tol_prod], na.rm=T), ppmD)
    }

    allpeaks[[i]] <- rbind(target_peaks[i,],library_peaks[matchi,])
    mcount <- append(mcount, 'match')
    lp_remain[matchi] <- NA
    idiff[,matchi] <- NA
    ppmdiff[,matchi] <- NA
  }

  c = length(allpeaks)+1

  for (i in lp_remain){
    if(is.na(i)){
      next
    }
    allpeaks[[c]] <- rbind(c(0, 0, 1,0), library_peaks[i,])
    mcount <- append(mcount, 'miss')
    c = c+1
  }

  allpeaksm <- do.call(rbind, allpeaks)

  wt <- as.numeric(allpeaksm[allpeaksm[,'type']==1,][,'w'])

  wl <- as.numeric(allpeaksm[allpeaksm[,'type']==2,][,'w'])

  cossim_out <- CosSim(wt, wl)

  percentage_match <- sum(mcount=='match')/length(mcount)


  return(c(cossim_out,percentage_match,sum(mcount=='match')))

}


# group1d <- function(v, thr){
#
#   v <- sort(v)
#   a <- ''
#   cl <-  1
#   c <- 1
#
#   for(i in 1:(length(v)-1)){
#     a[i] <- (1e6*(v[i+1]-v[i]))/v[i]
#
#     if(as.numeric(a[i])<thr){
#       cl[i+1] <- c
#     }else{
#       c <- c+1
#       cl[i+1] <- c
#     }
#
#   }
#   return(cl)
# }
#
#
# alignPeaks <- function(v, type, thr){
#
#   v <- sort(v)
#   a <- ''
#   p <- ''
#   c <- 1
#
#   for(i in 1:(length(v)-1)){
#     if (!ty[i]==ty[i+1]){
#       a[i] <- (1e6*(v[i+1]-v[i]))/v[i]
#       if (a[i]>thr){
#         p[i] = p
#         p[i+1] = p
#       }else{
#
#       }
#       print(paste(ty[i], ty[i+1], a[i]))
#     }else{
#       print(i)
#     }
#
#
#   }
#   return(cl)
# }


CosSim <- function(A,B) {
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}





