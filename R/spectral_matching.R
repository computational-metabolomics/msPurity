#' @title Spectral matching for LC-MSMS
#'
#' @description
#' The function works for LC-MS or DI-MS datasets.
#'
#' @param filePth character = Path of the file to be processed
#' @return  dataframe of the median mz, intensity, signal-to-noise ratio.
#' @examples
#' mzmlPth <- system.file("extdata", "dims", "mzML", "B02_Daph_TEST_pos.mzML", package="msPurityData")
#' avP <- averageSpectraSingle(mzmlPth)
#'
#' Perform purity assessments
#' pa <- purityA(xset@filepaths, cores=params$nSlaves, interpol = "linear", iwNorm = TRUE, ilim = 0.05)
#' pa <- frag4feature(pa, xset)
#' spectral_matching(pa, grp_peaklist, xset, out_dir='.',db_name)
#' @export
spectral_matching_lcmsms <- function(pa, xset, ra_thres=2, cores=1, pol='positive', ppm_tol_MSMS=10, ppm_tol_prec=5,
                                     topn=50, score_thres=0.6, out_dir='.', db_name=NA, grp_peaklist=NA, library_db_pth=NA,
                                     instrument_types=NA){
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

  target_db_pth <- export_2_sqlite(pa, grp_peaklist, xset, out_dir, db_name)

  match_2_library(target_db_pth,
                  library_db_pth,
                  ra_thres=ra_thres,
                  cores=cores,
                  pol=pol,
                  ppm_tol_MSMS=ppm_tol_MSMS,
                  ppm_tol_prec=ppm_tol_prec,
                  topn=topn,
                  instrument_types = instrument_types)

  get_topn_xcms_table(target_db_pth, topn=topn, score_f=score_thres)

  return('db_name'=db_name)
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

update_cn_order <- function(pk_name, fk_names, df){
  # primary key needs to be at the start
  # foreign keys at the end
  cn <- colnames(df)

  if (anyNA(fk_names)){
    columnorder <- c(pk_name, cn[!cn %in% pk_name])
  }else{
    columnorder <- c(pk_name, cn[!cn %in% c(pk_name, fk_names)], fk_names)
  }
  return(df[,columnorder])
}

export_2_sqlite <- function(pa, grp_peaklist, xset, out_dir, db_name, explicit=TRUE){

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

  file_query <- get_create_query(pk='fileid', fks=NA, table_name='fileinfo', df=filedf)

  sqr <- DBI::dbSendQuery(con, file_query)
  DBI::dbClearResult(sqr)

  DBI::dbWriteTable(con, name='fileinfo', value=filedf, row.names=F, append=T)

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

  fks <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'))
  c_peak_query <- get_create_query(pk='cid', fks=fks, table_name='c_peaks', df=c_peaks)

  sqr <- DBI::dbSendQuery(con, c_peak_query)
  DBI::dbClearResult(sqr)
  DBI::dbWriteTable(con, name='c_peaks', value=c_peaks, row.names=F, append=T)



  ###############################################
  # Add c_peak_groups (i.e. XCMS grouped peaks)
  ###############################################
  if (is.matrix(grp_peaklist)){
    grp_peaklist <- data.frame(grp_peaklist)
  }
  grp_peaklist <- update_cn_order(pk_name = 'grpid',fk_names= NA, df = grp_peaklist)

  grp_query <- get_create_query(pk='grpid', fks=NA, table_name='c_peak_groups', df=grp_peaklist)

  sqr <- DBI::dbSendQuery(con, grp_query)
  DBI::dbClearResult(sqr)
  DBI::dbWriteTable(con, name='c_peak_groups', value=grp_peaklist, row.names=F, append=T)

  ###############################################
  # Add s_peak_meta (i.e. scan information)
  ###############################################
  dropc <- c('filename')
  scan_info <- scan_info[,!colnames(scan_info) %in% dropc]
  scan_info <- update_cn_order(pk_name = 'pid',fk_names= 'fileid', df = scan_info)

  fks <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'))
  s_peak_meta_query <- get_create_query(pk='pid', fks=fks, table_name='s_peak_meta', df=scan_info )

  sqr <- DBI::dbSendQuery(con, s_peak_meta_query)
  DBI::dbClearResult(sqr)

  DBI::dbWriteTable(con, name='s_peak_meta', value=scan_info, row.names=F, append=T)

  ###############################################
  # Add s_peaks (i.e. the mz, i from each scan)
  ###############################################
  # ensure the filedf is in the same order as the scan_info file
  # other wise the s_peak_meta ids might not match up
  filedf <- filedf[as.numeric(unique(scan_info$fileid)),]

  scanpeaks_frag <- plyr::ddply(filedf, ~ fileid, function(x){
    mr <- mzR::openMSfile(as.character(x$filepth))
    scanpeaks <- mzR::peaks(mr)
    scans <- mzR::header(mr)
    names(scanpeaks) <- seq(1, length(scanpeaks))
    scanpeaks_df <- plyr::ldply(scanpeaks[scans$seqNum[scans$msLevel>1]], .id=T)
  })

  comb <- paste(scanpeaks_frag[,1], scanpeaks_frag[,2], sep=' ')
  scanpeaks_frag <- cbind(1:nrow(scanpeaks_frag), cumsum(!duplicated(comb)), scanpeaks_frag)

  colnames(scanpeaks_frag) <- c('sid','pid', 'fileid', 'scan', 'mz', 'i')

  scanpeaks_frag <- update_cn_order(pk_name = 'sid',fk_names= c('fileid','pid'), df = scanpeaks_frag)

  fks <- list('fileid'=list('new_name'='fileid', 'ref_name'='fileid', 'ref_table'='fileinfo'),
              'pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'))

  s_peak_query <- get_create_query(pk='sid', fks=fks, table_name='s_peaks', df=scanpeaks_frag)

  sqr <- DBI::dbSendQuery(con, s_peak_query)
  DBI::dbClearResult(sqr)

  DBI::dbWriteTable(con, name='s_peaks', value=scanpeaks_frag, row.names=F, append=T)

  ###############################################
  # Add MANY-to-MANY links for c_peak to c_peak_group
  ###############################################
  grpdf <- pa@grped_df
  c_peak_X_c_peak_group <- unique(grpdf[ ,c('grpid', 'cid')])
  c_peak_X_c_peak_group <- cbind('cXg_id'=1:nrow(c_peak_X_c_peak_group), c_peak_X_c_peak_group)

  fks <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'),
              'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks'))
  c_peak_X_c_peak_query <- get_create_query(pk='cXg_id', fks=fks, table_name='c_peak_X_c_peak_group', df=c_peak_X_c_peak_group)
  sqr <- DBI::dbSendQuery(con, c_peak_X_c_peak_query)
  DBI::dbClearResult(sqr)

  DBI::dbWriteTable(con, name='c_peak_X_c_peak_group', value=c_peak_X_c_peak_group, row.names=F, append=T)

  ###############################################
  # Add MANY-to-MANY links for c_peak to s_peak_meta
  ###############################################
  c_peak_X_s_peak_meta <- unique(grpdf[ ,c('pid', 'cid')])
  c_peak_X_s_peak_meta <- cbind('cXp_id'=1:nrow(c_peak_X_s_peak_meta), c_peak_X_s_peak_meta)

  fks <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'),
              'cid'=list('new_name'='cid', 'ref_name'='cid', 'ref_table'='c_peaks'))

  c_peak_X_s_peak_meta_query <- get_create_query(pk='cXp_id', fks=fks, table_name='c_peak_X_s_peak_meta', df=c_peak_X_s_peak_meta)
  sqr <- DBI::dbSendQuery(con, c_peak_X_s_peak_meta_query)
  DBI::dbClearResult(sqr)

  DBI::dbWriteTable(con, name='c_peak_X_s_peak_meta', value=c_peak_X_s_peak_meta, row.names=F, append=T)
  DBI::dbDisconnect(con)

  return(db_pth)



}

ra_calc <- function(x){
  x$ra <- (x$i/max(x$i))*100
  return(x)
}

custom_dbWriteTable <- function(pk_name, fks, df, table_name){

  df <- update_cn_order(pk_name=pk_name, fk_names=names(fks), df = df)
  query <- get_create_query(pk=pk_name, fks=fks, table_name=table_name, df=df)
  sqr <- DBI::dbSendQuery(con, df)
  DBI::dbWriteTable(con, name=table_name, value=library_meta_f, row.names=T, append=T)


}

get_target_spectra_list <- function(x){
  list(x, scan_info[scan_info$pid==unique(x$pid),]$precursorMZ)
}

match_2_library <- function(target_db_pth, library_db_pth, instrument_types=NA, mslevel=NA, mslevel_match=TRUE,
                            ra_thres=2, cores=1, pol, ppm_tol_MSMS=100, ppm_tol_prec=50, topn=5){

  if (anyNA(instrument_types)){
    instrument_types <- c('CE-ESI-TOF', 'ESI-ITFT', 'ESI-ITTOF', 'ESI-QTOF', 'LC-ESI-IT',
                          'LC-ESI-ITFT', 'LC-ESI-ITTOF','LC-ESI-QFT', 'LC-ESI-QIT', 'LC-ESI-QQ', 'LC-ESI-QTOF', 'LC-ESI-TOF')
  }

  #library_db <- 'W:\\users\\tnl495\\massbank\\massbank.db'
  # Get all of the target MS2 data
  conT <- DBI::dbConnect(RSQLite::SQLite(), target_db_pth)

  # Get all of the library peaks based on parameters given
  conL <- DBI::dbConnect(RSQLite::SQLite(), library_db_pth)

  target_spectra <- DBI::dbGetQuery(conT, 'SELECT * FROM s_peaks' )

  # get relative abundance
  target_spectra <- plyr::ddply(target_spectra, ~ pid, ra_calc)

  # Match each target spectra to the library spectra
  # Only keep matches above certain criteria

  ints_string <- paste("'",paste(instrument_types, collapse = "', '"), "'", sep='')
  library_meta_query <- sprintf("SELECT * FROM meta_info WHERE polarity = '%s' AND instrument_type IN (%s)", pol, ints_string)
  library_meta <- DBI::dbGetQuery(conL, library_meta_query)

  if (!is.na(mslevel)){
    library_meta <- library_meta[library_meta$ms_level==mslevel,]
  }

  #prec_ppm <- abs(1e6*outer(scan_info$precursorMZ, meta_f$precursor_mz, '-')/scan_info$precursorMZ)

  #pf <- apply(prec_ppm, 1, min)<=50
  #nrow(scan_info[pf,])
  #prec_ppm <- cbind(prec_ppm, )
  #prec_ppm[]

  library_meta_ids <- paste(library_meta$UID,collapse=", ")

  library_spectra <- DBI::dbGetQuery(conL, sprintf('SELECT * FROM spectra WHERE meta_info_uid IN (%s)', library_meta_ids) )

  # Will update this

  library_spectra <- plyr::ddply(library_spectra, ~ meta_info_uid, ra_calc)


  # Relative abundnace threshold
  target_spectra <- target_spectra[target_spectra$ra>=ra_thres,]
  #library_spectra <- library_spectra[library_spectra$ra>ra_thres,] # mass bank default does not do this filter

  # Flag for spectra type
  target_spectra$type <- 1
  library_spectra$type <- 2

  # weighted intensity value
  target_spectra$w <- (target_spectra$ra^0.5)*(target_spectra$mz^2)
  library_spectra$w <- (library_spectra$ra^0.5)*(library_spectra$mz^2)

  # Order by relative abundance for target spectra
  target_spectra <- target_spectra[order(target_spectra[,'pid'], -target_spectra[,'ra']), ]

  scan_info <- DBI::dbGetQuery(conT,'SELECT * FROM s_peak_meta' )

  # we loop through the target spectra as list
  target_spectra_list <- plyr::dlply(target_spectra, ~ pid, get_target_spectra_list)

  # Convert library spectra to matrix for speed
  library_spectra <- as.matrix(library_spectra[,c('mz','ra','type','w','meta_info_uid')])


  if (cores>1){
    operator <- foreach::'%dopar%'
    # The cores are split across the files. If enough cores available they
    # are also split across the clustering algorithm as well

    cl<-parallel::makeCluster(cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    parallel = TRUE

  }else{
    operator <- foreach::'%do%'
    parallel = FALSE
  }

  #tnm <- length(target_spectra_list)
  #pb <- txtProgressBar(max=tnm, style=3)
  #progress <- function(n) setTxtProgressBar(pb, n)
  #opts <- list(progress=progress)

  allmatches_l <- plyr::llply(target_spectra_list,.parallel = parallel, match_targets,
                        library_spectra=library_spectra,
                        library_precs=library_meta$precursor_mz,
                        ppm_tol_MSMS=ppm_tol_MSMS,
                        ppm_tol_prec=ppm_tol_prec)

  allmatches <- plyr::ldply(allmatches_l, .id = 'pid')

  if (!is.na(topn)){
    allmatches <- plyr::ddply(allmatches, ~ pid, get_topn)
  }

  library_meta_f <- library_meta[library_meta$UID %in% allmatches$lid, ]

  library_meta_f <- update_cn_order(pk_name='UID',fk_names= NA, df = library_meta_f)

  library_meta_f_query <- get_create_query(pk='UID', fks=NA, table_name='library_meta', df=library_meta_f)
  sqr <- DBI::dbSendQuery(conT, library_meta_f_query)

  DBI::dbWriteTable(conT, name='library_meta', value=library_meta_f, row.names=T, append=T)

  allmdf2 <- merge(allmatches, meta, by.x = 'lid', by.y='UID')

  allmdf2 <- merge(allmdf2, scan_info , by.x = 'pid', by.y='pid')

  allmdf2 <- allmdf2[order(allmdf2$pid),]

  dbWriteTable(conT, name='matches', value=allmdf2, row.names=T, append=T)

  dbDisconnect(conT)
  dbDisconnect(conL)

}

get_topn <- function(x){
  if(nrow(x)>topn){
    return(x[topn,])
  }else{
    return(x)
  }

}

get_topn_xcms_table <- function(target_db_pth, topn, score_f=0.3, frag_nm_f=1){

  #library_db <- 'W:\\users\\tnl495\\massbank\\massbank.db'
  # Get all of the target MS2 data
  conT <- dbConnect(RSQLite::SQLite(), target_db_pth)
  MA <- dbGetQuery(conT, 'SELECT * FROM matches')



  pids <- paste(MA$pid, collapse=", ")

  XLI <- dbGetQuery(conT, 'SELECT matches.*, xcms_msms_link.grpid, xcms_msms_link.filename, xcms_msms_link.precurMtchID
                           FROM xcms_msms_link
                           LEFT JOIN matches ON matches.pid=s_peak_meta.pid
                           LEFT JOIN c_peak_X_s_peak_meta ON matches.pid=s_peak_meta.pid
                           ')

  XLI <- XLI[XLI$score>=score_f & XLI$match>=frag_nm_f,]

  XLI <- XLI[order(XLI$grpid, -XLI$score),]

  # Only use topn of hits per scan
  XLI <- ddply(XLI, .(pid), function(x){
    if(nrow(x)>topn){
      return(x[topn,])
    }else{
      return(x)
    }
  })

  # Get the most common hit across multiple scans
  xcms_mtch <- ddply(XLI, .(grpid), function(x){

    # Get all the unique names matched from each scan
    pidl <- dlply(x, .(pid), function(y){unique(y$name)})

    # Get the number of scans for this feature
    pids <- unique(x$pid)

    snm <- length(pids)

    # get the 'best' match based on overlap between scans
    overlap_mtch <- overlap(pidl)

    # If joint match just use the first match
    if(length(overlap_mtch>1)){
      overlap_mtch <- overlap_mtch[1]
    }

    # Get the average match scores
    mx <- x[x$name==overlap_mtch,]

    mperc <- median(mx$perc_mtch)
    mscore <- median(mx$score)
    mmatch <- median(mx$match)

    # Get all the matches
    allmtch <- paste(unlist(pidl),collapse=", ")
    allpids <- paste(unlist(pids),collapse=", ")
    alllids <- paste(unlist(unique(mx$lid)),collapse=", ")
    out_v <- c(overlap_mtch, mscore, mmatch, mperc, allmtch, snm, allpids, alllids)

    return(out_v)
  })

  if(nrow(xcms_mtch)==0){
    print('NO MATCHES FOR XCMS')
    dbDisconnect(conT)
    return(0)
  }


  colnames(xcms_mtch) <- c('grpid', 'overlap_match_name', 'median_score',
                           'median_matches', 'median_perc', 'all_matches','number_of_scans','all_pids','all_lids')


  dbWriteTable(conT, name='xcms_match', value=xcms_mtch, row.names=F, append=T)

  return(xcms_mtch)


}


overlap <- function(l) {
  # Stolen from
  # http://codereview.stackexchange.com/questions/17905/compute-intersections-of-all-combinations-of-vectors-in-a-list-of-vectors-in-r
  results <- lapply(l, unique)

  # combinations of m elements of list l
  for (m in seq(along=l)[-1]) {
    # generate and iterate through combinations of length m
    for (indices in combn(seq(length(l)), m, simplify=FALSE)) {

      # make name by concatenating the names of the elements
      # of l that we're intersecting
      name_1 <- paste(names(l)[indices[-m]], collapse="_")
      name_2 <- names(l)[indices[m]]
      name <- paste(name_1, name_2, sep="_")

      results[[name]] <- intersect(results[[name_1]], results[[name_2]])

    }
  }

  for(i in length(results):1){
    if(identical(results[[i]], character(0))){
      next
    }else{
      return(results[[i]])
    }
  }
}



match_targets <- function(target_peaks_list, library_spectra, ppm_tol_MSMS=100, ra_diff=10, library_precs, ppm_tol_prec=50){

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
  library_meta_ids <- unique(library_spectra[,'meta_info_uid'])

  # Get ppm difference between precursor from library and target
  ppmdiff_prec <- as.vector(abs(1e6*outer(target_prec, library_precs, '-')/target_prec))

  # Filter the library spectra that does not meet the precursor tolerance check
  library_spectra_red <- library_spectra[library_spectra[,'meta_info_uid'] %in% library_meta_ids[ppmdiff_prec<=ppm_tol_prec],]

  # Exit if no peaks left (if vector it means it is just 1 row)
  if(is.vector(library_spectra_red)){
    lnms <- names(library_spectra_red)
    library_spectra_red <- matrix(library_spectra_red, nrow=1, ncol=length(library_spectra_red))
    colnames(library_spectra_red) <- lnms
  }else if (nrow(library_spectra_red)==0){
    return(NULL)
  }

  # calculate ppm difference of all peaks (I am presuming this is faster than doing it library spectra at a time...)
  ppmdiff_all <- abs(1e6*outer(target_peaks[,'mz'], library_spectra_red[,'mz'], '-')/library_spectra_red[,'mz'])

  # Calculate the change in percentage of RA of everything
  idiff_all <- abs(outer(target_peaks[,'ra'], library_spectra_red[,'ra'], '-')/library_spectra_red[,'ra'])*100

  # Get index for difference matrices
  gidx <- c(0,cumsum(table(library_spectra_red[,'meta_info_uid'])))

  lnm <- length(unique(library_spectra_red[,'meta_info_uid']))


  l = list()

  library_meta_info_ids <- unique(library_spectra_red[,'meta_info_uid'])

  for(j in 1:lnm){

    uid <- library_meta_info_ids[j]

    mtch<- matchi(ppmdiff = ppmdiff_all[,(gidx[j]+1):gidx[j+1]],
                  idiff = idiff_all[,(gidx[j]+1):gidx[j+1]],
                  library_peaks = library_spectra_red[library_spectra_red[,'meta_info_uid']==uid,],
                  target_peaks = target_peaks,
                  ppm_tol_MSMS=ppm_tol_MSMS,
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


matchi <-function(library_peaks, target_peaks, ppmdiff, idiff, ppm_tol_MSMS=50, ra_diff=10){

  if (!any(ppmdiff<ppm_tol_MSMS)){
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

  # Follwoing the pMatch-hammer method for peak matching but with slight variation that we also check the percentage difference for
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

    if (bm>ppm_tol_MSMS){
      allpeaks[[i]] <- rbind(target_peaks[i,],c(0, 0, 2,0))
      mcount <- append(mcount, 'miss')
      next
    }

    # First check to see if there is a matching intensity value within ra_diff (default 10%)
    intenc <- iD[ppmD<ppm_tol_MSMS & iD<ra_diff & !is.na(ppmD) & !is.na(ra_diff)]

    if (!identical(unname(intenc), numeric(0))){
      matchi <- match(min(intenc, na.rm = T), iD)
    }else{
      matchi <- match(min(ppmD[ppmD<ppm_tol_MSMS], na.rm=T), ppmD)
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


group1d <- function(v, thr){

  v <- sort(v)
  a <- ''
  cl <-  1
  c <- 1

  for(i in 1:(length(v)-1)){
    a[i] <- (1e6*(v[i+1]-v[i]))/v[i]

    if(as.numeric(a[i])<thr){
      cl[i+1] <- c
    }else{
      c <- c+1
      cl[i+1] <- c
    }

  }
  return(cl)
}


alignPeaks <- function(v, type, thr){

  v <- sort(v)
  a <- ''
  p <- ''
  c <- 1

  for(i in 1:(length(v)-1)){
    if (!ty[i]==ty[i+1]){
      a[i] <- (1e6*(v[i+1]-v[i]))/v[i]
      if (a[i]>thr){
        p[i] = p
        p[i+1] = p
      }else{

      }
      print(paste(ty[i], ty[i+1], a[i]))
    }else{
      print(i)
    }


  }
  return(cl)
}


CosSim <- function(A,B) {
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}





