#' @title Spectral matching
#'
#' @description
#' Perform spectral matching to spectral libraries using dot product cosine on a LC-MS/MS dataset and link to XCMS features.
#'
#' @param query_db_pth character; Path of the database of targets (queries) that will be searched against the library spectra. Generated
#'                       either from frag4feature or from create_database functions.
#' @param library_db_pth character [optional]; path to library spectral SQLite database.
#'                                             Defaults to msPurityData package data.
#' @param ra_thres_q numeric; Relative abundance threshold for target (query) spectra
#'                      (Peaks below this RA threshold will be excluded)
#' @param ra_thres_l numeric; Relative abundance threshold for library spectra
#' @param cores numeric; Number of cores to use
#' @param pol character; Polarity ['positive' or 'negative']
#' @param ppm_tol_prod numeric; PPM tolerance to match to product
#' @param ppm_tol_prec numeric; PPM tolerance to match to precursor
#' @param score_thres numeric; Dot product cosine score threshold
#' @param ra_w numeric; Relative abundance weight for spectra
#' @param mz_w numeric; mz weight for spectra
#' @param match_alg character; Can either use dot product cosine (dpc) or match factor (mf) for spectral matching. Defaults to dpc
#' @param spectra_type_q character; Type of fragmentation spectra from query to match with "scans" =  all individual scans,
#'                                  "av_intra" = averaged spectra (intra), "av_inter" = averaged spectra (inter), "av_all" = averaged all
#'                                   spectra ignoring inter-intra relationships
#' @param topn numeric [optional]; Only use top n matches
#' @param db_name character [optional]; Name of the result database
#'                                   (e.g. can use CAMERA peaklist)
#' @param instrument_types vector [optional]; Vector of instrument types, defaults to all
#' @param library_sources vector [optional]; Vector of library sources. Default option is for massbank only but the 'lipidblast'
#'                                    library is also available
#' @param scan_ids vector [optional]; Vector of unique scan ids calculated from msPurity "pid". These scans will on
#'                        used for the spectral matching. All scans will be used if set to NA
#' @param rt_range vector [optional]; Vector of rention time range to filter the library spectra (rtmin, rtmax). Default is to ignore
#'                                    retention time range
#' @param rttol numeric [optional]; Tolerance in time range between the Library and Query database retention time (in seconds) NA to ignore
#' @param pa purityA object [optional]; If target_db_pth set to NA, a new database can be created using pa, xset and grp_peaklist
#' @param xset xcms object [optional]; If target_db_pth set to NA, a new database can be created using pa, xset and grp_peaklist
#' @param grp_peaklist dataframe [optional]; If target_db_pth set to NA, a new database can be created using pa, xset and grp_peaklist
#' @param out_dir character [optional]; If target_db_pth set to NA, Out directory for the SQLite result database
#' @param target_db_pth character [deprecated]; The query database path (use query_db_pth for future use)
#' @param ra_thres_t numeric [deprecated]; The relative abundance threshold for the query spectra (use ra_thres_q for future use)
#'
#' @return list of database details and dataframe summarising the results for the xcms features
#' @examples
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- frag4feature(pa, xset)
#' pa <- averageAllFragSpectra(pa)
#' db_path <- create_database(pa, xset)
#' result <- spectral_matching(db_path, spectra_type_q="av_all")
#' @export
spectral_matching <- function(query_db_pth, ra_thres_l=0, ra_thres_q=2, cores=1, pol='positive', ppm_tol_prod=10, ppm_tol_prec=5,
                                     score_thres=0.6, topn=NA,  db_name=NA, library_db_pth=NA,
                                     instrument_types=NA, library_sources='massbank', scan_ids=NA,
                                     pa=NA, xset=NA, grp_peaklist=NA, out_dir='.', ra_w=0.5, mz_w=2,
                                     spectra_type_q="scans", ra_thres_t=NA, target_db_pth=NA, rt_range=c(NA, NA), rttol=NA,
                                     match_alg='dpc'){
  message("Running msPurity spectral matching function for LC-MS(/MS) data")

  if (!is.na(ra_thres_t)){
    message("ra_thres_t argument has been deprecated and will be remove in future versions of msPurity,
            use ra_thres_q for future use")
    ra_thres_q <- ra_thres_t
  }

  if (!is.na(target_db_pth)){
    message("Please note that target_db_pth argument has been deprecated and will be remove in future versions of msPurity,
            use query_db_pth for future use")
    query_db_pth <- target_db_pth
  }



  ########################################################
  # Export the target data into sqlite database
  ########################################################
  if (is.na(query_db_pth)){
    query_db_pth <- create_database(pa=pa, xset=xset, out_dir=out_dir, grp_peaklist=grp_peaklist, db_name=db_name)
  }

  if (is.na(library_db_pth)){
    library_db_pth <- system.file("extdata", "library_spectra", "library_spectra.db", package="msPurityData")
  }


  ########################################################
  # Perform the spectral matching
  ########################################################
  message("Performing spectral matching")
  matched <- match_2_library(query_db_pth,
                  library_db_pth,
                  ra_thres_q=ra_thres_q,
                  ra_thres_l=ra_thres_l,
                  cores=cores,
                  pol=pol,
                  ppm_tol_prod=ppm_tol_prod,
                  ppm_tol_prec=ppm_tol_prec,
                  topn=topn,
                  instrument_types = instrument_types,
                  library_sources = library_sources,
                  scan_ids = scan_ids,
                  ra_w = ra_w,
                  mz_w = mz_w,
                  spectra_type_q = spectra_type_q,
                  rt_range = rt_range,
                  rttol = rttol,
                  match_alg=match_alg

                  )

  ########################################################
  # Create a summary table for xcms grouped objects
  ########################################################
  if (matched){
    message("Summarising LC features annotations")
    xcms_summary_df <- get_xcms_sm_summary(query_db_pth, topn=topn, score_f=score_thres,spectra_type_q=spectra_type_q)
  }else{
    xcms_summary_df <- NA
  }


  return(list('result_db_pth' = query_db_pth, 'xcms_summary_df' = xcms_summary_df))
}



ra_calc <- function(x){
  x$ra <- (x$i/max(x$i))*100
  return(x)
}

get_query_spectra_list_s_peak <- function(x, scan_info){
  si = scan_info[scan_info$pid==unique(x$pid),]
  list(x, si$precursorMZ, si$retentionTime)
}

get_query_spectra_list_c_peak_group <- function(x, c_peak_groups){
  av_peaks = x[,c('mz','ra','type','w','grpid','grpid_fileid')]
  cgi <- c_peak_groups[c_peak_groups$grpid==unique(x$grpid),]
  list(av_peaks, cgi$mz, cgi$rt)
}

match_2_library <- function(query_db_pth, library_db_pth, instrument_types=NA, mslevel=NA, mslevel_match=TRUE,
                            ra_thres_q=2, ra_thres_l=0, cores=1, pol, ppm_tol_prod=100, ppm_tol_prec=50, topn=5,
                            library_sources=NA, scan_ids=NA, ra_w=0.5, mz_w=2,
                            spectra_type_q="scans", rt_range=NA, rttol=NA, match_alg='dpc'){




  ########################################################
  # Get target (query) spectra
  ########################################################
  # Get all of the target MS2 data
  conQ <- DBI::dbConnect(RSQLite::SQLite(), query_db_pth)

  # Get scan peaks for target spectra
  if (spectra_type_q=="scans"){
    query_spectra <- DBI::dbGetQuery(conQ, 'SELECT * FROM s_peaks ' )
  }else if (spectra_type_q=="av_intra"){
    query_spectra <- DBI::dbGetQuery(conQ, 'SELECT * FROM av_peaks WHERE  method="intra" AND pass_flag=1' )
  }else if (spectra_type_q=="av_inter"){
    query_spectra <- DBI::dbGetQuery(conQ, 'SELECT * FROM av_peaks WHERE  method="inter" AND pass_flag=1' )
  }else if (spectra_type_q=="av_all"){
    query_spectra <- DBI::dbGetQuery(conQ, 'SELECT * FROM av_peaks WHERE  method="all" AND pass_flag=1' )
  }

  if (!anyNA(scan_ids) && spectra_type_q=="scans"){
    query_spectra <- query_spectra[query_spectra[,'pid'] %in% scan_ids,]
  }

  if (spectra_type_q=="scans"){
    query_spectra$grpid_fileid <- paste(NA, query_spectra$fileid, sep='_')
  }else{
    query_spectra$grpid_fileid <- paste(query_spectra$grpid, query_spectra$fileid, sep='_')
  }



  # Calculate relative abundance
  if (spectra_type_q=="scans"){
    query_spectra <- plyr::ddply(query_spectra, ~ pid, ra_calc)
  }else if (spectra_type_q=="av_intra"){

    query_spectra <- plyr::ddply(query_spectra, ~ grpid_fileid, ra_calc)
  }else{
    query_spectra <- plyr::ddply(query_spectra, ~ grpid, ra_calc)
  }


  # Relative abundance threshold
  query_spectra <- query_spectra[query_spectra$ra>=ra_thres_q,]

  # Flag for spectra type
  query_spectra$type <- 1


  ########################################################
  # Get library spectra
  ########################################################
  # Get all of the library peaks based on parameters given
  conL <- DBI::dbConnect(RSQLite::SQLite(), library_db_pth)

  # Match each target spectra to the library spectra
  # Only keep matches with certain criteria certain criteria

  if (anyNA(library_sources)){
    library_meta_query <- sprintf("SELECT * FROM library_spectra_meta WHERE lower(polarity) = lower('%s')", pol)
  }else{
    l_source_str <- paste("'",paste(library_sources, collapse = "', '"), "'", sep='')
    library_meta_query <- sprintf("SELECT m.*, s.name AS source_name FROM library_spectra_meta as m
                                    LEFT JOIN library_spectra_source AS s ON s.id=m.library_spectra_source_id
                                    WHERE lower(m.polarity) = lower('%s') AND s.name IN (%s)", pol, l_source_str)

  }




  if (!anyNA(instrument_types)){
    # instrument_types <- c('CE-ESI-TOF', 'ESI-ITFT', 'ESI-ITTOF', 'ESI-QTOF', 'LC-ESI-IT',
    #                       'LC-ESI-ITFT', 'LC-ESI-ITTOF','LC-ESI-QFT', 'LC-ESI-QIT', 'LC-ESI-QQ', 'LC-ESI-QTOF', 'LC-ESI-TOF')
    ints_string <- paste("'",paste(instrument_types, collapse = "', '"), "'", sep='')
    library_meta_query <- paste(library_meta_query, sprintf(" AND instrument_type IN (%s)", ints_string))

  }


  if (!anyNA(rt_range)){
    library_meta_query <- paste(library_meta_query, sprintf(" AND retention_time > %f AND retention_time < %f", rt_range[1], rt_range[2]))
  }


  library_meta <- DBI::dbGetQuery(conL, library_meta_query)

  if (nrow(library_meta)==0){
    message('No library spectra matching criteria')
    return(NULL)
  }

  if (!is.na(mslevel)){
    library_meta <- library_meta[library_meta$ms_level==mslevel,]
  }

  if (nrow(library_meta)==0){
    message('No library spectra matching criteria')
    return(NULL)
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
  query_spectra$w <- (query_spectra$ra^ra_w)*(query_spectra$mz^mz_w)
  library_spectra$w <- (library_spectra$ra^ra_w)*(library_spectra$mz^mz_w)


  ########################################################
  # Get target (query) spectra in format ready for matching
  ########################################################
  # Order by relative abundance for target spectra
  if (spectra_type_q=="scans"){
    query_spectra <- query_spectra[order(query_spectra[,'pid'], -query_spectra[,'ra']), ]
    scan_info <- DBI::dbGetQuery(conQ,'SELECT * FROM s_peak_meta' )
    # we loop through the target spectra as list
    query_spectra_list <- plyr::dlply(query_spectra, ~ pid, get_query_spectra_list_s_peak,
                                      scan_info=scan_info)



  }else{
    # get the group peaks (but also make sure we have the full min and max peakwidth)
    #c_peak_groups <- DBI::dbGetQuery(conQ,
    #                                      'SELECT  g.*, MIN(c.rtmin) as rtmin_full,  MAX(c.rtmax) as rtmax_full FROM c_peak_groups as g
    #                                          LEFT JOIN c_peak_X_c_peak_group as cxg ON
    #                                            cxg.grpid=g.grpid
    #                                          LEFT JOIN c_peaks as c ON
    #                                            cxg.cid=c.cid
    #                                        WHERE cxg.bestpeak=1
    #                                        GROUP BY cxg.grpid;')
    c_peak_groups <- DBI::dbGetQuery(conQ, 'SELECT * FROM c_peak_groups')

    # we loop through the target spectra as list
    if (spectra_type_q=='av_intra'){

      query_spectra <- query_spectra[order(query_spectra[,'grpid_fileid'], -query_spectra[,'ra']), ]
      query_spectra_list <- plyr::dlply(query_spectra, ~ grpid_fileid, get_query_spectra_list_c_peak_group,
                                        c_peak_groups=c_peak_groups)
    }else{
      query_spectra <- query_spectra[order(query_spectra[,'grpid'], -query_spectra[,'ra']), ]
      query_spectra_list <- plyr::dlply(query_spectra, ~ grpid, get_query_spectra_list_c_peak_group,
                                        c_peak_groups=c_peak_groups)
    }



  }


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

  allmatches_l <- plyr::llply(query_spectra_list, .parallel = parallel, match_targets,
                        library_spectra=library_spectra,
                        library_precs=library_meta$precursor_mz,
                        library_rt=library_meta$retention_time,
                        ppm_tol_prod=ppm_tol_prod,
                        ppm_tol_prec=ppm_tol_prec,
                        rttol=rttol,
                        match_alg=match_alg)



  if (spectra_type_q=="scans"){
    allmatches <- plyr::ldply(allmatches_l, .id = 'pid')
  }else if (spectra_type_q=="av_intra"){
    allmatches <- plyr::ldply(allmatches_l, .id = 'grpid_fileid')
    grpid_fileid_df <- data.frame(stringr::str_split_fixed(allmatches$grpid_fileid, "_", 2))
    colnames(grpid_fileid_df) <- c('grpid', 'fileid')
    allmatches <- cbind(allmatches, grpid_fileid_df)

  }else{
    allmatches <- plyr::ldply(allmatches_l, .id = 'grpid')
  }



  ########################################################
  # Submit to database
  ########################################################
  if (nrow(allmatches)>0){
    if (!is.na(topn)){
      if (spectra_type_q=="scans"){
        allmatches <- plyr::ddply(allmatches, ~ pid, get_topn)
      }else if (spectra_type_q=="av_intra"){
        allmatches <- plyr::ddply(allmatches, ~ grpid_fileid, get_topn)
      }else{
        allmatches <- plyr::ddply(allmatches, ~ grpid, get_topn)
      }
    }

    library_meta_f <- library_meta[library_meta$id %in% allmatches$lid, ]

    colnames(library_meta_f)[which(colnames(library_meta_f)=='id')] = 'lid'

    # get all compound details as well
    compound_details <- DBI::dbGetQuery(conL, sprintf('SELECT  DISTINCT c.* FROM library_spectra_meta AS m
                                               LEFT JOIN metab_compound AS
                                               c on c.inchikey_id=m.inchikey_id
                                               WHERE m.id IN (%s)', paste(unique(allmatches$lid), collapse=",")) )
    compound_details <- data.frame(lapply(compound_details, as.character), stringsAsFactors=FALSE)
    custom_dbWriteTable(name_pk = 'inchikey_id', fks = NA,
                        df=compound_details, table_name = 'metab_compound', con = conQ, pk_type='TEXT')



    fk_l = list('inchikey_id'=list('new_name'='inchikey_id', 'ref_name'='inchikey_id', 'ref_table'='library_spectra_meta'))
    custom_dbWriteTable(name_pk = 'lid', fks = fk_l,
                        df=library_meta_f, table_name = 'library_spectra_meta', con = conQ)

    allmatches$mid <- 1:nrow(allmatches)

    fks_lid <- list('lid'=list('new_name'='lid', 'ref_name'='lid', 'ref_table'='library_spectra_meta'))

    if (spectra_type_q=="scans"){
      fks_q <- list('pid'=list('new_name'='pid', 'ref_name'='pid', 'ref_table'='s_peak_meta'))
    }else{
      fks_q <- list('grpid'=list('new_name'='grpid', 'ref_name'='grpid', 'ref_table'='c_peak_groups'))
    }

    custom_dbWriteTable(name_pk = 'mid', fks = append(fks_lid, fks_q),
                        df=allmatches, table_name = 'matches', con = conQ)






    matched = TRUE
  }else{
    matched = FALSE
  }

  DBI::dbDisconnect(conQ)
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

get_xcms_sm_summary <- function(query_db_pth, topn=NA, score_f=0.3, frag_nm_f=1,spectra_type_q='scans'){

  # Get all of the target fragmentation data
  conQ <- DBI::dbConnect(RSQLite::SQLite(), query_db_pth)

  if (spectra_type_q=="scans"){
    XLI <- DBI::dbGetQuery(conQ, 'SELECT * FROM c_peak_groups
        LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.grpid=c_peak_groups.grpid
        LEFT JOIN c_peaks on c_peaks.cid=cXg.cid
        LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.cid=c_peaks.cid
        LEFT JOIN s_peak_meta ON cXs.pid=s_peak_meta.pid
        LEFT JOIN matches ON matches.pid=s_peak_meta.pid
        LEFT JOIN library_spectra_meta ON matches.lid=library_spectra_meta.lid
        WHERE matches.score IS NOT NULL')
  }else{
    XLI <- DBI::dbGetQuery(conQ, 'SELECT * FROM c_peak_groups
        LEFT JOIN matches ON matches.grpid=c_peak_groups.grpid
        LEFT JOIN library_spectra_meta ON matches.lid=library_spectra_meta.lid
        WHERE matches.score IS NOT NULL')
  }

  XLI <- XLI[XLI$score>=score_f & XLI$match>=frag_nm_f,]

  XLI <- XLI[order(XLI$grpid, -XLI$score),]

  # Only use topn of hits per scan
  if (!is.na(topn)){
    if (spectra_type_q){
      XLI <- plyr::ddply(XLI, ~ pid, check_topn, topn=topn)
    }else{
      XLI <- plyr::ddply(XLI, ~ grpid, check_topn, topn=topn)
    }

  }


  # Summaries the annotation hits
  xcms_mtch <- plyr::ddply(XLI, ~ grpid, get_ann_summary, spectra_type_q=spectra_type_q)

  if(nrow(xcms_mtch)==0){
    message('NO MATCHES FOR XCMS')
    DBI::dbDisconnect(conQ)
    return(0)
  }

  DBI::dbWriteTable(conQ, name='xcms_match', value=xcms_mtch, row.names=F, append=T)
  DBI::dbDisconnect(conQ)
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

get_ann_summary <- function(x, spectra_type_q){
  # get the 'best' match based on the best scored compounds with
  # the same name
  med_results <- plyr::ddply(x, ~ name, median_match_results)

  med_results <- med_results[order(med_results$best_median_score, decreasing = TRUE),]
  colnames(med_results)[1] <- 'best_name'
  best_match <- med_results[1,]

  unique_names <- med_results$name

  # Get all the matches
  allnames <- paste(unlist(x$name),collapse=", ")
  if (spectra_type_q=='scans'){
    allpids <- paste(unlist(x$pid),collapse=", ")
    snm <- length(unique(unlist(allpids)))

  }else{
    allpids <- NA
    snm <- NA
  }

  alllids <- paste(unlist(x$lid),collapse=", ")
  allscores <- paste(unlist(round(x$score,3)),collapse=", ")



  out_v <- c(unlist(best_match), 'all_names'=allnames, 'all_pids'=allpids,
             'all_lids'=alllids, 'all_scores'=allscores, 'numb_unique_scans'=snm, 'spectra_type'= spectra_type_q)

}

match_targets <- function(target_peaks_list, library_spectra, ppm_tol_prod=100, ra_diff=10, library_precs, library_rt, ppm_tol_prec=50, rttol=NA, match_alg){

  # Get target peaks and precursor of target
  target_peaks <- target_peaks_list[[1]]
  target_prec <- target_peaks_list[[2]]
  target_rt <- target_peaks_list[[3]]

  tcnm = c('mz','ra','type','w')

  # Convert target peaks to matrix for speed
  if(is.vector(target_peaks) | nrow(target_peaks)==1){
    target_peaks <- target_peaks[tcnm]
    target_peaks <- matrix(as.numeric(target_peaks), nrow=1, ncol=4)
    colnames(target_peaks) <- tcnm
  }else{
    target_peaks <- as.matrix(target_peaks[,tcnm])
  }

  # get unique ids of library_spectra precursor 'meta' info
  library_meta_ids <- unique(library_spectra[,'library_spectra_meta_id'])

  # Get ppm difference between precursor from library and target
  out_prec <- outer(target_prec, library_precs, '-')*1e6
  colprec <- library_precs[col(out_prec)] # so we can divide by column if library
  ppmdiff_prec <- as.vector(abs(out_prec/colprec))

  # Filter the library spectra that does not meet the precursor tolerance check
  library_spectra_red <- library_spectra[library_spectra[,'library_spectra_meta_id'] %in% library_meta_ids[ppmdiff_prec<=ppm_tol_prec],]

  # get difference from rention time
  if(!is.na(rttol)){

    rt_diff <- abs(target_rt-library_rt)



    # Filter the library spectra that does not meet the precursor tolerance check
    library_spectra_red <- library_spectra_red[library_spectra_red[,'library_spectra_meta_id'] %in% library_meta_ids[rt_diff<=rttol],]
  }

  # Exit if no peaks left (if vector it means it is just 1 row)
  if(is.vector(library_spectra_red)){
    lnms <- names(library_spectra_red)
    library_spectra_red <- matrix(library_spectra_red, nrow=1, ncol=length(library_spectra_red))
    colnames(library_spectra_red) <- lnms
  }else if (nrow(library_spectra_red)==0){
    return(NULL)
  }

  # calculate ppm error of all peaks (I am presuming this is faster than doing it library spectra at a time...)
  out_peaks_mz <- outer(target_peaks[,'mz'], library_spectra_red[,'mz'], '-')*1e6
  colmz <- library_spectra_red[,'mz'][col(out_peaks_mz)] # so we can calculate column wide
  ppmdiff_all <- abs(out_peaks_mz/colmz)

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
                  ra_diff=ra_diff,
                  match_alg=match_alg)

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


matchi <-function(library_peaks, target_peaks, ppmdiff, idiff, ppm_tol_prod=50, ra_diff=10, match_alg){

  if (!any(ppmdiff<ppm_tol_prod)){
    return(c(NA,NA,NA))
  }

  if(is.vector(target_peaks) || nrow(target_peaks)==1){
      ppmdiff <- matrix(ppmdiff, nrow=1)
      idiff <- matrix(idiff, nrow=1)
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

  # need to ensure target_peaks, ppmdiff and idiff are all in mtarix format

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
      matchi <- match(min(intenc, na.rm = TRUE), iD)
    }else{
      matchi <- match(min(ppmD[ppmD<ppm_tol_prod], na.rm=TRUE), ppmD)
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

  wtm <- allpeaksm[allpeaksm[,'type']==1,]
  if (is.vector(wtm)){
    wt <- as.numeric(wtm['w'])
  }else{
    wt <- as.numeric(wtm[,'w'])
  }

  wlm <- allpeaksm[allpeaksm[,'type']==2,]
  if (is.vector(wlm)){
    wl <- as.numeric(wlm['w'])
  }else{
    wl <- as.numeric(wlm[,'w'])
  }
  if (match_alg=='dpc'){
    sim_out <- CosSim(wt, wl)
  }else if (match_alg=='mf'){
    print('MATCH FACTOR')
    sim_out <- match_factor(wt, wl)
  }


  percentage_match <- sum(mcount=='match')/length(mcount)


  return(c(sim_out,percentage_match,sum(mcount=='match')))

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





