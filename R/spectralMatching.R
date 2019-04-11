#' @title Spectral matching
#'
#' @description
#' Perform spectral matching to spectral libraries using dot product cosine on a LC-MS/MS dataset and link to XCMS features.
#'
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
#' pa <- filterFragSpectra(pa, allfrag=TRUE)
#' pa <- averageAllFragSpectra(pa)
#' q_dbPth <- createDatabase(pa, xset)
#' result <- spectralMatching(q_dbPth)
#' @importFrom magrittr %>%
#' @name %>%
#' @export
spectralMatching <- function(
                             q_dbPth,
                             l_dbPth=NA,

                             q_purity=NA,
                             q_ppmProd=10,
                             q_ppmPrec=5,
                             q_raThres=0,
                             q_pol='positive',
                             q_instrumentTypes=NA,
                             q_instrument=NA,
                             q_sources=NA,
                             q_spectraType="av_all",
                             q_pids=NA,
                             q_rtrange=c(NA, NA),
                             q_spectraFilter=TRUE,
                             q_xcmsGroups=NA,

                             l_purity=NA,
                             l_ppmProd=10,
                             l_ppmPrec=5,
                             l_raThres=0,
                             l_pol='positive',
                             l_instrumentTypes=NA,
                             l_instrument=NA,
                             l_sources=NA,
                             l_spectraType=NA,
                             l_pids=NA,
                             l_rtrange=c(NA, NA),
                             l_spectraFilter=NA,
                             l_xcmsGroups=NA,


                             usePrecursors=TRUE,
                             score_thres=0.6,
                             raW=0.5,
                             mzW=2,
                             rttol=NA,
                             match_alg='dpc',
                             cores=1,
                             out_dir='.',
                             createNew=FALSE){
  message("Running msPurity spectral matching function for LC-MS(/MS) data")
  if(createNew){
    nfn <- file.path(out_dir,paste('lcmsms_data', format(Sys.time(), "%Y-%m-%d-%I%M%S"), 'spectral_matching_result.sqlite', sep="-"))
    file.copy(from = q_dbPth, to=nfn)
    q_dbPth <- nfn
  }


  if (is.na(l_dbPth)){
    l_dbPth <- system.file("extdata", "library_spectra", "library_spectra.db", package="msPurityData")
  }

  ########################################################
  # Filter the query dataset
  ########################################################
  message("Filter query dataset")
  q_con <- DBI::dbConnect(RSQLite::SQLite(), q_dbPth)

  q_speakmeta <- filterSMeta(purity =q_purity,
              raThres= q_raThres,
              pol = q_pol,
              instrumentTypes = q_instrumentTypes,
              instrument = q_instrument,
              sources = q_sources,
              pids = q_pids,
              rtrange = q_rtrange,
              con = q_con,
              xcmsGroups = q_xcmsGroups,
              spectraType = q_spectraType)

  q_speaks <- getScanPeaksSqlite(q_con, q_spectraFilter)

  ########################################################
  # Filter the library dataset
  ########################################################
  message("Filter library dataset")
  l_con <- DBI::dbConnect(RSQLite::SQLite(), l_dbPth)

  l_speakmeta <- filterSMeta(purity = l_purity,
                             raThres = l_raThres,
                             pol = l_pol,
                             instrumentTypes = l_instrumentTypes,
                             instrument = l_instrument,
                             sources = l_sources,
                             pids = l_pids,
                             rtrange = l_rtrange,
                             con = l_con,
                             xcmsGroups = l_xcmsGroups,
                             spectraType = l_spectraType)

  l_speaks <- getScanPeaksSqlite(l_con, l_spectraFilter)

  ########################################################
  # Loop through the query dataset and spectra match
  # against the library spectra
  ########################################################
  # can't parallize dpylr without non cran package
  # Go back to using good old plyr
  q_fpids <- dplyr::pull(q_speakmeta, pid)


  if (cores>1){
    cl<-parallel::makeCluster(cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    parallel = TRUE
  }else{
    parallel = FALSE
  }
  message('align and match')

  matched <- plyr::adply(q_fpids, 1, queryVlibrary, .parallel = parallel,
                                      q_speakmeta=q_speakmeta,
                                      q_speaks=q_speaks,
                                      l_speakmeta=l_speakmeta,
                                      l_speaks=l_speaks,
                                      q_ppmPrec=q_ppmPrec,
                                      q_ppmProd=q_ppmProd,
                                      l_ppmPrec=l_ppmPrec,
                                      l_ppmProd=l_ppmProd,
                         .progress="text")

  if (nrow(matched)==0){
    message('No matches found')
    return(NULL)
  }

  print(matched)

  matched$mid <- 1:nrow(matched)
  custom_dbWriteTable(name_pk = 'mid', df=matched, fks=NA, table_name = 'sm_matches', con = q_con)

  ########################################################
  # Create a summary table for xcms grouped objects
  ########################################################
  if (DBI::dbExistsTable(q_con, "c_peak_groups")){
    # check if the query is from an XCMS object
    message("Summarising LC features annotations")

    xcmsMatchedResults <- getXcmsSmSummary(q_con, matched,spectraType=q_spectraType)
  }else{
    xcmsMatchedResults <- NA
  }

  DBI::dbDisconnect(q_con)
  DBI::dbDisconnect(l_con)


  return(list('q_dbPth' = q_dbPth, 'xcmsMatchedResults' = xcmsMatchedResults))
}

getScanPeaksSqlite <- function(con, spectraFilter, spectraType){
  if (DBI::dbExistsTable(con, "s_peaks")){
    speaks <- con %>% dplyr::tbl("s_peaks")
  }else if (DBI::dbExistsTable(con, "library_spectra")) {
    # old sqlite format
    speaks <- con %>% dplyr::tbl("library_spectra")
  }else{
    stop('No spectra available')
  }

  if (!is.na(spectraFilter)){
    return(speaks %>% dplyr::filter(pass_flag==TRUE))
  }else{
    return(speaks)
  }

  if (spectraType){
    return(speaks %>% dplyr::filter(type==spectraType))
  }else{
    return(speaks)
  }



}

getXcmsSmSummary <- function(con, matched, scoreF=0, fragNmF=1, spectraType='scans'){

  if (spectraType=='scans'){
    sqlStmt <- sprintf("SELECT * FROM c_peak_groups
                           LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.grpid=c_peak_groups.grpid
                           LEFT JOIN c_peaks on c_peaks.cid=cXg.cid
                           LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.cid=c_peaks.cid
                           LEFT JOIN s_peak_meta ON cXs.pid=s_peak_meta.pid
                           WHERE s_peak_meta.pid in (%s)", paste(unique(matched$qpid), collapse=','))

  }else{
    sqlStmt <- sprintf("SELECT cpg.*, spm.pid FROM c_peak_groups AS cpg
                           LEFT JOIN s_peak_meta AS spm ON cpg.grpid=spm.grpid
                           WHERE spm.pid in (%s)", paste(unique(matched$qpid), collapse=','))

  }
  xcmsGroupedPeaks <- DBI::dbGetQuery(con, sqlStmt)
  xcmsMatchedResults <- merge(xcmsGroupedPeaks, matched, by.x='pid', by.y='qpid')
  if(nrow(xcmsMatchedResults)==0){
    message('NO MATCHES FOR XCMS')
    DBI::dbDisconnect(q_con)
    return(0)
  }

  xcmsMatchedResults <- xcmsMatchedResults[order(xcmsMatchedResults$grpid, -as.numeric(xcmsMatchedResults$dpc)),]


  DBI::dbWriteTable(con, name='xcms_match', value=xcmsMatchedResults, row.names=F, append=T)

  return(xcmsMatchedResults)


}

filterSMeta <- function(purity=NA,
                        raThres=0,
                        pol='positive',
                        instrumentTypes=NA,
                        instrument=NA,
                        sources=NA,
                        xcmsGroups=NA,
                        pids=NA,
                        rtrange=c(NA, NA),
                        spectraType=NA,
                        con){

  if (DBI::dbExistsTable(con, "s_peak_meta")){
    speakmeta <- con %>% dplyr::tbl("s_peak_meta")
  }else if (DBI::dbExistsTable(con, "library_spectra_meta")) {
    # old sqlite format
    speakmeta <- con %>% dplyr::tbl("library_spectra_meta")
  }else{
    stop('No meta data for spectra available')
  }

  if (!is.na(pol)){
    speakmeta <- speakmeta %>% dplyr::filter(lower(polarity) == lower(pol))
  }

  if (!anyNA(instrumentTypes)  && !is.na(instrumentTypes)){
    speakmeta <- speakmeta %>% dplyr::filter(instrument_type %in% instrumentTypes || instrument_type %in% instrumentTypes)
  }else if (!anyNA(instrumentTypes)){
    speakmeta <- speakmeta %>% dplyr::filter(instrument_type %in% instrumentTypes)
  }else if (!anyNA(instrument)){
    speakmeta <- speakmeta %>% dplyr::filter(instrument %in% instrument)
  }

  if(!anyNA(sources) & DBI::dbExistsTable(con, "library_spectra_source")){
    sourcetbl <- con %>% dplyr::tbl("library_spectra_source")
    speakmeta <- speakmeta %>% dplyr::left_join(sourcetbl,  by=c('library_spectra_source_id'='id'),suffix = c("", ".y")) %>%
              dplyr::filter(name.y %in% sources)
  }

  if(!is.na(purity)){
    speakmeta <- speakmeta %>% dplyr::filter(inPurity > purity)
  }
  if(!anyNA(pids)){
    if ("pid" %in% names(speakmeta %>% dplyr::collect())){
      speakmeta <- speakmeta %>% dplyr::filter(pid %in% pids)
    }else{
      speakmeta <- speakmeta %>% dplyr::filter(id %in% pids)
    }

  }

  if(!anyNA(rtrange)){
    speakmeta <- speakmeta %>% dplyr::filter(retention_time > rtrange[1] && retention_time < rtrange[2])
  }

  if (!anyNA(xcmsGroups) & DBI::dbExistsTable(con, "c_peak_groups")){

    XLI <- DBI::dbGetQuery(con, paste0('SELECT cpg.grpid, spm.pid FROM s_peak_meta AS spm
                           LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.pid=spm.pid
                           LEFT JOIN c_peaks AS cp on cp.cid=cXs.cid
                           LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.cid=cp.cid
                           LEFT JOIN c_peak_groups AS cpg ON cXg.grpid=cpg.grpid
                           WHERE cpg.grpid in (', paste(xcmsGroups, collapse = ','), ')'))

    speakmeta <- speakmeta %>% dplyr::filter(grpid %in% xcmsGroups || pid %in% XLI$pid)

  }

  if (!is.na(spectraType)){
    if (spectraType=='av_all'){
      spectraType = 'all'
    }
    speakmeta <- speakmeta %>% dplyr::filter(spectrum_type == spectraType)
  }



  return(speakmeta)


}

filterPrecursors <- function(q_pid, q_speakmeta, l_speakmeta, q_ppmPrec, l_ppmPrec){

  return(l_speakmetaFiltered)

}



queryVlibrary <- function(q_pid, q_speakmeta, q_speaks, l_speakmeta, l_speaks,
                          q_ppmPrec, q_ppmProd, l_ppmPrec, l_ppmProd, usePrecursors){

  # filter precursors
  # get meta and peaks assoicated with pid

  q_speakmetai <- q_speakmeta %>% dplyr::filter(pid==q_pid) %>% dplyr::collect()
  q_precMZ <- q_speakmetai$precursor_mz
  q_speaksi <- q_speaks %>% dplyr::filter(pid==q_pid)  %>% dplyr::collect()


  if (nrow(q_speaksi)==0){
    return(NULL)
  }

  # Check if precursors are within tolerance
  # We have ppm tolerances for both the library and the query
  q_precMZlo = q_precMZ - ((q_precMZ*0.000001)*q_ppmPrec)
  q_precMZhi = q_precMZ + ((q_precMZ*0.000001)*q_ppmPrec)

  # Search against the range for the library
  l_fspeakmeta <- l_speakmeta %>%
     dplyr::filter( (q_precMZhi >= precursor_mz - ((precursor_mz*0.000001)*l_ppmPrec))
                    &
                    (precursor_mz + ((precursor_mz*0.000001)*l_ppmPrec) >= q_precMZlo)) %>%
                    #summarise(id)  %>% # need to change pid
                    dplyr::collect()

  if(nrow(l_fspeakmeta)==0){
    return(NULL)
  }

  if ('pid' %in% colnames(l_fspeakmeta )){
    l_fpids <- l_fspeakmeta$pid

  }else{
    l_fpids <- l_fspeakmeta$id
  }


  searched <- plyr::adply(l_fpids , 1, queryVlibrarySingle,
                            q_speaksi=q_speaksi,
                            l_speakmeta=l_speakmeta,
                            l_speaks=l_speaks,
                            q_ppmProd=q_ppmProd,
                            l_ppmProd=l_ppmProd

              )

  colnames(searched)[1] <- 'lpid'
  searched$qpid <- q_pid





  return(searched)



}


queryVlibrarySingle <- function(l_pid, q_speaksi, l_speakmeta, l_speaks, q_ppmProd, l_ppmProd){

  if ('pid' %in% colnames(l_speaks)){
    l_speaksi <- l_speaks %>% dplyr::filter(pid==l_pid) %>% dplyr::collect() # need to change to pid
    l_speakmetai <- data.frame(l_speakmeta %>% dplyr::filter(pid==l_pid) %>% dplyr::collect())
  }else{
    l_speaksi <- l_speaks %>% dplyr::filter(library_spectra_meta_id==l_pid) %>% dplyr::collect()
    l_speakmetai <- data.frame(l_speakmeta %>% dplyr::filter(id==l_pid) %>% dplyr::collect())
  }

  # ensure we have the relative abundance
  l_speaksi$ra <- (l_speaksi$i/max(l_speaksi$i))*100
  q_speaksi$ra <- (q_speaksi$i/max(q_speaksi$i))*100

  am <- alignAndMatch(q_speaksi, l_speaksi, q_ppmProd, l_ppmProd)

  return(c(am, 'accession'=l_speakmetai$accession,
               'name'=l_speakmetai$name))
}


alignAndMatch <- function(q_speaksi, l_speaksi, q_ppmProd, l_ppmProd){
  ###################
  # Align
  ###################

  q_speaksi$w <- (q_speaksi$mz^2) * (q_speaksi$ra^0.5)
  l_speaksi$w <- (l_speaksi$mz^2) * (l_speaksi$ra^0.5)
  q_speaksi <- data.frame(q_speaksi)
  l_speaksi <- data.frame(l_speaksi)

  aligned <-align(q_speaksi, l_speaksi, l_ppmProd, q_ppmProd, raDiffThres=10)

  ##################
  # Match
  ##################
  dpcOut <- dpc(aligned$q, aligned$l)

  rl <- aligned$l[!aligned$q==0]
  rq <- aligned$q[!aligned$q==0]
  rdpcOut <- dpc(rq, rl)

  cdpcOut <- cdpc(aligned$q, aligned$l)

  return(c('dpc'=dpcOut, 'rdpc'=rdpcOut, 'cdpc'=cdpcOut,
           'mcount'=aligned$mcount, 'allcount'=aligned$total,
           'mpercent'=aligned$percMatchAll))

}





mzComparePPMrange <- function(mz1, mz2, ppm1, ppm2){
  mz1Lo = round(mz1 - ((mz1*0.000001)*ppm1), 10)
  mz1Up = round(mz1 + ((mz1*0.000001)*ppm1), 10)

  # get ranges of the "blank" to remove
  mz2Lo = round(mz2 - ((mz2*0.000001)*ppm2), 10)
  mz2Up = round(mz2 + ((mz2*0.000001)*ppm2), 10)

  return(overlap(mz1Lo, mz1Up, mz2Lo, mz2Up))
}


align <-function(q_speaksi,l_speaksi, l_ppmProd=100, q_ppmProd=100, raDiffThres=10){

  # Following the pMatch-hammer method for peak matching but with slight variation that we also check the percentage difference for
  # the relative intensity as well
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3928522/

  # get query and library peaks in the correct format
  q_speaksi$ql <-'query'
  l_speaksi$ql <-'library'

  q_speaksi <- q_speaksi[,c('mz','ra','ql','w')]
  l_speaksi <- l_speaksi[,c('mz','ra','ql','w')]

  # Get matrix of ppm bool, ppm differences and relative abundance differences
  idiff <- ppmdiff <- ppmbool <- matrix(0L, nrow = nrow(q_speaksi), ncol = nrow(l_speaksi))

  for(x in 1:nrow(q_speaksi)){
    for(y in 1:nrow(l_speaksi)){
      # Get ranges of the "sample"
      idiff[x, y] = abs(q_speaksi$ra[x] - l_speaksi$ra[y])   ## need to update to RA!!!
      ppmbool[x, y] = mzComparePPMrange(q_speaksi$mz[x], l_speaksi$mz[y], q_ppmProd, l_ppmProd)
      ppmdiff[x, y] = abs(1e6*(q_speaksi$mz[x]-abs(l_speaksi$mz[y]))/l_speaksi$mz[y])
    }
  }

  # Loop through the differences for the target to library
  allpeaks <-  list()
  mcount = NULL

  # we keep a vector detailing which library peaks we have already matched
  lp_remain <- 1:ncol(ppmbool)

  for (i in 1:nrow(ppmbool)){
    # get ppm diff and intensity diff
    ppmB <- ppmbool[i,]
    ppmD <- ppmdiff[i,]
    iD <- idiff[i,]

    if(sum(ppmB, na.rm = TRUE)==0){
      allpeaks[[i]] <- rbind(q_speaksi[i,],c(0, 0, 'library',0))
      mcount <- append(mcount, 'miss')
      next
    }

    # First check to see if there is a matching intensity value within ra_diff (default 10%)
    intenc <- iD[ppmD==1 & iD<raDiffThres & !is.na(ppmD) & !is.na(raDiffThres)]

    if (!identical(unname(intenc), numeric(0))){
      matchi <- match(min(intenc, na.rm = TRUE), iD)
    }else{
      # if there are identical matches for the ppm range searched. Then pick the match
      # with the smallest ppm difference. This doesn't happen that often, but to
      # keep consistent when it does
      matchi <- match(min(ppmD[ppmB==1], na.rm=TRUE), ppmD)
    }

    allpeaks[[i]] <- rbind(q_speaksi[i,],l_speaksi[matchi,])
    mcount <- append(mcount, 'match')
    lp_remain[matchi] <- NA
    idiff[,matchi] <- NA
    ppmdiff[,matchi] <- NA
    ppmbool[,matchi] <- NA
  }

  cp = length(allpeaks)+1

  # get ramining query empty values
  for (i in lp_remain){
    if(is.na(i)){
      next
    }
    allpeaks[[cp]] <- rbind(c(0, 0, 'query',0), l_speaksi[i,])
    mcount <- append(mcount, 'miss')
    cp = cp+1
  }

  allpeaksm <- data.frame(do.call(rbind, allpeaks))

  wqm <- allpeaksm[allpeaksm[,'ql']=='query',]
  if (is.vector(wqm)){
    wq <- as.numeric(wqm['w'])
  }else{
    wq <- as.numeric(wqm[,'w'])
  }

  wlm <- allpeaksm[allpeaksm[,'ql']=='library',]

  if (is.vector(wlm)){
    wl <- as.numeric(wlm['w'])
  }else{
    wl <- as.numeric(wlm[,'w'])
  }


  pmatch <- sum(mcount=='match')/length(mcount)
  return(list(q=wq, l=wl, mcount=sum(mcount=='match'), total=length(mcount), percMatchAll=pmatch))
}









