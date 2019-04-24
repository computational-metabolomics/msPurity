#' @title Spectral matching for LC-MS/MS datasets
#' @aliases spectralMatching
#' @description
#' **General**
#'
#' Perform spectral matching to spectral libraries for an LC-MS/MS dataset.
#'
#' The spectral matching is performed from a **Query** SQLite spectral-database against a **Library** SQLite spectral-database.
#'
#' The SQLite schema of the spectral database can be detailed Schema details can be found
#' [here](https://bioconductor.org/packages/release/bioc/vignettes/msPurity/inst/doc/msPurity-spectral-datatabase-schema.html).
#'
#' The query spectral-database in most cases should contain be the "unknown" spectra database generated the msPurity
#' function createDatabase as part of a msPurity-XCMS data processing workflow.
#'
#' The library spectral-database in most cases should contain the "known" spectra from either public or user generated resources.
#' The library SQLite database by default contains data from MoNA including Massbank, HMDB, LipidBlast and GNPS.
#' A larger database can be downloaded from [here](https://github.com/computational-metabolomics/msp2db/releases).
#' To create a user generated library SQLite database the following tool can be used to generate a SQLite database
#' from a collection of MSP files: [msp2db](https://github.com/computational-metabolomics/msp2db/releases).
#' It should be noted though, that as long as the schema of the spectral-database is as described here, then any database can be used
#' for either the library or query -  even allowing for the same database to be used.
#'
#' The spectral matching functionality has four main components, spectral filtering, spectral alignment, spectral matching,
#' and summarising the results.
#'
#' Spectral filtering is simply filtering both the library and query spectra to be search against (e.g. choosing
#' the library source, instrument, retention time, precursor PPM tolerance etc).
#'
#' The spectral alignment stage involves aligning the query peaks to the library peaks. The approach used is similar
#' to modified pMatch algorithm described in Zhou et al 2015.
#'
#' The spectral matching of the aligned spectra is performed against a combined intensity and m/z weighted vector - created for both
#' the query and library spectra (wq and wl). See below:
#'
#' \deqn{w=intensity^x * mz^y}
#'
#' Where x and y represent weight factors, defaults to *x*=0.5 and *y*=2 as per MassBank. These can be adjusted by
#' the user though.
#'
#' The aligned weighted vectors are then matched using dot product cosine, reverse dot product cosine and the composite dot product.
#' See below for dot product cosine equation.
#'
#' \deqn{dpc =  wq * wl / \sqrt{\sum wq^2} * \sqrt{\sum wl^2}}
#'
#' See the vigenttes for more details regarding matching algorithms used.
#'
#' **Example LC-MS/MS processing workflow**
#'
#'  * Purity assessments
#'    +  (mzML files) -> purityA -> (pa)
#'  * XCMS processing
#'    +  (mzML files) -> xcms.xcmsSet -> xcms.merge -> xcms.group -> xcms.retcor -> xcms.group -> (xset)
#'  * Fragmentation processing
#'    + (xset, pa) -> frag4feature -> filterFragSpectra -> averageAllFragSpectra -> createDatabase -> **spectralMatching** -> (sqlite spectral database)
#'
#'
#'
#' @param q_dbPth character; Path of the database of queries that will be searched against the library spectra. Generated from createDatabase
#' @param l_dbPth character; path to library spectral SQLite database. Defaults to msPurityData package data.
#'
#' @param q_purity character; Precursor ion purity threshold for the query spectra
#' @param q_ppmProd numeric; ppm tolerance for query product
#' @param q_ppmPrec numeric; ppm tolerance for query precursor
#' @param q_raThres numeric; Relative abundance threshold for query spectra
#' @param q_pol character; Polarity of query spectra ('positive', 'negative', NA).
#' @param q_instrumentTypes vector; Instrument types for query spectra.
#' @param q_instruments vector; Instruments for query spectra (note that this is used in combination with q_instrumentTypes - any
#'                              spectra matching either q_instrumentTypes or q_instruments will be used).
#' @param q_sources vector; Sources of query spectra (e.g. massbank, hmdb).
#' @param q_spectraTypes character; Spectra types of query spectra to perfrom spectral matching e.g. ('scans', 'av_all', 'intra', 'inter')
#' @param q_pids vector; pids for query spectra (correspond to column 'pid; in s_peak_meta)
#' @param q_rtrange vector; retention time range (in secs) of query spectra, first value mininum time and second value max e.g. c(0, 10) is between 0 and 10 seconds
#' @param q_spectraFilter boolean; For query spectra, if prior filtering performed with msPurity, flag peaks will be removed from spectral matching
#' @param q_xcmsGroups vector; XCMS group ids for query spectra
#' @param q_accessions vector; accession ids to filter query spectra
#'
#' @param l_purity character; Precursor ion purity threshold for the library spectra (uses interpolated purity - inPurity)
#' @param l_ppmProd numeric; ppm tolerance for library product
#' @param l_ppmPrec numeric; ppm tolerance for library precursor
#' @param l_raThres numeric; Relative abundance threshold for library spectra
#' @param l_pol character; Polarity of library spectra ('positive', 'negative', NA)
#' @param l_instrumentTypes vector; Instrument types for library spectra.
#' @param l_instruments vector; Instruments for library spectra (note that this is used in combination with q_instrumentTypes - any
#'                              spectra matching either q_instrumentTypes or q_instruments will be used).
#' @param l_sources vector; Sources of library spectra (e.g. massbank, hmdb).
#' @param l_spectraTypes vector; Spectra type of library spectra to perfrom spectral matching with e.g. ('scans', 'av_all', 'intra', 'inter')
#' @param l_pids vector; pids for library spectra (correspond to column 'pid; in s_peak_meta)
#' @param l_rtrange vector; retention time range (in secs) of library spectra, first value mininum time and second value max e.g. c(0, 10) is between 0 and 10 seconds
#' @param l_spectraFilter boolean; For library spectra, if prior filtering performed with msPurity, flag peaks will be removed from spectral matching
#' @param l_xcmsGroups vector; XCMS group ids for library spectra
#' @param l_accessions vector; accession ids to filter library spectra
#'
#' @param usePrecursors boolean; If TRUE spectra will be filtered by similarity of precursors based on ppm range defined by l_ppmPrec and q_ppmPrec
#' @param raW numeric; Relative abundance weight for spectra (default to 0.5 as determined by massbank for ESI data)
#' @param mzW numeric; mz weight for spectra (default to 2 as determined by massbank for ESI data)
#' @param rttol numeric ; Tolerance in time range between the library and query spectra retention time
#'
#' @param cores numeric; Number of cores to use
#' @param updateDb boolean; Update the Query SQLite database with the results
#' @param copyDb boolean; If updating the database - perform on a copy rather thatn the original query database
#' @param outPth character; If copying the database - the path of the new database file
#'
#' @return Returns a list containing the following elements
#'
#' **q_dbPth**
#'
#' Path of the query database (this will have been updated with the annotation results if updateDb argument used)
#'
#' **matchedResults**
#'
#' All matched results from the query spectra to the library spectra. Contains the following columns
#'
#' * dpc - dot product cosine of the match
#' * rdpc - reverse dot product cosine of the match
#' * cdpc - composite dot product cosine of the match
#' * mcount - number of matching peaks
#' * allcount - total number of peaks across both query and library spectra
#' * mpercent - percentage of matching peaks across both query and library spectra
#' * accession -  accession of library match
#' * name - name of library match
#' * inchikey - inchikey of library match
#' * lpid - pid in database of library match
#' * qpid - pid in database of query match
#' * mid - id of the match
#'
#' **xcmsMatchedResults**
#'
#' If the qeury spectra had XCMS based chromotographic peaks tables (e.g c_peak_groups, c_peaks) in the sqlite database - it will
#' be possible to summarise the matches for each XCMS grouped feature. The dataframe contains the following columns
#'
#' * pid - pid in database of query match
#' * grpid - grpid of the XCMS grouped feature for query match
#' * mz - derived from XCMS grouped feature
#' * mzmin - derived from XCMS grouped feature
#' * mzmax - derived from XCMS grouped feature
#' * rt - derived from XCMS grouped feature
#' * rtmin - derived from XCMS grouped feature
#' * rtmax - derived from XCMS grouped feature
#' * npeaks - derived from XCMS grouped feature
#' * grp_name - derived from XCMS grouped feature
#' * dpc - dot product cosine of the match
#' * rdpc - reverse dot product cosine of the match
#' * cdpc - composite dot product cosine of the match
#' * mcount - number of matching peaks
#' * allcount - total number of peaks across both query and library spectra
#' * mpercent - percentage of matching peaks across both query and library spectra
#' * accession -  accession of library match
#' * name - name of library match
#' * inchikey - inchikey of library match
#' * lpid - pid in database of library match
#' * mid - id of the match
#'
#' @return list of database details and dataframe summarising the results for the xcms features
#'
#' @examples
#' #msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' #xset <- xcms::xcmsSet(msmsPths)
#' #xset <- xcms::group(xset)
#' #xset <- xcms::retcor(xset)
#' #xset <- xcms::group(xset)
#'
#' #pa  <- purityA(msmsPths)
#' #pa <- frag4feature(pa, xset)
#' #pa <- filterFragSpectra(pa, allfrag=TRUE)
#' #pa <- averageAllFragSpectra(pa)
#' #q_dbPth <- createDatabase(pa, xset)
#' q_dbPth <- system.file("extdata", "tests", "db", "createDatabase_example.sqlite", package="msPurity")
#' result <- spectralMatching(q_dbPth, q_xcmsGroups = c(12, 27), cores=1, l_accessions=c('CCMSLIB00000577898','CE000616'))
#'
#' @importFrom magrittr %>%
#' @md
#' @export
spectralMatching <- function(
                             q_dbPth,
                             l_dbPth=NA,

                             q_purity=NA,
                             q_ppmProd=10,
                             q_ppmPrec=5,
                             q_raThres=NA,
                             q_pol='positive',
                             q_instrumentTypes=NA,
                             q_instruments=NA,
                             q_sources=NA,
                             q_spectraTypes='av_all',
                             q_pids=NA,
                             q_rtrange=c(NA, NA),
                             q_spectraFilter=TRUE,
                             q_xcmsGroups=NA,
                             q_accessions=NA,

                             l_purity=NA,
                             l_ppmProd=10,
                             l_ppmPrec=5,
                             l_raThres=NA,
                             l_pol='positive',
                             l_instrumentTypes=NA,
                             l_instruments=NA,
                             l_sources=NA,
                             l_spectraTypes=NA,
                             l_pids=NA,
                             l_rtrange=c(NA, NA),
                             l_spectraFilter=FALSE,
                             l_xcmsGroups=NA,
                             l_accessions=NA,
                             usePrecursors=TRUE,
                             raW=0.5,
                             mzW=2,
                             rttol=NA,

                             cores=1,
                             updateDb=FALSE,
                             copyDb=FALSE,
                             outPth='sm_result.sqlite'){
  message("Running msPurity spectral matching function for LC-MS(/MS) data")
  if(updateDb && copyDb){
    file.copy(from = q_dbPth, to=outPth)
    q_dbPth <- outPth
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
              pol = q_pol,
              instrumentTypes = q_instrumentTypes,
              instruments = q_instruments,
              sources = q_sources,
              pids = q_pids,
              rtrange = q_rtrange,
              con = q_con,
              xcmsGroups = q_xcmsGroups,
              spectraTypes = q_spectraTypes,
              accessions = q_accessions)

  ########################################################
  # Filter the library dataset
  ########################################################
  message("Filter library dataset")
  l_con <- DBI::dbConnect(RSQLite::SQLite(), l_dbPth)

  l_speakmeta <- filterSMeta(purity = l_purity,
                             raThres = l_raThres,
                             pol = l_pol,
                             instrumentTypes = l_instrumentTypes,
                             instruments = l_instruments,
                             sources = l_sources,
                             pids = l_pids,
                             rtrange = l_rtrange,
                             con = l_con,
                             xcmsGroups = l_xcmsGroups,
                             spectraTypes = l_spectraTypes,
                             accessions = l_accessions)



  ########################################################
  # Loop through the query dataset and spectra match
  # against the library spectra
  ########################################################
  # can't parallize dpylr without non cran package
  # Go back to using good old plyr


  q_fpids <- pullPid(q_speakmeta)

  l_fpids <- pullPid(l_speakmeta)


  message('aligning and matching')


  # Check cores and choose if parallel or not (do or dopar)
  if(cores>1){
    cl<-parallel::makeCluster(cores, type = "SOCK")
    doSNOW::registerDoSNOW(cl)
    parallel = TRUE
  } else{
    parallel = FALSE
  }

  # run parallel (or not) using foreach
  matched <- plyr::adply(q_fpids, 1, queryVlibrary, l_pids=l_fpids,
                                                    q_dbPth=q_dbPth,
                                                    l_dbPth=l_dbPth,
                                                    q_ppmPrec=q_ppmPrec,
                                                    q_ppmProd=q_ppmProd,
                                                    l_ppmPrec=l_ppmPrec,
                                                    l_ppmProd=l_ppmProd,
                                                    l_spectraFilter=l_spectraFilter,
                                                    q_spectraFilter=q_spectraFilter,
                                                    l_raThres=l_raThres,
                                                    q_raThres=q_raThres,
                                                    usePrecursors=usePrecursors,
                                                    mzW=mzW,
                                                    raW=raW,
                                                    rttol=rttol,
                          .parallel=parallel)

  # run parallel (or not) using foreach
  # matched <- operator(foreach::foreach(i = 1:2,
  #                                      .packages = c('dbplyr', 'magrittr', 'dplyr')),
  #                                      queryVlibrary(q_pid = 1,
  #                                                    q_speakmeta=q_speakmeta,
  #                                                    q_speaks=q_speaks,
  #                                                    l_speakmeta=l_speakmeta,
  #                                                    l_speaks=l_speaks,
  #                                                    q_ppmPrec=q_ppmPrec,
  #                                                    q_ppmProd=q_ppmProd,
  #                                                    l_ppmPrec=l_ppmPrec,
  #                                                    l_ppmProd=l_ppmProd,
  #                                                    raW=raW,
  #                                                    mzW=mzW,
  #                                                    rttol=rttol,
  #                                                    usePrecursors=usePrecursors)
  #                     )

  if(cores>1){
    parallel::stopCluster(cl)
  }




  if (nrow(matched)==0){
    message('No matches found')
    return(NULL)
  }


  # remove the plyr id column
  matched <- matched[,!names(matched)=='X1']
  matched$mid <- 1:nrow(matched)

  nmCols <- c("dpc","rdpc","cdpc","mcount", "allcount", "mpercent", "lpid", "qpid", "mid")
  matched[,nmCols] <- as.numeric(as.character(unlist(matched[,nmCols])))


  if (updateDb){

    custom_dbWriteTable(name_pk = 'mid', df=matched, fks=NA, table_name = 'sm_matches', con = q_con)

    if (DBI::dbExistsTable(l_con, "metab_compound") ){
      # Schema needs to be updated to be more generic

      if (DBI::dbExistsTable(l_con, "library_spectra_meta") ){
        l_compounds <- DBI::dbGetQuery(l_con, sprintf('SELECT  DISTINCT c.* FROM library_spectra_meta AS m
                                                          LEFT JOIN metab_compound AS
                                                          c on c.inchikey_id=m.inchikey_id
                                                          WHERE m.id IN (%s)', paste(unique(matched$lpid), collapse=",")) )
      }else{
        l_compounds <- DBI::dbGetQuery(l_con, sprintf('SELECT  DISTINCT c.* FROM s_peak_meta AS m
                                                          LEFT JOIN metab_compound AS
                                                          c on c.inchikey_id=m.inchikey_id
                                                          WHERE m.pid IN (%s)', paste(unique(matched$lpid), collapse=",")) )

      }

      if(nrow(l_compounds)>0){
        l_compounds <- data.frame(lapply(l_compounds, as.character), stringsAsFactors=FALSE)

        if(DBI::dbExistsTable(q_con, "metab_compound")){
          q_compounds <- q_con %>% dplyr::tbl("metab_compound")
          q_compounds <- q_compounds %>% dplyr::collect()
          q_inchi <- q_compounds$inchikey_id
          if(length(q_inchi)>0){
            l_compounds <- l_compounds[!l_compounds$inchikey_id %in%  q_inchi,]
          }
          if(length(l_compounds$inchikey_id[!is.na(l_compounds$inchikey_id)])>0){
            DBI::dbWriteTable(q_con, name="metab_compound", value=l_compounds, row.names=FALSE, append=TRUE)
          }
        }else{
          custom_dbWriteTable(name_pk = 'inchikey_id', fks = NA,
                              df=l_compounds, table_name = 'metab_compound', con = q_con, pk_type='TEXT')
        }
      }

      if (DBI::dbExistsTable(l_con, "library_spectra_meta") ){
        library_spectra_meta <- DBI::dbGetQuery(l_con, sprintf('SELECT  * FROM library_spectra_meta AS m
                                                          WHERE m.id IN (%s)', paste(unique(matched$lpid), collapse=",")) )
        pk = 'id'
      }else{
        library_spectra_meta <- DBI::dbGetQuery(l_con, sprintf('SELECT  * FROM s_peak_meta AS m
                                                          WHERE m.pid IN (%s)', paste(unique(matched$lpid), collapse=",")) )
        pk = 'pid'
      }

      #fk_l = list('inchikey_id'=list('new_name'='inchikey_id', 'ref_name'='inchikey_id', 'ref_table'='metab_compound'))

      custom_dbWriteTable(name_pk = pk, fks = NA,
                          df=library_spectra_meta, table_name = 'l_s_peak_meta', con = q_con)

    }




  }

  ########################################################
  # Create a summary table for xcms grouped objects
  ########################################################
  if (DBI::dbExistsTable(q_con, "c_peak_groups")){
    # check if the query is from an XCMS object
    message("Summarising LC features annotations")

    xcmsMatchedResults <- getXcmsSmSummary(q_con, matched,spectraTypes=q_spectraTypes)
    if (updateDb){
      DBI::dbWriteTable(q_con, name='xcms_match', value=xcmsMatchedResults, row.names=F, append=T)
    }

  }else{
    xcmsMatchedResults <- NA
  }

  DBI::dbDisconnect(q_con)
  DBI::dbDisconnect(l_con)


  return(list('q_dbPth' = q_dbPth, 'matchedResults' = matched,
              'xcmsMatchedResults' = xcmsMatchedResults))
}

# filterPid <- function(sp, pids){
#   nms <- names(sp %>% dplyr::collect())
#   if ("library_spectra_meta_id" %in% nms){
#     sp <- sp %>% dplyr::filter(library_spectra_meta_id %in% pids)
#   }else if ("pid" %in% nms){
#     sp <- sp %>% dplyr::filter(pid %in% pids)
#   }else{
#     sp <- sp %>% dplyr::filter(id %in% pids)
#   }
#   return(sp)
# }

pullPid <- function(sp, pids){
  tble <- sp %>% dplyr::collect()
  nms <- colnames(tble)

  if ("pid" %in% nms){
    pids <- tble$pid
  }else{
    pids <- tble$id
  }

  return(pids)
}

getScanPeaksSqlite <- function(con, spectraFilter=TRUE, spectraTypes=NA, raThres=NA, pids=NA){


  if (DBI::dbExistsTable(con, "s_peaks")){
    speaks <- con %>% dplyr::tbl("s_peaks")
  }else if (DBI::dbExistsTable(con, "library_spectra")) {
    # old sqlite format
    speaks <- con %>% dplyr::tbl("library_spectra")
  }else{
    stop('No spectra available')
  }

  cn <- getPeakCols(con)

  if ('pid' %in% cn$name  && !anyNA(pids)){
    speaks <- speaks %>% dplyr::filter(pid %in% pids)
  }else if("library_spectra_meta_id" %in% cn$name  && !anyNA(pids)){
    speaks <- speaks %>% dplyr::filter(library_spectra_meta_id %in% pids)
  }

  if ('pass_flag' %in% cn$name && spectraFilter ){
     speaks <- speaks %>% dplyr::filter(pass_flag==TRUE)
  }

  if ('type' %in% cn$name && !anyNA(spectraTypes)){
    speaks <- speaks %>% dplyr::filter(type %in% spectraType)
  }

  if ('ra' %in% cn$name && !is.na(raThres)){
    speaks <- speaks %>% dplyr::filter(ra>raThres)
  }

  return(speaks)



}

getXcmsSmSummary <- function(con, matched, scoreF=0, fragNmF=1, spectraTypes='scans'){

  if ('scans' %in% spectraTypes){
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


  return(xcmsMatchedResults)


}

filterSMeta <- function(purity=NA,
                        raThres=0,
                        pol='positive',
                        instrumentTypes=NA,
                        instruments=NA,
                        sources=NA,
                        xcmsGroups=NA,
                        pids=NA,
                        rtrange=c(NA, NA),
                        spectraTypes=NA,
                        accessions=NA,
                        con){

  # get column names
  # PRAGMA table_info();
  meta_cn <- getMetaCols(con)

  speakmeta <- getSmeta(con, pids)

  if('accession' %in% meta_cn$name && !anyNA(accessions)){
    speakmeta <- speakmeta %>% dplyr::filter(accession %in% accessions)
  }
  #print('accession')
  #print(speakmeta)

  if('inPurity' %in% meta_cn$name && !is.na(purity)){
    speakmeta <- speakmeta %>% dplyr::filter(inPurity > purity)
  }
  #print('purity')
  #print(speakmeta)

  if ('polarity' %in% meta_cn$name && !is.na(pol)){
    speakmeta <- speakmeta %>% dplyr::filter(lower(polarity) == lower(pol))
  }

  #print('polarity')
  #print(speakmeta)


  if ('instrument_type' %in% meta_cn$name &&  'instrument' %in% meta_cn$name &&  !anyNA(instrumentTypes)  && !anyNA(instruments)){
    speakmeta <- speakmeta %>% dplyr::filter(instrument_type %in% instrumentTypes || instrument %in% instruments)
  }else if ('instrument_type' %in% meta_cn$name && !anyNA(instrumentTypes)){
    speakmeta <- speakmeta %>% dplyr::filter(instrument_type %in% instrumentTypes)
  }else if ('instrument' %in% meta_cn$name && !anyNA(instruments)){
    speakmeta <- speakmeta %>% dplyr::filter(instrument %in% instruments)
  }

  #print('instruments')
  #print(speakmeta)


  if(!anyNA(sources)){
    if (DBI::dbExistsTable(con, "library_spectra_source")){
      sourcetbl <- con %>% dplyr::tbl("library_spectra_source")
      speakmeta <- speakmeta %>% dplyr::left_join(sourcetbl,  by=c('library_spectra_source_id'='id'),suffix = c("", ".y")) %>%
        dplyr::filter(name.y %in% sources)

    }else if (DBI::dbExistsTable(con, "source")){
      sourcetbl <- con %>% dplyr::tbl("source")
      speakmeta <- speakmeta %>% dplyr::left_join(sourcetbl,  by=c('sourceid'='id'),suffix = c("", ".y")) %>%
        dplyr::filter(name.y %in% sources)
    }
  }

  #print('sources')
  #print(speakmeta)



  if('retention_time' %in% meta_cn$name && !anyNA(rtrange)){
    speakmeta <- speakmeta %>% dplyr::filter(retention_time >= rtrange[1] & retention_time <= rtrange[2])
  }

  #print('rtrange')
  #print(speakmeta)


  if ( !anyNA(xcmsGroups) && DBI::dbExistsTable(con, "c_peak_groups")){

    XLI <- DBI::dbGetQuery(con, paste0('SELECT cpg.grpid, spm.pid FROM s_peak_meta AS spm
                           LEFT JOIN c_peak_X_s_peak_meta AS cXs ON cXs.pid=spm.pid
                           LEFT JOIN c_peaks AS cp on cp.cid=cXs.cid
                           LEFT JOIN c_peak_X_c_peak_group AS cXg ON cXg.cid=cp.cid
                           LEFT JOIN c_peak_groups AS cpg ON cXg.grpid=cpg.grpid
                           WHERE cpg.grpid in (', paste(xcmsGroups, collapse = ','), ')'))

    xcmsGroups <- as.character(xcmsGroups)

    # doesn't work with database calls on travis (have to split into to filters)
    # speakmeta <- speakmeta %>% dplyr::filter(grpid %in% xcmsGroups | pid %in% XLI$pid)
    metaGrpPids <- speakmeta %>% dplyr::filter(grpid %in% xcmsGroups) %>% dplyr::pull(pid)
    allGrpPids <- unique(c(XLI$pid,metaGrpPids))
    speakmeta <- speakmeta %>% dplyr::filter(pid %in% allGrpPids)



  }

  #print('xcms groups')
  #print(speakmeta)


  if ('spectrum_type' %in% meta_cn$name && !is.na(spectraTypes)){
    if ('av_all' %in% spectraTypes){
      spectraTypes[spectraTypes=='av_all'] = 'all'
    }



    speakmeta <- speakmeta %>% dplyr::filter(spectrum_type %in% spectraTypes)
  }


  return(speakmeta)


}

filterPrecursors <- function(q_pid, q_speakmeta, l_speakmeta, q_ppmPrec, l_ppmPrec){

  return(l_speakmetaFiltered)

}

getSmeta <- function(con, pids=NA){
  if (DBI::dbExistsTable(con, "s_peak_meta")){
    speakmeta <- con %>% dplyr::tbl("s_peak_meta")
    if (!anyNA(pids)){
      speakmeta <- speakmeta  %>%  dplyr::filter(pid %in% pids)
    }

  }else if (DBI::dbExistsTable(con, "library_spectra_meta")) {
    # old sqlite format
    speakmeta <- con %>% dplyr::tbl("library_spectra_meta")
    if (!anyNA(pids)){

      speakmeta <- speakmeta  %>% dplyr::filter(id %in% pids)
    }

  }else{
    stop('No meta data for spectra available')
  }
  return(speakmeta)
}

getMetaCols <- function(con){
  if(DBI::dbExistsTable(con, "s_peak_meta")){
    meta_cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(s_peak_meta)')
  }else{
    meta_cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(library_spectra_meta)')
  }
  return(meta_cn)
}

getPeakCols <- function(con){
  if(DBI::dbExistsTable(con, "s_peak_meta")){
    cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(s_peaks)')
  }else{
    cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(library_spectra)')
  }
  return(cn)
}

queryVlibrary <- function(q_pid, l_pids, q_dbPth, l_dbPth, q_ppmPrec, q_ppmProd, l_ppmPrec, l_ppmProd,
                          l_spectraFilter, q_spectraFilter, l_raThres, q_raThres, usePrecursors, mzW, raW, rttol){

  q_con <- DBI::dbConnect(RSQLite::SQLite(), q_dbPth)
  l_con <- DBI::dbConnect(RSQLite::SQLite(), l_dbPth)


  q_speakmetai <- getSmeta(q_con, q_pid) %>% dplyr::collect()

  q_speaksi <- getScanPeaksSqlite(q_con, spectraFilter=q_spectraFilter, pids=q_pid) %>% dplyr::collect()

  # if no peaks, then skip
  if (nrow(q_speaksi)==0){
    return(NULL)
  }



  l_speakmeta <- getSmeta(l_con, l_pids) %>% dplyr::collect()
  l_speaks <- getScanPeaksSqlite(l_con, spectraFilter=l_spectraFilter, pids=l_pids) %>% dplyr::collect()

  if(usePrecursors){

    q_precMZ <- q_speakmetai$precursor_mz
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
  }else{
    l_fspeakmeta <- l_speakmeta %>% dplyr::collect()
  }


  if(!is.na(rttol)){

    l_fspeakmeta <- l_fspeakmeta %>% dplyr::filter(abs(retention_time-q_speakmetai)<rttol)
  }

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
                            l_ppmProd=l_ppmProd,
                            raW=raW,
                            mzW=mzW

              )


  searched$qpid <- q_pid


  DBI::dbDisconnect(q_con)
  DBI::dbDisconnect(l_con)



  return(searched)



}


queryVlibrarySingle <- function(l_pid, q_speaksi, l_speakmeta, l_speaks, q_ppmProd, l_ppmProd, raW, mzW){




  if ('pid' %in% colnames(l_speaks)){
    l_speaksi <- l_speaks %>% dplyr::filter(pid==l_pid) %>% dplyr::collect()
    l_speakmetai <- data.frame(l_speakmeta %>% dplyr::filter(pid==l_pid) %>% dplyr::collect())
  }else{
    l_speaksi <- l_speaks %>% dplyr::filter(library_spectra_meta_id==l_pid) %>% dplyr::collect()
    l_speakmetai <- data.frame(l_speakmeta %>% dplyr::filter(id==l_pid) %>% dplyr::collect())
  }


  # ensure we have the relative abundance
  l_speaksi$ra <- (l_speaksi$i/max(l_speaksi$i))*100
  q_speaksi$ra <- (q_speaksi$i/max(q_speaksi$i))*100


  am <- alignAndMatch(q_speaksi, l_speaksi, q_ppmProd, l_ppmProd, raW, mzW)

  if ('pid' %in% colnames(l_speakmetai)){
    lpids <- l_speakmetai$pid
  }else{
    lpids <- l_speakmetai$id
  }

  return(c(am, 'accession'=l_speakmetai$accession,
               'name'=l_speakmetai$name,
               'inchikey'=l_speakmetai$inchikey_id,
               'lpid'=lpids
           ))
}


alignAndMatch <- function(q_speaksi, l_speaksi, q_ppmProd, l_ppmProd, raW, mzW){
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









