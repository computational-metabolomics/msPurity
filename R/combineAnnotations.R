#' @title Combine Annotations
#'
#' @description
#'  Combine the annotation results from msPurity spectral matching, MetFrag, Sirius CSI:FingerID and probmetab
#'  based on weighted scores for each technique aligning each annotation by inchikey and XCMS grouped feature.
#'
#'  Advised that the tool is run with a local compound database that includes (example found here xx)
#'
#'  Can be run without a local compound database but will take a very long time to finish.
#'
#' @param sm_resultPth character;
#' @param metfrag_resultPth character;
#' @param sirius_csi_resultPth character;
#' @param probmetab_resultPth character;
#' @param weights list;
#' @param localCompoundDbPth character;
#' @param outPth character;
#'
#' @examples
#' metfrag_resultPth <- system.file("extdata", "tests", "external_annotations", "metfrag.tsv", package="msPurity")
#' # run the standard spectral matching workflow to get the sm_resultPth
#' sm_resultPth <- system.file("extdata","tests", "sm", "spectralMatching_result.sqlite", package="msPurity")
#' localCompoundDbPth <- system.file("extdata", "tests", "db", "compounds_18July2019_0319.sqlite", package="msPurity")
#' combined <- combineAnnotations(sm_resultPth,
#'                                metfrag_resultPth,
#'                                outPth=file.path(tempdir(), 'combined.sqlite'),
#'                                localCompoundDbPth=localCompoundDbPth)
#' @return purityA object with slots for fragmentation-XCMS links
#' @export
combineAnnotations <- function(sm_resultPth,
                               localCompoundDbPth,
                               metfrag_resultPth=NA,
                               sirius_csi_resultPth=NA,
                               probmetab_resultPth=NA,
                               weights=list('sm'=0.2,
                                            'metfrag'=0.15,
                                            'sirius_csifingerid'=0.15,
                                            'cfm'=0.15,
                                            'probmetab'=0.05,
                                            'biosim'=0.3
                               ),
                               outPth=NA

                               ){

  if(!is.na(outPth)){
    file.copy(sm_resultPth, outPth)
    sqlitePth <- outPth
  }else{
    sqlitePth <- sm_resultPth
  }

  # sm, metfrag, sirius, probmetab, lipidsearch, mzcloud, bs
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlitePth)

  if (!is.null(localCompoundDbPth)){
    con_comp <- DBI::dbConnect(RSQLite::SQLite(), localCompoundDbPth)
  }else{
    con_comp <- NULL
  }

  #############################################
  # Add all the annotation results
  #############################################
  message('Add Metfrag results')
  addMetFragResults(metfrag_resultPth, con)
  message('Add Sirius CSI:FingerID results')
  addSiriusResults(sirius_csi_resultPth, con, silentRestErrors=silentRestErrors, con_comp=con_comp)
  message('Add Probmetab results')
  addProbmetabResults(probmetab_resultPth, con, silentRestErrors=silentRestErrors, con_comp=con_comp)
  # Add lipidsearch result TODO
  # Add mzcloud result TODO

  ####################################
  # Update the metab_compound
  ####################################
  # We now have inchikeys for all the annotations. For each inchikey we get all the associated informations
  message('Update the compound table')
  metab_compounds <- RSQLite::dbGetQuery(con, 'SELECT * FROM metab_compound')
  metab_compounds <- updateColName(metab_compounds, 'inchikey_id', 'inchikey')

  # Get all details for inchikeys (using local compoundDB - currently we only have the inchikeys)
  metab_compounds_m <- inchikeySelect(metab_compounds$inchikey, con = con_comp, all=TRUE)

  # Add anything that is not found (e.g. will be a few cases where pubchem won't have inchikey but the spectral matching result will have)
  missing_comps <- metab_compounds[!metab_compounds$inchikey %in% metab_compounds_m$inchikey,]

  missing_comps <- cbind(missing_comps, splitInchis(missing_comps$inchikey))
  colnames(missing_comps)[colnames(missing_comps)=='pubchem_id'] = 'pubchem_cids'
  # Remove rows we won't be using from missing_comps
  keep_cols = c('inchikey', 'inchikey1', 'inchikey2', 'inchikey3', 'pubchem_cids',
                  'exact_mass', 'molecular_formula', 'name')
  missing_comps <- missing_comps[,keep_cols]
  findx <- sapply(missing_comps, is.factor)
  missing_comps[findx] <- lapply(missing_comps[findx], as.character)
  metab_compounds_m <- dplyr::bind_rows(metab_compounds_m, missing_comps)


  # remove the old table and add the modified table
  DBI::dbExecute(con, 'DROP TABLE metab_compound')
  DBI::dbWriteTable(conn=con, name='metab_compound', value=metab_compounds_m, row.names=FALSE)

  # Calculate final score and rank
  # Build up the table by looping through xcms group
  message('Combine annotations')
  c_peak_groups <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups')
  combined_scores <- plyr::ddply(c_peak_groups, ~grpid, combineScoresGrp, weights=weights, con=con)

  DBI::dbWriteTable(conn=con, name='combined_annotations', value=combined_scores, row.names=FALSE)

  return(getAnnotationSummary(con))

}



addMetFragResults <- function(metfrag_resultPth, con, con_comp){
  # Add metfrag details
  if (!is.na(metfrag_resultPth) && file.exists(metfrag_resultPth)){
    DBI::dbWriteTable(conn=con, name='metfrag_results', value=metfrag_resultPth, sep='\t', header=T)

    DBI::dbExecute(con,
                     'INSERT INTO metab_compound (inchikey_id)
                   SELECT DISTINCT m.InChIKey
                   FROM metfrag_results AS m
                   LEFT JOIN
                   metab_compound AS c ON m.InChIKey = c.inchikey_id
                   WHERE c.inchikey_id IS NULL')

  }
}



addSiriusResults <- function(sirius_csi_resultPth, con, silentRestErrors, con_comp=NULL){
  # Add sirius details
  if (!is.na(sirius_csi_resultPth) && file.exists(sirius_csi_resultPth)){
    DBI::dbWriteTable(conn=con, name='sirius_csifingerid_results', value=sirius_csi_resultPth, sep='\t', header=T, row.names=T, nrows = 4)

    inchikey2ds <- DBI::dbGetQuery(con, "SELECT DISTINCT inchikey2D FROM sirius_csifingerid_results")

    # get all inchikeys from local database
    compDetails <- inchikeySelect(inchikey2ds$inchikey2D, inchikeyCol = 'inchikey1', con=con_comp)

    # check current inchikeys and add
    new_compDetails <- insertNewInchikeys(con, con_comp, compDetails)

    ##################
    # calculate Sirius score
    ##################
    # Loop through sqlite database by grpid.
    uids <- DBI::dbGetQuery(con, "SELECT DISTINCT grpid FROM sirius_csifingerid_results")

    rs <- DBI::dbSendQuery(con, "SELECT rowid, Score FROM sirius_csifingerid_results WHERE grpid = ?")

    bounded_score <- plyr::adply(uids$grpid, 1, getBoundedSiriusScore, rs=rs)
    DBI::dbClearResult(rs)
    # caculate minmax normalised value (save)
    # Add new column
    DBI::dbExecute(con, "ALTER TABLE sirius_csifingerid_results ADD COLUMN bounded_score text")
    bounded_score$bounded_score <- as.numeric(bounded_score$bounded_score)
    bounded_score <- bounded_score[,c('rowid', 'bounded_score')]

    DBI::dbExecute(con, "UPDATE sirius_csifingerid_results SET bounded_score = :bounded_score  WHERE rowid = :rowid",
                   params=bounded_score)

  }
}

getBoundedSiriusScore <- function(x, rs){
  DBI::dbBind(rs, unique(x))
  scores <- DBI::dbFetch(rs)
  if(nrow(scores)==1){
    bounded_score <- 1
  }else{
    bounded_score <- negMinMax(abs(scores$Score))
  }
  scores$bounded_score <- bounded_score
  return(scores)

}

negMinMax <- function(x){
  xn <- (x-min(x))/(max(x)-min(x))
  return(abs(xn-1))

}



addProbmetabResults <- function(probmetab_resultPth, con, silentRestErrors, con_comp=NULL){
  if (!is.na(probmetab_resultPth) && file.exists(probmetab_resultPth)){
    addProbmetab(probmetab_resultPth, con)

    # Fetch in chunks
    kegg_cids <- DBI::dbGetQuery(con, "SELECT DISTINCT mpc FROM probmetab_results")

    kegg_cid_str <- vecToStringStmt(kegg_cids$mpc)
    compDetails = DBI::dbGetQuery(con_comp, sprintf("SELECT kegg_cid, inchikey FROM kegg WHERE kegg_cid IN (%s)", kegg_cid_str))
    DBI::dbWriteTable(con, name='kegg', value=compDetails, row.names=FALSE)

    # Check if any of the same inchikeys have been added
    new_compDetails <- insertNewInchikeys(con, con_comp, compDetails)

  }


}


addProbmetab <- function(pth, con){
  if (!is.null(pth)){

    df <- utils::read.table(pth,  header = TRUE, sep='\t', stringsAsFactors = FALSE,  comment.char = "")

    nmap <- DBI::dbGetQuery(con, 'SELECT grpid, grp_name FROM c_peak_groups')


    df$grpid <- nmap$grpid[match(df$name, nmap$grp_name)]


    start <- T
    for (i in 1:nrow(df)){

      x <- df[i,]

      if(is.na(x$proba) | x$proba =='NA'){
        next
      }
      mpc <- gsub("cpd:", "", x$mpc)
      mpc <- stringr::str_split(mpc, ';')


      proba <- stringr::str_split(x$proba, ';')

      for (j in 1:length(mpc[[1]])){

        row <-  c(x$grpid, x$propmz, mpc[[1]][j], proba[[1]][j])

        if (start){
          df_out <- data.frame(t(row), stringsAsFactors=F)
          start <- F
        }else{
          df_out <- data.frame(rbind(df_out, row), stringsAsFactors=F)
        }
      }

    }

    colnames(df_out) <- c('grpid', 'propmz', 'mpc', 'proba')
    DBI::dbWriteTable(con, name='probmetab_results', value=df_out, row.names=FALSE)

  }
}


getAnnotationSummary <- function(con){
  meta_cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(l_s_peak_meta)')
  if ("pid" %in% meta_cn){
    pid_id <- 'pid'
  }else{
    pid_id <- 'id'
  }
  sql_stmt <- sprintf('
                      SELECT
                      cpg.grpid,
                      cpg.grp_name,
                      cpg.mz,
                      cpg.rt,
                      mc.inchikey,
                      mc.inchikey1,
                      mc.inchikey2,
                      mc.inchikey3,
                      mc.name,
                      mc.exact_mass,
                      mc.molecular_formula,
                      mc.pubchem_cids,
                      mc.kegg_cids,
                      mc.brite,
                      mc.kegg_drug,
                      GROUP_CONCAT(DISTINCT(cast(spm.acquisitionNum  as INTEGER) )) AS fragmentation_acquistion_num,
                      ROUND(AVG(spm.inPurity),3) AS mean_precursor_ion_purity,
                      l.accession,
                      ca.sirius_score,
                      ca.sirius_wscore,
                      ca.metfrag_score,
                      ca.metfrag_wscore,
                      ca.sm_score,
                      ca.sm_wscore,
                      ca.probmetab_score,
                      ca.probmetab_wscore,
                      ca.wscore,
                      ca.rank
                      FROM combined_annotations AS ca
                      LEFT JOIN
                      metab_compound AS mc ON ca.inchikey = mc.inchikey
                      LEFT JOIN
                      l_s_peak_meta AS l ON l.%s = ca.sm_lpid
                      LEFT JOIN
                      c_peak_groups AS cpg ON cpg.grpid = ca.grpid
                      LEFT JOIN
                      c_peak_X_c_peak_group AS cpXcpg ON cpXcpg.grpid = cpg.grpid
                      LEFT JOIN
                      c_peaks AS cp ON cp.cid = cpXcpg.cid
                      LEFT JOIN
                      c_peak_X_s_peak_meta AS cpXspm ON cpXspm.cid = cp.cid
                      LEFT JOIN
                      s_peak_meta AS spm ON spm.pid = cpXspm.pid
                      GROUP BY ca.grpid, ca.inchikey
                      ORDER BY ca.grpid, ca.rank
                      ', pid_id)



  return(DBI::dbGetQuery(conn=con, sql_stmt))


}


combineScoresGrp <- function(c_peak_group, weights, con){

  # get sirius data
  grpid <- unique(c_peak_group$grpid)

  if(DBI::dbExistsTable(con, 'sirius_csifingerid_results')){
    sirius <- DBI::dbGetQuery(con, 'SELECT
                              s.rowid AS sirius_id,
                              s.bounded_score AS sirius_score,
                              mc.inchikey
                              FROM metab_compound AS mc
                              LEFT JOIN
                              sirius_csifingerid_results AS s
                              ON s.inchikey2D = mc.inchikey1 WHERE s.grpid= :grpid', params=list('grpid'=grpid))
    # can have multiple hits for each inchikey (e.g. multiple scans)
    sirius <- plyr::ddply(sirius, ~inchikey, getBestScore, 'sirius_score')


    sirius$sirius_wscore <- as.numeric(sirius$sirius_score)*weights$sirius_csifingerid

  }else{
    sirius <- data.frame(inchikey=NA, sirius_score=NA, sirius_wscore=NA)
  }


  if(DBI::dbExistsTable(con, 'metfrag_results')){
    # get metfrag
    metfrag <- DBI::dbGetQuery(con, 'SELECT
                               m.rowid as metfrag_id,
                               m.Score AS metfrag_score,
                               mc.inchikey
                               FROM metab_compound AS mc
                               LEFT JOIN
                               metfrag_results AS m
                               ON m.InChIKey = mc.inchikey WHERE m.grpid=:grpid', params=list('grpid'=grpid))
    # can have multiple hits for each inchikey (e.g. multiple scans)
    metfrag <- plyr::ddply(metfrag, ~inchikey, getBestScore, 'metfrag_score')

    metfrag$metfrag_wscore <- as.numeric(metfrag$metfrag_score)*weights$metfrag

  }else{
    metfrag <- data.frame(inchikey=NA, metfrag_score=NA, metfrag_wscore=NA)

  }

  meta_cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(l_s_peak_meta)')
  if ("pid" %in% meta_cn){
    pid_id <- 'pid'
  }else{
    pid_id <- 'id'
  }

  # get spectral matching
  sm <- DBI::dbGetQuery(con, sprintf('SELECT sm.lpid AS sm_lpid,
                                     sm.mid AS sm_mid,
                                     sm.dpc AS sm_score,
                                     mc.inchikey
                                     FROM metab_compound AS mc
                                     LEFT JOIN
                                     l_s_peak_meta AS l ON l.inchikey_id = mc.inchikey
                                     LEFT JOIN
                                     xcms_match AS sm ON l.%s = sm.lpid
                                     WHERE sm.grpid = :grpid', pid_id), params=list('grpid'=grpid))

  # can have multiple hits for each inchikey (e.g. multiple scans and library spectra for each inchikey)
  # we just take the best hit
  sm <- plyr::ddply(sm, ~inchikey, getBestScore, 'sm_score')
  sm$sm_wscore <- as.numeric(sm$sm_score)*weights$sm


  # get probmetab
  if(DBI::dbExistsTable(con, 'probmetab_results')){
    probmetab <- DBI::dbGetQuery(con, paste('SELECT p.rowid as probmetab_id,
                                            p.proba AS probmetab_score,
                                            mc.inchikey
                                            FROM metab_compound AS mc
                                            LEFT JOIN
                                            kegg as k ON k.inchikey=mc.inchikey
                                            LEFT JOIN
                                            probmetab_results AS p ON p.mpc = k.kegg_cid
                                            WHERE p.grpid= ', grpid, sep=' '))
    # can have multiple hits for each inchikey (e.g. multiple scans)
    probmetab <- plyr::ddply(probmetab, ~inchikey, getBestScore, 'probmetab_score')
    probmetab$probmetab_wscore <- as.numeric(probmetab$probmetab_score)*weights$probmetab
  }else{
    probmetab <- data.frame(inchikey=NA, probmetab_score=NA, probmetab_wscore=NA)
  }

  score_list <- list(sirius, metfrag, sm, probmetab)

  # combine all
  combined_df <- Reduce(function(...) merge(..., all=TRUE, by='inchikey'), score_list)

  # Get biosim score (if available)
  # todo

  combined_df <- combined_df[!is.na(combined_df$inchikey),]
  combined_df[is.na(combined_df)] <- 0
  combined_df$wscore <- rowSums(combined_df[, c('sirius_wscore', 'metfrag_wscore', 'sm_wscore', 'probmetab_wscore')])

  # Order
  combined_df <- combined_df[order(combined_df$wscore, decreasing=TRUE),]

  if (nrow(combined_df)==0){
    return(NULL)
  }else{

    # Add rank
    combined_df$rank <- as.numeric(factor(rank(-combined_df$wscore)))
  }

  return(combined_df)
}


splitInchis <- function(inchikeys){
  inchisplit <- data.frame(stringr::str_split_fixed(inchikeys, '-', 3))
  colnames(inchisplit) <- c('inchikey1', 'inchikey2', 'inchikey3')
  return(inchisplit)
}


updateCol<- function(df, old, new){

  df[,old][!is.na(df[,new])] <- df[,new][!is.na(df[,new])]

  return(df)

}

getBestScore <- function(x, nm){
  # get the best hit
  x <- x[x[,nm]==max(x[,nm]),]

  # if multiple with the best hit just use the top hit
  return(x[1,])
}

updateColName <- function(df, old, new){
  colnames(df)[colnames(df)==old] <- new
  return(df)
}

vecToStringInsertStmt <- function(vec){
  return(paste0("('", unique(vec), "')", sep="", collapse=","))
}

inchikeyInsert <- function(inchikeys, con){
  inchi_str <- vecToStringInsertStmt(inchikeys)
  insert_stmt <- sprintf('INSERT INTO metab_compound (inchikey_id) VALUES %s', inchi_str)
  DBI::dbExecute(con, insert_stmt)
}

vecToStringInStmt <- function(vec){
  return(paste0("'", unique(vec), "'", sep="", collapse=","))
}

inchikeySelect <- function(inchikeys, con, inchikeyCol='inchikey', outCols='inchikey', all=FALSE){
  inchi_str <- vecToStringInStmt(inchikeys)
  if (all){
    sql_stmt <- sprintf('SELECT * FROM metab_compound WHERE %s IN (%s)', inchikeyCol, inchi_str)
  }else{
    sql_stmt <- sprintf('SELECT %s FROM metab_compound WHERE %s IN (%s)', outCols, inchikeyCol, inchi_str)
  }

  return(RSQLite::dbGetQuery(con, sql_stmt))
}

insertNewInchikeys <- function(con, con_comp, compDetails){
  # Check for new compounds
  matching_inchi <- inchikeySelect(compDetails$inchikey, inchikeyCol = 'inchikey_id', con=con, outCols="inchikey_id")

  new_compDetails <- compDetails[!compDetails$inchikey %in% matching_inchi$inchikey_id,,drop=FALSE]

  new_compDetails <- new_compDetails[!duplicated(new_compDetails$inchikey),,drop=FALSE]

  if (nrow(new_compDetails)>0){
    # Add new full inchikeys to metab_compound (full details can be added later)
    inchikeyInsert(new_compDetails$inchikey, con)
  }

  return(new_compDetails)
}
