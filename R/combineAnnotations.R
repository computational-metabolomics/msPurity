#' @title Combine Annotations
#'
#' @description
#'  Combine the annotation results from msPurity spectral matching, MetFrag, Sirius CSI:FingerID, probmetab and any generic
#'  MS1 lookup software (e.g. results from the BEAMS software)
#'
#'  The annotation results are then aligned by inchikey and XCMS grouped feature.
#'
#'  The tool has to be run with a local compound database (available on request - contact t.n.lawson@@bham.ac.uk)
#'
#'
#' @param sm_resultPth character; Path to the msPurity SQLite database used for spectral matching
#' @param metfrag_resultPth character; Path to the tsv table of metfrag results
#' @param sirius_csi_resultPth character; Path to the tsv table of Sirius CSI:Finger ID results
#' @param probmetab_resultPth character; Path to the tsv table of Probmetab results
#' @param ms1_lookup_resultPth character; Path to generic tsv table of MS1 lookup results
#' @param ms1_lookup_dbSource character; Source of the compound database used for ms1_lookup
#'                            (currently only supports HMDB, KEGG or PubChem)
#' @param ms1_lookup_checkAdducts boolean; Check if adducts match to those found in CAMERA (requires the database to have been created with CAMERA object)
#' @param ms1_lookup_keepAdducts vecotr; Keep only adducts found from the MS1 lookup that are found in this vector
#' @param weights list;
#' @param compoundDbPth character; Path to local compound database with pubchem, hmdb, KEGG and metab_compound summary table
#'                                 (full database available on request - contact t.n.lawson@@bham.ac.uk).
#'                                 This is only applicable if using "compoundDbType sqlite")
#' @param compoundDbType character; Database type for compound database can be either (sqlite, postgres or mysql)
#' @param compoundDbName character; Database name (only applicable for postgres and mysql)
#' @param compoundDbHost character; Database host (only applicable for postgres and mysql)
#' @param compoundDbPort character; Database port (only applicable for postgres and mysql)
#' @param compoundDbUser character; Database user (only applicable for postgres and mysql)
#' @param compoundDbPass character; Database pass (only applicable for postgres and mysql) - Note this is not secure!
#' @param summaryOutput boolean; If a summary dataframe is to be created
#'
#'
#' @param outPth character;
#'
#' @examples
#' metfrag_resultPth <- system.file("extdata", "tests", "external_annotations",
#'                                       "metfrag.tsv", package="msPurity")
#' # run the standard spectral matching workflow to get the sm_resultPth
#' sm_resultPth <- system.file("extdata","tests", "sm",
#'                         "spectralMatching_result.sqlite", package="msPurity")
#' compoundDbPth <- system.file("extdata", "tests", "db",
#'                           "metab_compound_subset.sqlite", package="msPurity")
#' combined <- combineAnnotations(sm_resultPth,
#'                              metfrag_resultPth,
#'                              outPth=file.path(tempdir(), 'combined.sqlite'),
#'                              compoundDbPth=compoundDbPth)
#' @return purityA object with slots for fragmentation-XCMS links
#' @export
combineAnnotations <- function(sm_resultPth,
                               compoundDbPth,
                               metfrag_resultPth=NA,
                               sirius_csi_resultPth=NA,
                               probmetab_resultPth=NA,
                               ms1_lookup_resultPth=NA,
                               ms1_lookup_dbSource='hmdb',
                               ms1_lookup_checkAdducts=FALSE,
                               ms1_lookup_keepAdducts=c('[M+H]+', '[M-H]-'),
                               weights=list('sm'=0.3,
                                            'metfrag'=0.2,
                                            'sirius_csifingerid'=0.2,
                                            'probmetab'=0,
                                            'ms1_lookup'=0.05,
                                            'biosim'=0.25
                               ),
                               compoundDbType='sqlite',
                               compoundDbName=NA,
                               compoundDbHost=NA,
                               compoundDbPort=NA,
                               compoundDbUser=NA,
                               compoundDbPass=NA,
                               outPth=NA,
                               summaryOutput=True

){

  if(!is.na(outPth)){
    file.copy(sm_resultPth, outPth)
    sqlitePth <- outPth
  }else{
    sqlitePth <- sm_resultPth
  }

  # sm, metfrag, sirius, probmetab, lipidsearch, mzcloud, bs
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlitePth)

  if (compoundDbType=='sqlite'){
    con_comp <- DBI::dbConnect(RSQLite::SQLite(), compoundDbPth)
  }else if (compoundDbType=='postgres'){
    con_comp <- DBI::dbConnect(RPostgres::Postgres(),
                               user=compoundDbUser,
                               password=compoundDbPass,
                               dbname=compoundDbname,
                               host=compoundDbHost,
                               pass=compoundDbPass)
  }else if (compoundDbType=='mysql'){
    con_comp <- DBI::dbConnect(RMySQL::MySQL(),
                               user=compoundDbUser,
                               password=compoundDbPass,
                               dbname=compoundDbName,
                               host=compoundDbHost,
                               pass=compoundDbPass)
  }

  #############################################
  # Add all the annotation results
  #############################################
  message('Add Metfrag results')
  addMetFragResults(metfrag_resultPth, con)
  message('Add Sirius CSI:FingerID results')
  addSiriusResults(sirius_csi_resultPth, con, con_comp=con_comp)
  message('Add Probmetab results')
  addProbmetabResults(probmetab_resultPth, con, con_comp=con_comp)
  message('Add Generic MS1 Lookup results')
  addGenericMS1LookupResults(ms1_lookup_resultPth,
                             ms1_lookup_dbSource,
                             ms1_lookup_checkAdducts,
                             ms1_lookup_keepAdducts,
                             con, con_comp=con_comp)
  # Add lipidsearch result TODO
  # Add mzcloud result TODO

  ####################################
  # Update the metab_compound
  ####################################
  # We now have inchikeys for all the annotations. For each inchikey we get all the associated informations
  message('Update the compound table')
  metab_compounds <- DBI::dbGetQuery(con, 'SELECT * FROM metab_compound')
  metab_compounds <- updateColName(metab_compounds, 'inchikey_id', 'inchikey')

  # Get all details for inchikeys (using local compoundDB - currently we only have the inchikeys)
  metab_compounds_m <- inchikeySelect(metab_compounds$inchikey, con = con_comp, all=TRUE)

  # Add anything that is not found (e.g. will be a few cases where pubchem won't have inchikey but the spectral matching result will have)
  missing_comps <- metab_compounds[!metab_compounds$inchikey %in% metab_compounds_m$inchikey,]
  missing_comps <- cbind(missing_comps, splitInchis(missing_comps$inchikey))
  colnames(missing_comps)[colnames(missing_comps)=='pubchem_id'] = 'pubchem_cids'

  keep_cols = c('inchikey', 'inchikey1', 'inchikey2', 'inchikey3', 'pubchem_cids',
                'exact_mass', 'molecular_formula', 'name')
  missing_comps <- missing_comps[,keep_cols]
  findx <- sapply(missing_comps, is.factor)
  missing_comps[findx] <- lapply(missing_comps[findx], as.character)
  missing_comps$exact_mass <- as.numeric(missing_comps$exact_mass)
  metab_compounds_m <- dplyr::bind_rows(metab_compounds_m, missing_comps)

  # remove the old table and add the modified table
  DBI::dbExecute(con, 'DROP TABLE metab_compound')
  DBI::dbWriteTable(conn=con, name='metab_compound', value=metab_compounds_m, row.names=FALSE)

  # Calculate final score and ranks
  # Build up the table by looping through xcms group
  message('Combine annotations')
  c_peak_groups <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups')
  
  DBI::dbDisconnect(con)
  combined_scores <- plyr::ddply(c_peak_groups, ~grpid, combineScoresGrp, weights=weights, sqlitePth=sqlitePth)
  
  message('Add combined scores to database')
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlitePth)
  DBI::dbWriteTable(conn=con, name='combined_annotations', value=combined_scores, row.names=FALSE)
  
  message('Get summary output')
  if (summaryOutput){
     DBI::dbDisconnect(con)
     con <- DBI::dbConnect(RSQLite::SQLite(), sqlitePth)
     anns <- getAnnotationSummary(con)

     DBI::dbDisconnect(con)
     DBI::dbDisconnect(con_comp)
     return(anns)
  }else{
     return(NULL)
  }

}

connect2db <- function(pth,type='sqlite',user=NA,pass=NA,dbname=NA,host=NA,port=NA){

  if (type=='sqlite'){
    con <- DBI::dbConnect(RSQLite::SQLite(), pth)
  }else if (type=='postgres'){
    con <- DBI::dbConnect(RPostgres::Postgres(),
                          user=user,
                          password=pass,
                          dbname=dbname,
                          host=host,
                          port=port,
                          pass=pass)
  }else if (type=='mysql'){
    con <- DBI::dbConnect(RMySQL::MySQL(),
                          user=user,
                          password=pass,
                          dbname=dbname,
                          host=host,
                          port=port,
                          pass=pass)
  }
  return(con)
}


addMetFragResults <- function(metfrag_resultPth, con, con_comp){
  # Add metfrag details
  if (!is.na(metfrag_resultPth) & file.exists(metfrag_resultPth) & length(count.fields(metfrag_resultPth, sep = "\t")>1)){
    if length(count.fields(metfrag_resultPth, sep = "\t"))
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



addSiriusResults <- function(sirius_csi_resultPth, con, con_comp=NULL){
  # Add sirius details
  if (!is.na(sirius_csi_resultPth) & file.exists(sirius_csi_resultPth) & length(count.fields(sirius_csi_resultPth, sep = "\t")>1)){
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
    if ('Score' %in% colnames(scores)){
      score <- abs(scores$Score)
    }else{
      score <- abs(scores$score)
    }
    bounded_score <- negMinMax(score)
  }
  scores$bounded_score <- bounded_score
  return(scores)

}

negMinMax <- function(x){
  xn <- (x-min(x))/(max(x)-min(x))
  return(abs(xn-1))

}



addProbmetabResults <- function(probmetab_resultPth, con, con_comp){
  if (!is.na(probmetab_resultPth) && file.exists(probmetab_resultPth)){
    addProbmetab(probmetab_resultPth, con)

    # Fetch in chunks
    kegg_cids <- DBI::dbGetQuery(con, "SELECT DISTINCT mpc FROM probmetab_results")

    kegg_cid_str <- vecToStringInStmt(kegg_cids$mpc)
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

adductCheckMS1Lookup <- function(row, keepAdducts, cameraAdducts){

  if ((row$adduct %in% keepAdducts) || (row$adduct %in% cameraAdducts[cameraAdducts$grpid==row$grpid,]$name)){
    return(row)
  }else{
    return(NULL)
  }

}

addGenericMS1LookupResults <- function(ms1_lookup_resultPth, ms1_lookup_dbSource, ms1_lookup_checkAdducts,
                                       ms1_lookup_keepAdducts, con, con_comp){
  if (!is.na(ms1_lookup_resultPth) & file.exists(ms1_lookup_resultPth) & length(count.fields(ms1_lookup_resultPth, sep = "\t")>1)){
    # Read in table

    ms1_lookup_result <- utils::read.table(ms1_lookup_resultPth,  header = TRUE, sep='\t', stringsAsFactors = FALSE,  comment.char = "",
                                           quote=NULL)

    if (!'grpid' %in% ms1_lookup_result){
      nmap <- DBI::dbGetQuery(con, 'SELECT grpid, grp_name FROM c_peak_groups')
      ms1_lookup_result <- merge(ms1_lookup_result, nmap, by.x='name', by.y='grp_name')
    }


    ms1_lookup_result$ms1_lookup_id <- 1:nrow(ms1_lookup_result)


    # Get the relevant adduct details from the CAMERA results (if we want to check if the adducts used
    # for the MS1 lookup are compatible)
    if(ms1_lookup_checkAdducts){
      cameraAdducts <- DBI::dbGetQuery(con, 'SELECT aa.grpid, ar.name FROM adduct_annotations
                                       AS aa LEFT JOIN adduct_rules AS ar ON aa.rule_id=ar.rule_id')
      # for each row of ms1_lookup_results, check if adduct is present adductKeeps,
      # if it is then keep, if not check if adduct present in cameraAdducts, if not the discard
    }else{
      cameraAdducts <- data.frame(matrix(ncol=2, nrow=0))
      colnames(cameraAdducts) <- c('grpid', 'adduct')
    }



    # Check if the adducts are from the "safe" list in keepAdducts, or if they match with camera adduct
    # annotations

    if (ms1_lookup_checkAdducts || (!is.null(ms1_lookup_keepAdducts) || !is.na(ms1_lookup_keepAdducts))){
      ms1_lookup_resultFiltered <- plyr::ddply(ms1_lookup_result , ~ ms1_lookup_id, adductCheckMS1Lookup,
                                               keepAdducts=ms1_lookup_keepAdducts, cameraAdducts=cameraAdducts)

    }else{
      ms1_lookup_resultFiltered <- ms1_lookup_result
    }

    # Get the compounds ids and relevant inchikeys
    compound_ids <- unique(ms1_lookup_resultFiltered$compound_id)
    compound_ids_str <- vecToStringInStmt(compound_ids)

    if (ms1_lookup_dbSource=='hmdb'){
      tableNm <- 'hmdb'
      columnNm <- 'hmdb_id'
    }else if (ms1_lookup_dbSource=='kegg'){
      tableNm <- 'kegg'
      columnNm <- 'kegg_cid'
    }else if (ms1_lookup_dbSource=='pubchem'){
      tableNm <- 'pubchem'
      columnNm <- 'pubchem_cid'
    }
    # Search the compound database for inchikeys
    compDetails = DBI::dbGetQuery(con_comp,
                                  sprintf("SELECT inchikey, %s FROM %s WHERE %s IN (%s)", columnNm, tableNm, columnNm, compound_ids_str))
    # save to database as reference
    DBI::dbWriteTable(con, name=tableNm, value=compDetails, row.names=FALSE)

    # update the results
    ms1_lookup_resultFiltered <- merge(ms1_lookup_resultFiltered, compDetails, by.x='compound_id', by.y=columnNm, all.x=TRUE)


    # Default score of 1 given for all results
    ms1_lookup_resultFiltered$ms1_lookup_score <- 1

    # Add the filtered results to the table
    DBI::dbWriteTable(con, name='ms1_lookup_results', value=ms1_lookup_resultFiltered, row.names=FALSE)

    # Check if any of the same inchikeys have been added
    new_compDetails <- insertNewInchikeys(con, con_comp, compDetails)

  }


}


getAnnotationSummary <- function(con){
  meta_cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(l_s_peak_meta)')
  if ("pid" %in% meta_cn){
    pid_id <- 'pid'
  }else{
    pid_id <- 'id'
  }

  cpg_cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(c_peak_group)')

  if ("adduct" %in% cpg_cn){
    camera_adduct <- 'cpg.adduct AS camera_adducts,'
  }else{
    camera_adduct <- ''
  }

  sql_stmt <- sprintf('
                      SELECT
                      cpg.grpid,
                      cpg.grp_name,
                      cpg.mz,
                      cpg.rt,
                      %s
                      mc.inchikey,
                      mc.inchi,
                      mc.inchikey,
                      mc.inchikey1,
                      mc.inchikey2,
                      mc.inchikey3,
                      mc.name,
                      mc.exact_mass,
                      mc.molecular_formula,
                      mc.pubchem_cids,
                      mc.kegg_cids,
                      mc.kegg_brite,
                      mc.kegg_drugs,
                      mc.hmdb_ids,
                      mc.hmdb_bio_custom_flag,
                      mc.hmdb_drug_flag,
                      mc.biosim_max_count,
                      mc.biosim_hmdb_ids,
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
                      ca.ms1_lookup_score,
                      ca.ms1_lookup_wscore,
                      mc.biosim_max_score,
                      ca.biosim_wscore,
                      ca.wscore,
                      ca.rank,
                      ca.adduct_overall
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
                      ', camera_adduct, pid_id)



  return(DBI::dbGetQuery(conn=con, sql_stmt))


}

combineRow <- function(row){
  row[row==0] = ''
  return(paste(unique(unique(row[row != ''])), collapse = ','))
}


combineScoresGrp <- function(c_peak_group, weights, sqlitePth){
  
  # connect to database 
  con <- connect2db(sqlitePth)
  
  # get sirius data
  grpid <- unique(c_peak_group$grpid)

  if(DBI::dbExistsTable(con, 'sirius_csifingerid_results')){
    # check if adduct column present
    sirius_cols <- DBI::dbGetQuery(con, 'PRAGMA table_info(sirius_csifingerid_results)')

    if ('adduct' %in% sirius_cols$name){
      adduct_colstr <- 's.adduct AS sirius_adduct'
    }else{
      adduct_colstr <- "'' AS sirius_adduct"
    }

    sirius <- DBI::dbGetQuery(con, sprintf('SELECT
                                           s.rowid AS sirius_id,
                                           s.bounded_score AS sirius_score,
                                           %s,
                                           mc.inchikey
                                           FROM metab_compound AS mc
                                           LEFT JOIN
                                           sirius_csifingerid_results AS s
                                           ON s.inchikey2D = mc.inchikey1 WHERE s.grpid= :grpid', adduct_colstr),
                              params=list('grpid'=grpid))
    # can have multiple hits for each inchikey (e.g. multiple scans)
    sirius <- plyr::ddply(sirius, ~inchikey, getBestScore, 'sirius_score')


    sirius$sirius_wscore <- as.numeric(sirius$sirius_score)*weights$sirius_csifingerid

  }else{
    sirius <- data.frame(inchikey=NA, sirius_score=NA, sirius_wscore=NA, sirius_adduct=NA)
  }


  if(DBI::dbExistsTable(con, 'metfrag_results')){
    metfrag_cols <- DBI::dbGetQuery(con, 'PRAGMA table_info(metfrag_results)')
    if ('adduct' %in% metfrag_cols$name){
      adduct_colstr <- 'm.adduct AS metfrag_adduct'
    }else{
      adduct_colstr <- "'' AS metfrag_adduct"
    }

    # get metfrag
    metfrag <- DBI::dbGetQuery(con, sprintf('SELECT
                                            m.rowid as metfrag_id,
                                            m.Score AS metfrag_score,
                                            %s,
                                            mc.inchikey
                                            FROM metab_compound AS mc
                                            LEFT JOIN
                                            metfrag_results AS m
                                            ON m.InChIKey = mc.inchikey WHERE m.grpid=:grpid', adduct_colstr),
                               params=list('grpid'=grpid))
    # can have multiple hits for each inchikey (e.g. multiple scans)
    metfrag <- plyr::ddply(metfrag, ~inchikey, getBestScore, 'metfrag_score')

    metfrag$metfrag_wscore <- as.numeric(metfrag$metfrag_score)*weights$metfrag

  }else{
    metfrag <- data.frame(inchikey=NA, metfrag_score=NA, metfrag_wscore=NA, metfrag_adduct=NA)

  }

  meta_cn <- DBI::dbGetQuery(con, 'PRAGMA table_info(l_s_peak_meta)')
  if ("pid" %in% meta_cn){
    pid_id <- 'pid'
  }else{
    pid_id <- 'id'
  }

  # get spectral matching
  sm <- DBI::dbGetQuery(con, sprintf('SELECT sm.lpid AS sm_lpid,
                                     sm.library_precursor_type AS sm_adduct,
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

  if(DBI::dbExistsTable(con, 'ms1_lookup_results')){
    ms1_lookup <- DBI::dbGetQuery(con, 'SELECT
                                  m.ms1_lookup_id,
                                  m.ms1_lookup_score,
                                  m.adduct AS ms1_lookup_adduct,
                                  mc.inchikey
                                  FROM metab_compound AS mc
                                  LEFT JOIN
                                  ms1_lookup_results AS m
                                  ON m.inchikey = mc.inchikey WHERE m.grpid=:grpid',
                                  params=list('grpid'=grpid))

    # can have multiple hits for each inchikey (e.g. multiple scans)
    ms1_lookup <- plyr::ddply(ms1_lookup, ~inchikey, getBestScore, 'ms1_lookup_score')
    ms1_lookup$ms1_lookup_wscore <- as.numeric(ms1_lookup$ms1_lookup_score)*weights$ms1_lookup
  }else{
    ms1_lookup <- data.frame(inchikey=NA, ms1_lookup_score=NA, ms1_lookup_wscore=NA, ms1_lookup_adduct=NA)
  }


  score_list <- list(sirius, metfrag, sm, probmetab, ms1_lookup)

  # combine all
  combined_df <- Reduce(function(...) merge(..., all=TRUE, by='inchikey'), score_list)

  # Get biosim score (if available)
  biosim_score <- DBI::dbGetQuery(con, 'SELECT inchikey, biosim_max_score FROM metab_compound')
  combined_df <- merge(combined_df, biosim_score, by='inchikey')
  combined_df$biosim_wscore <- combined_df$biosim_max_score * weights$biosim

  combined_df <- combined_df[!is.na(combined_df$inchikey),]
  combined_df[is.na(combined_df)] <- 0
  combined_df$wscore <- rowSums(combined_df[, c('sirius_wscore', 'metfrag_wscore', 'sm_wscore',
                                                'probmetab_wscore', 'probmetab_wscore',
                                                'ms1_lookup_wscore', 'biosim_wscore')])

  # Order
  combined_df <- combined_df[order(combined_df$wscore, decreasing=TRUE),]



  # in most cases the adduct will be the same but there is a chance a different naming was used
  # or that there are adduct types that have give the same neutral mass
  combined_df$adduct_overall <- apply(combined_df[,c('sirius_adduct',
                                                     'metfrag_adduct', 'ms1_lookup_adduct',
                                                     'sm_adduct')], 1, combineRow)

  if (nrow(combined_df)==0){
    return(NULL)
  }else{

    # Add rank
    combined_df$rank <- as.numeric(factor(rank(-combined_df$wscore)))
  }
                        
  DBI::dbDisconnect(con)

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

  return(DBI::dbGetQuery(con, sql_stmt))
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
