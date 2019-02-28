#' @title Combine Annotations
#'
#' @description
#'  Combine the annotation results from msPurity spectral matching, MetFrag, Sirius CSI:FingerID and probmetab
#'  based on weighted scores for each technique aligning each annotation by inchikey and XCMS grouped feature.
#'
#' @param sqlitePth character;
#' @param metfrag_resultPth character;
#' @param sirius_csi_resultPth character;
#' @param probmetab_resultPth character;
#' @param weights list;
#' @param silentRestErrors boolean;
#'
#' @examples
#' metfrag_resultPth <- system.file("extdata", "external_annotations", "metfrag.tsv", package="msPurity")
#' sirius_csi_resultPth <- system.file("extdata", "external_annotations", "sirus_csifingerid.tsv", package="msPurity")
#' probmetab_resultPth <- system.file("extdata", "external_annotations", "probmetab.tsv", package="msPurity")
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths)
#' pa <- frag4feature(pa, xset)
#' pa <- averageAllFragSpectra(pa, minfrac = 1)
#' dbpth <- create_database(pa, xset)
#' result <- spectral_matching(dbpth, spectra_type_q="av_all")
#' combineAnnotations(dbpth, metfrag_resultPth, sirius_csi_resultPth, probmetab_resultPth)
#' @return purityA object with slots for fragmentation-XCMS links
#' @export
combineAnnotations <- function(sqlitePth,
                               metfrag_resultPth=NA,
                               sirius_csi_resultPth=NA,
                               probmetab_resultPth=NA,
                               weights=list('sm'=0.4,'metfrag'=0.2,'sirius_csifingerid'= 0.2,
                                            'probmetab'=0.1
                               ),
                               silentRestErrors=FALSE
                               ){

  # sm, metfrag, sirius, probmetab, lipidsearch, mzcloud, bs
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlitePth)

  #############################################
  # Add all the annotation results
  #############################################
  message('Add Metfrag results')
  addMetFragResults(metfrag_resultPth, con)
  message('Add Sirius CSI:FingerID results')
  addSiriusResults(sirius_csi_resultPth, con, silentRestErrors=silentRestErrors)
  message('Add Probmetab results')
  addProbmetabResults(probmetab_resultPth, con, silentRestErrors=silentRestErrors)
  # Add lipidsearch result TODO
  # Add mzcloud result TODO


  ####################################
  # Update the metab_compound (and associated compound source tables e.g. KEGG, HMDB)
  ####################################
  # We now have inchikeys for all the annotations. For each inchikey we get all the associated
  # information from KEGG, Pubchem, LipidMaps and some additional IDS
  message('Update the compound table')
  # Check new metab_compound entries - and get relevant information from
  metab_compounds <- RSQLite::dbGetQuery(con, 'SELECT * FROM metab_compound')

  metab_compounds <- updateColName(metab_compounds, 'inchikey_id', 'inchikey')

  # pubchem
  metab_compounds_m <- addPubchemToMetabCompound(metab_compounds, con, silentRestErrors=silentRestErrors)

  # KEGG
  metab_compounds_m <- addKeggToMetabCompound(metab_compounds_m, con, silentRestErrors=silentRestErrors)

  # LipidMaps
  # Takes too long
  #lipidmaps_ids <- webchem::cts_convert(metab_compounds_m$inchikey, from = 'inchikey', to='lipidmaps')

  # HMDB
  # Takes too long
  #hmdb_ids <- webchem::cts_convert(metab_compounds_m$inchikey, from = 'inchikey', to='Human Metabolome Database')

  # PAN pesticide database
  # Takes too long
  #pan_pesticide <- webchem::pan_query(metab_compounds_m$name  , match='best')

  # Alan Wood's Compendium of Pesticide Common Names
  # Takes too long
  #aw_pesticide <- webchem::aw_query(metab_compounds_m$name , type='commonname')

  #message('Create biological similarity')
  # Calculate compound "biological" similarity score TODO
  #metab_compounds_m$biological_similarity <- 0
  #metab_compounds_m$biological_similarity_top_match <- NA

  # addd New database with all the new information
  inchisplit <- data.frame(stringr::str_split_fixed(metab_compounds_m$inchikey, '-', 3))
  colnames(inchisplit) <- c('inchikey1', 'inchikey2', 'inchikey3')
  metab_compounds_m <- cbind(metab_compounds_m, inchisplit)
  DBI::dbExecute(con, 'DROP TABLE metab_compound')
  DBI::dbWriteTable(conn=con, name='metab_compound', value=metab_compounds_m, row.names=FALSE)

  # Calculate final score and rank
  # Build up the table by looping through xcms group
  message('Combine annotations')
  c_peak_groups <- DBI::dbGetQuery(con, 'SELECT * FROM c_peak_groups')
  combined_scores <- plyr::ddply(c_peak_groups, ~grpid,combineScoresGrp, weights=weights, con=con)

  DBI::dbWriteTable(conn=con, name='combined_annotations', value=combined_scores, row.names=FALSE)

  return(combined_scores)



}

combineScoresGrp <- function(c_peak_group, weights, con){

  # get sirius data
  grpid <- unique(c_peak_group$grpid)

  sirius <- DBI::dbGetQuery(con, 'SELECT
                            s.rowid AS sirius_id,
                            s.bounded_score AS sirius_score,
                            mc.inchikey
                            FROM metab_compound AS mc
                            LEFT JOIN
                            sirius_csifingerid_results AS s
                            ON s.InChIkey2D = mc.inchikey1 WHERE s.grpid= :grpid', params=list('grpid'=grpid))
  # can have multiple hits for each inchikey (e.g. multiple scans)
  sirius <- plyr::ddply(sirius, ~inchikey, getBestScore, 'sirius_score')


  sirius$sirius_wscore <- as.numeric(sirius$sirius_score)*weights$sirius_csifingerid

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

  # get spectral matching
  sm <- DBI::dbGetQuery(con, 'SELECT sm.lid AS sm_lid,
                        sm.mid AS sm_mid,
                        sm.score AS sm_score,
                        mc.inchikey
                        FROM metab_compound AS mc
                        LEFT JOIN
                        library_spectra_meta AS l ON l.inchikey_id = mc.inchikey
                        LEFT JOIN
                        matches AS sm ON l.lid = sm.lid
                        WHERE sm.grpid = :grpid', params=list('grpid'=grpid))

  # can have multiple hits for each inchikey (e.g. multiple scans and library spectra for each inchikey)
  # we just take the best hit
  sm <- plyr::ddply(sm, ~inchikey, getBestScore, 'sm_score')
  sm$sm_wscore <- as.numeric(sm$sm_score)*weights$sm


  # get probmetab
  probmetab <- DBI::dbGetQuery(con, paste('SELECT p.rowid as probmetab_id,
                                          p.proba AS probmetab_score,
                                          mc.inchikey
                                          FROM metab_compound AS mc
                                          LEFT JOIN
                                          kegg as k ON k.inchikey=mc.inchikey
                                          LEFT JOIN
                                          probmetab_results AS p ON p.mpc = k.kegg_id
                                          WHERE p.grpid= ', grpid, sep=' '))
  # can have multiple hits for each inchikey (e.g. multiple scans)
  probmetab <- plyr::ddply(probmetab, ~inchikey, getBestScore, 'probmetab_score')

  probmetab$probmetab_wscore <- as.numeric(probmetab$probmetab_score)*weights$probmetab

  score_list <- list(sirius, metfrag, sm, probmetab)

  # combine al
  combined_df <- Reduce(function(d1, d2, d3, d4) merge(d1, d2, by = "inchikey", all = TRUE), score_list)

  combined_df[is.na(combined_df)] <- 0
  combined_df$wscore <- rowSums(combined_df[, c('sirius_wscore', 'metfrag_wscore', 'metfrag_wscore', 'sm_wscore', 'probmetab_wscore')])

  # add rank
  combined_df <- combined_df[order(combined_df$wscore, decreasing=TRUE),]

  combined_df$rank[order(combined_df$wscore, decreasing = FALSE)]

  return(combined_df)
}


updateCol<- function(df, old, new, type='bothna'){

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

addPubchemToMetabCompound <- function(metab_compounds, con, silentRestErrors){
  # Pubchem
  pubchem_details <- plyr::adply(metab_compounds$inchikey, 1,getPubchemDetailsExtended, silentRestErrors)


  DBI::dbWriteTable(con, name='pubchem', value=pubchem_details[,c('inchikey', 'cid')], row.names=FALSE)

  # filter out any duplicate pubchem entries (keep the first hit for the metab_compound database)
  pubchem_details <- pubchem_details[!duplicated(pubchem_details$inchikey),]

  metab_compounds_m <- merge(metab_compounds, pubchem_details[,c('inchikey','name', 'other_names', 'smiles', 'molecular_formula',
                                                                 'molecular_weight', 'exact_mass', 'iupac_prefered')],
                             by = "inchikey", suffixes = c('', '_pubchem'), all = TRUE)

  metab_compounds_m <- updateCol(metab_compounds_m, 'name', 'iupac_prefered')
  metab_compounds_m <- updateCol(metab_compounds_m, 'name', 'name_pubchem')

  metab_compounds_m <- updateCol(metab_compounds_m, 'other_names', 'other_names_pubchem')
  metab_compounds_m <- updateCol(metab_compounds_m, 'smiles', 'smiles_pubchem')
  metab_compounds_m <- updateCol(metab_compounds_m, 'molecular_formula', 'molecular_formula_pubchem')
  metab_compounds_m <- updateCol(metab_compounds_m, 'molecular_weight', 'molecular_weight_pubchem')
  metab_compounds_m <- updateCol(metab_compounds_m, 'exact_mass', 'exact_mass_pubchem')

  return(metab_compounds_m)
}

addKeggToMetabCompound <- function(metab_compounds, con, silentRestErrors){
  # KEGG (this will give both Drug ID and Compound ID ) - this will give the very basic
  # info if the KEGG compound is a drug or not
  keggids <- plyr::adply(metab_compounds$inchikey,1, getKeggFromInchi, silentRestErrors=silentRestErrors)
  keggids <- keggids[!keggids$kegg_id=='unknown',]


  DBI::dbWriteTable(con, name='kegg', value=keggids[,c('inchikey', 'kegg_id')], row.names=FALSE)
  #


  kegg_details <- plyr::ddply(keggids, ~inchikey, function(x){
    df <- data.frame(t(sapply(x$kegg_id,  getKeggDetail)))
    if (all(is.na((df$kegg_did)))){
      df$drug <- NA
    }else{
      dids <- df$kegg_did[!is.na(df$kegg_did)]
      df$drug <- paste(dids, collapse=',')
    }
    return(df)

  })

  # we don't want to include the drug details
  kegg_details <- kegg_details[is.na(kegg_details$kegg_did),]

  # if more than one kegg compound per inchi we just record the first one (all others)
  # are recorded in the kegg-inchi table
  kegg_details <- kegg_details[!duplicated(kegg_details$inchikey),]

  colnames(kegg_details)[colnames(kegg_details)=='name'] = 'kegg_name'
  colnames(kegg_details)[colnames(kegg_details)=='other_names'] = 'kegg_other_names'

  metab_compounds_m <- merge(metab_compounds, kegg_details[,c('inchikey','drug', 'brite1', 'brite2', 'kegg_name', 'kegg_other_names')],
                             by= "inchikey", all = TRUE)

  metab_compounds_m$kegg_name <- as.character(metab_compounds_m$kegg_name)
  metab_compounds_m$kegg_other_names <- as.character(metab_compounds_m$kegg_other_names)


  metab_compounds_m <- updateCol(metab_compounds_m, 'name', 'kegg_name')
  metab_compounds_m <- updateCol(metab_compounds_m, 'other_names', 'kegg_other_names')


  return(metab_compounds_m)


}

addMetFragResults <- function(metfrag_resultPth, con){
  # Add metfrag details
  if (!is.na(metfrag_resultPth) & file.exists(metfrag_resultPth)){
    DBI::dbWriteTable(conn=con, name='metfrag_results', value=metfrag_resultPth, sep='\t', header=T)
    #  any new inchikeys
    DBI::dbExecute(con,
                   'INSERT INTO metab_compound (inchikey_id)
                   SELECT DISTINCT m.InChIKey
                   FROM metfrag_results AS m
                   LEFT JOIN
                   metab_compound AS c ON m.InChIKey = c.inchikey_id
                   WHERE c.inchikey_id IS NULL')
  }

  }

addSiriusResults <- function(sirius_csi_resultPth, con, silentRestErrors){
  # Add sirius details
  if (!is.na(sirius_csi_resultPth) & file.exists(sirius_csi_resultPth)){
    DBI::dbWriteTable(conn=con, name='sirius_csifingerid_results', value=sirius_csi_resultPth, sep='\t', header=T, row.names=T, nrows = 4)

    # Fetch in chunks
    inchikey2ds <- DBI::dbGetQuery(con, "SELECT DISTINCT InChIKey2D FROM sirius_csifingerid_results")



    pubchemDetails <- plyr::adply(inchikey2ds$InChIkey2D, 1,getPubchemDetails, silentRestErrors=silentRestErrors)

    # filter out pubchem details we already have. And add new
    # Shouldnt be required as we haven't added any pubchem yet
    # rs <- RSQLite::dbSendQuery(con, sprintf('SELECT * FROM pubchem WHERE "cid" IN (%s)', paste(pubchemDetails$cid, collapse=',')))
    # matching_cid <- DBI::dbFetch(rs)
    # DBI::dbClearResult(rs)
    # pubchemDetails[pubchemDetails$cid %in% matching_cid$cid,]
    #custom_dbWriteTable(name_pk = 'cid', fks=NA,
    #                    table_name = 'pubchem', df=pubchemDetails, con=con)
    inchi_str <- paste("'", unique(pubchemDetails$inchikey), "'", collapse=",", sep='')
    # filter out pubchem details we already have.
    rs <- RSQLite::dbSendQuery(con, sprintf('SELECT inchikey_id FROM metab_compound WHERE "inchikey_id" IN (%s)', inchi_str))
    matching_inchi <- DBI::dbFetch(rs)
    DBI::dbClearResult(rs)
    new_pubchemDetails <- pubchemDetails[!pubchemDetails$inchikey %in% matching_inchi$inchikey_id,]

    new_pubchemDetails = new_pubchemDetails[!duplicated(new_pubchemDetails$inchikey),]

    if (nrow(new_pubchemDetails)>0){
      new_inchi_str <- paste("('", unique(new_pubchemDetails$inchikey), "')", collapse=",", sep='')
      insert_stmt <- sprintf('INSERT INTO metab_compound (inchikey_id) VALUES %s', new_inchi_str)
      # Add new full inchikeys to metab_compound
      rs <- DBI::dbExecute(con, insert_stmt)

    }
    ##################
    # calculate Sirius score
    ##################
    # Loop through sqlite database by UID.
    uids <- DBI::dbGetQuery(con, "SELECT DISTINCT UID FROM sirius_csifingerid_results")

    rs <- DBI::dbSendQuery(con, "SELECT rowid, Score FROM sirius_csifingerid_results WHERE UID = ?")

    bounded_score <- plyr::adply(uids$UID, 1, getBoundedSiriusScore, rs=rs)
    DBI::dbClearResult(rs)
    # caculate minmax normalised value (save)
    # Add new column
    DBI::dbExecute(con, "ALTER TABLE sirius_csifingerid_results ADD COLUMN bounded_score text")

    DBI::dbExecute(con, "UPDATE sirius_csifingerid_results SET bounded_score = :bounded_score  WHERE rowid = :rowid",
                   params=bounded_score[,c('rowid', 'bounded_score')])

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
addProbmetabResults <- function(probmetab_resultPth, con, silentRestErrors){
  if (!is.na(probmetab_resultPth) & file.exists(probmetab_resultPth)){
    addProbmetab(probmetab_resultPth, con)

    # Fetch in chunks
    kegg_cids <- DBI::dbGetQuery(con, "SELECT DISTINCT mpc FROM probmetab_results")


    # Should keep top level of Brite (e.g. )
    #kegg_details <- plyr::adply(kegg_cids, 1, getKeggDetail)

    inchikey_col <- plyr::adply(kegg_cids, 1, getInchiFromKeggCid, silentRestErrors=silentRestErrors)

    #kegg_alldetails <- data.frame(merge(kegg_details, inchikey_col[,-1]))
    inchi_str <- paste("'", inchikey_col$inchikey, "'",collapse = ',', sep='')
    sql_stmt <- sprintf('SELECT inchikey_id FROM metab_compound WHERE "inchikey_id" IN (%s)', inchi_str)

    matching_inchi <- RSQLite::dbGetQuery(con, sql_stmt)

    new_inchi <- inchikey_col[!inchikey_col$inchikey %in% matching_inchi$inchikey_id,]

    if (nrow(new_inchi)>0){
      new_inchi_str <- paste("('", unique(new_inchi$inchikey), "')", collapse=",", sep='')
      insert_stmt <- sprintf('INSERT INTO metab_compound (inchikey_id) VALUES %s', new_inchi_str)
      # Add new full inchikeys to metab_compound
      DBI::dbExecute(con, insert_stmt)

    }

    # need to add "unknown" compounds


  }

}

getKeggDetail <-  function(kegg_id){
  kegg_results <- KEGGREST::keggGet(kegg_id)
  kegg_result <- kegg_results[[1]]

  # cid
  kegg_cid <- unname(unlist(kegg_result$ENTRY['Compound']))
  kegg_did <- unname(unlist(kegg_result$ENTRY['Drug']))

  # name
  name <- as.character(kegg_result$NAME[1])

  # other names
  if (length(kegg_result$NAME)>1){
    other_names <- paste(kegg_result$NAME[2:length(kegg_result$NAME)], collapse=' ')
  }else{
    other_names <- NA
  }

  ending <- substr(name, nchar(name), nchar(name))

  if (ending==';'){
    name = substr(name,1,nchar(name)-1)
  }


  # formula
  mf <- kegg_result$FORMULA

  # exact mass
  exact_mass <- kegg_result$EXACT_MASS

  # Mol weight
  molw <- kegg_result$MOL_WEIGHT

  # Brite 1
  if (length(kegg_result$BRITE)>0){
    brite1 <- kegg_result$BRITE[1]
  }else{
    brite1 <- NA
  }

  if (length(kegg_result$BRITE)>1){
    brite2 <- kegg_result$BRITE[2]
  }else{
    brite2 <- NA
  }

  row <- c('kegg_cid'=kegg_cid, 'kegg_did'=kegg_did,'name'=name, 'other_names'=other_names, 'molecular_formula'=mf,
           'exact_mass'=exact_mass, 'molecular_weight'=molw, 'brite1'=brite1, 'brite2'=brite2)


  return(row)



}

getInchiFromKeggCid <- function(kegg_cid, silentRestErrors){
  # Get KEGG details

  inchikey_result  <- NULL
  attempt <- 1
  while( is.null(inchikey_result ) && attempt <= 3 ) {
    attempt <- attempt + 1
    rest_stmt <- sprintf("https://cts.fiehnlab.ucdavis.edu/rest/score/KEGG/%s/biological", kegg_cid)

    tryCatch({
      inchikey_result <-  jsonlite::fromJSON(rest_stmt)
    }, error = function(e) {
      if (silentRestErrors){
        message('REST-API ERROR,', e)
      }

    })

  }

  if (is.null(inchikey_result ) || (inchikey_result$result$score[1]==0)){
    inchikey = paste('UNKNOWN-KEGG-CID', kegg_cid, 'uid', uuid::UUIDgenerate(use.time = NA), sep='-')
  }else{
    inchikey = inchikey_result$result$InChIKey
  }

  return(data.frame(c('kegg_cid'=kegg_cid,'inchikey'=inchikey)))

}

getKeggFromInchi <- function(inchikey, silentRestErrors){
  pat <- "UNKNOWN-KEGG-CID-(C[0-9]+).*"
  if (grepl(pat, inchikey)){
    kegg_unknown <- sub(pat, "\\1", inchikey[grepl(pat, inchikey)])
    dfout <- data.frame('kegg_id'= kegg_unknown, 'inchikey'= inchikey)
    return(dfout)
  }


  inchikey# Get KEGG details
  kegg_ids = NA
  kegg_result  <- NULL
  attempt <- 1
  while( is.null(kegg_result ) && attempt <= 3 ) {
    attempt <- attempt + 1
    rest_stmt <- sprintf("https://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/KEGG/%s", inchikey)

    tryCatch({
      kegg_result <-  jsonlite::fromJSON(rest_stmt)
    }, error = function(e) {
      if (silentRestErrors){
        message('REST-API ERROR,', e)
      }
    })


    if (length(kegg_result$results[[1]])==0){
      kegg_ids = 'unknown'
    }else{
      kegg_ids = as.character(unlist(kegg_result$results[[1]]))
    }

    dfout <- data.frame('kegg_id'=kegg_ids, 'inchikey'= inchikey)
    return(dfout)

  }
}

addProbmetab <- function(pth, con){
  if (!is.null(pth)){

    df <- read.table(pth,  header = TRUE, sep='\t', stringsAsFactors = FALSE,  comment.char = "")

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


getPubchemDetails <- function(inchikey, silentRestErrors){


  compounds <- NULL
  attempt <- 1
  while( is.null(compounds) && attempt <= 3 ) {
    attempt <- attempt + 1
    parse_str <- sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/json", inchikey)

    tryCatch({
      compounds <-  jsonlite::fromJSON(parse_str)
    }, error = function(e) {
      if (silentRestErrors){
        message('REST-API ERROR,', e)
      }
    })
  }

  if (is.null(compounds)){
    return(NULL)
  }

  details <- plyr::ldply(compounds$PC_Compounds$props, getPubchemPropAll)

  cids <- compounds$PC_Compounds$id$id$cid
  details$cid <- cids


  inchisplit <- stringr::str_split_fixed(details$inchikey, "-", 3)
  colnames(inchisplit) <- c('inchikey1', 'inchikey2', 'inchikey3')
  alldetails <- cbind(details, inchisplit)

  return(details)
}


getPubchemDetailsExtended <- function(inchikey, silentRestErrors){
  details <- getPubchemDetails(inchikey, silentRestErrors=silentRestErrors)
  syn <- getSynonyms(details$cid)
  if (!is.null(syn)){
    details <- merge(syn, details)
  }
  details$name[is.na(details$name)] <- details$iupac_prefered[is.na(details$name)]

  return(details)
}

getSynonyms <- function(cids, silentRestErrors=FALSE){
  synonyms <- NULL
  attempt <- 1
  while( is.null(synonyms) && attempt <= 3 ) {
    attempt <- attempt + 1
    parse_str <- sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/json',
                         paste(cids, collapse=','))

    tryCatch({
      synonyms <-  jsonlite::fromJSON(parse_str)
    }, error = function(e) {
      if (silentRestErrors){
        message('REST-API ERROR,', e)
      }
    })
  }

  if (is.null(synonyms)){
    return(NULL)
  }

  top5df <- plyr::ldply(synonyms$InformationList$Information$Synonym, shortName)

  top5df$cid <- synonyms$InformationList$Information$CID
  colnames(top5df)[1:2] <- c('name', 'other_names')

  return(top5df)
}

shortName <- function(x){

  if(length(x)==0){
    return(c(NA, NA))
  }else if(length(x)<5){
    max_indx = length(x)
  }else{
    max_indx = 5
  }

  return(c(x[1], paste(x[2:max_indx], collapse=',')))

}

getPubchemPropAll <- function(x){
  return(c(  'inchikey' = getPubchemProp(x, 'InChIKey'),
             'exact_mass' = getPubchemProp(x, 'Mass', 'Exact', 'fval'),
             'molecular_formula' = getPubchemProp(x, 'Molecular Formula'),
             'molecular_weight' = getPubchemProp(x, 'Molecular Weight', 'fval'),
             'iupac_systematic' = getPubchemProp(x, 'IUPAC Name', 'Systematic'),
             'iupac_prefered' = getPubchemProp(x, 'IUPAC Name', 'Preferred'),
             'smiles' = getPubchemProp(x, 'SMILES', 'Canonical')
  ))
}

getPubchemProp <- function(x, label, name=NA, type='sval'){
  urn <- x$urn
  value <- x$value
  prop = NA

  if ((is.na(name)) & (label %in% urn$label)){
    prop <- value[,type][urn$label==label]
  }else if ((label %in% urn$label) & (name %in% urn$name) ){
    prop <- value[,type][(urn$label==label) & (urn$name==name)]
  }

  return(prop)

}




