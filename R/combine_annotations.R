#' @title Combine Annotations
#'
#' @description
#'
#' @examples
#'
#' @export

combineAnnotations <- function(sqlitePth,
                               metfrag_resultPth=NA,
                               sirius_csi_resultPth=NA,
                               probmetab_resultPth=NA,
                               lipidsearch_resultPth=NA,
                               mzcloud_resultPth=NA,
                               xset,
                               weights=c(0.4, 0.2, 0.2, 0.05, 0, 0, 0.05),
                               avFrag="av_all"){

  metfrag_resultPth <- './inst/extdata/external_annotations/metfrag.tsv'
  sirius_csi_resultPth <- './inst/extdata/external_annotations/sirus_csifingerid.tsv'
  probmetab_resultPth <- './inst/extdata/external_annotations/probmetab.tsv'
  # sm, metfrag, sirius, probmetab, lipidsearch, mzcloud, bs
  con <- DBI::dbConnect(RSQLite::SQLite(), db_path)

  # Add metfrag details
  if (!is.na(metfrag_resultPth) & file.exists(metfrag_resultPth)){
    DBI::dbWriteTable(conn=con, name='metfrag_results', value=metfrag_resultPth, sep='\t', header=T)
    #  any new inchikeys
    dbExecute(con,
      'INSERT INTO metab_compound (inchikey_id)
          SELECT DISTINCT m.InChIKey
      FROM metfrag_results AS m
          LEFT JOIN
              metab_compound AS c ON m.InChIKey = c.inchikey_id
      WHERE c.inchikey_id IS NULL')
  }

  # Add metfrag details
  if (!is.na(sirius_csi_resultPth) & file.exists(sirius_csi_resultPth)){
    DBI::dbWriteTable(conn=con, name='sirius_csifingerid_results', value=sirius_csi_resultPth, sep='\t', header=T)
    # any new inchikeys
    # get all pubchemids

    metfrag_inchikeys  <- dbExecute(con,
      'INSERT INTO metab_compound (inchikey_id)
          SELECT DISTINCT m.InChIKey
      FROM metfrag_results AS m
          LEFT JOIN
              metab_compound AS c ON m.InChIKey = c.inchikey_id
      WHERE c.inchikey_id IS NULL')
  }





  # Add any new inchikeys

  # Add probmetab results
  addProbmetab(probmetab_resultPth, xset,  con)

  # Add any new inchikeys



  # Add lipidsearch result
  # Add any new inchikeys

  # Add mzcloud result
  # Add any new inchikeys

  # Check new compound entries - and get relevant information from pubchem

  # Calculate compound "biological" similarity score

  # Calculate final score and rank

}







addProbmetab <- function(pth, xset, con){
  if (!is.null(pth)){

    df <- read.table(pth,  header = TRUE, sep='\t', stringsAsFactors = FALSE,  comment.char = "")
    df$grp_id <- match(df$name, xcms::groupnames(xset))
    start <- T
    for (i in 1:nrow(df)){

      x <- df[i,]

      if(is.na(x$proba) | x$proba =='NA'){
        next
      }

      mpc <- stringr::str_split(x$mpc, ';')
      proba <- stringr::str_split(x$proba, ';')

      for (j in 1:length(mpc[[1]])){

        row <-  c(x$grp_id, x$propmz, mpc[[1]][j], proba[[1]][j])

        if (start){
          df_out <- data.frame(t(row), stringsAsFactors=F)
          start <- F
        }else{
          df_out <- data.frame(rbind(df_out, row), stringsAsFactors=F)
        }
      }

    }

    colnames(df_out) <- c('grp_id', 'propmz', 'mpc', 'proba')
    DBI::dbWriteTable(con, name='probmetab_results', value=df_out, row.names=FALSE)

  }
}



