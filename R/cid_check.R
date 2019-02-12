install.packages("httr")

#Require the package so you can use it
require("httr")

install.packages("jsonlite")

#Require the package so you can use it
require("jsonlite")


r <- NULL
attempt <- 1
while( is.null(r) && attempt <= 3 ) {
  attempt <- attempt + 1
  try(
    r <- extractPubChemData(inchikey)
  )
}

'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/679/synonyms/TXT'

extractCids <- function(inchikey){
  compounds <-  fromJSON(sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/%s/json', inchikey))

  details <- plyr::ldply(compounds$PC_Compounds$props, getPubchemPropAll)

  cids <- compounds$PC_Compounds$id$id$cid
  details$cid <- cids

  merge(getSynonyms(cids), details, by=cid)

  head(details)


  return(cids)
}

getSynonyms <- function(cids){
  synonyms <-  fromJSON(sprintf('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/%s/synonyms/json',
                                paste(cids, collapse=',')))

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

  return(c(x[1], paste(x[1:max_indx], collapse=',')))

}

getPubchemProp <- function(x, label, name=NA){
  urn <- x$urn
  value <- x$value
  prop = NA
  if (is.na(name)){
    if (label %in% urn$label){
      prop <- value$sval[urn$label==label]
    }
  }else{
    if ((label %in% urn$label) & (name %in% urn$name) ){
      prop <- value$sval[(urn$label==label) & (urn$name==name)]
    }
  }
  return(prop)

}

getPubchemPropAll <- function(x){
  return(c(  'inchikey' = getPubchemProp(x, 'InChIKey'),
             'exact_mass' = getPubchemProp(x, 'Mass', 'Exact'),
             'mf' = getPubchemProp(x, 'Molecular Formula'),
             'mw' = getPubchemProp(x, 'Molecular Weight'),
             'iupac_systematic' = getPubchemProp(x, 'IUPAC Name', 'Systematic'),
             'iupac_prefered' = getPubchemProp(x, 'IUPAC Name', 'Preferred')
  ))
}


