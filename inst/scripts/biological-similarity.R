#' #' @title Biological similarity based SVM model
#' #'
#' #' @description
#' #'  Search a query list of compounds against a list of biological compounds and determine how similar to
#' #'  biological structure (uses hmdb)
#' #'  https://www.bioconductor.org/packages/devel/bioc/vignettes/ChemmineR/inst/doc/ChemmineR.html
#' #'
#' #'  Biological compounds currently defaults to hmdb
#' #'
#' #' @examples
#' #' @export
#' #'
#' bioSim <- function(query_inchikeys){
#'
#' }
#'
#' bioSimSdfSingle <- function(querySdfPth,
#'                             model=NULL,
#'                             modelDistances=NULL,
#'                             popSd=NULL,
#'                             popMean=NULL,
#'                             librarySdfPth=NULL,
#'                             klim=0.5,
#'                             topn=10){
#'   # THIS DOES NOT WORK FOR LARGE MOLECULES - Always seems to crash!
#'   # get all sdf for molecules
#'   # this part should be outside of function for speed
#'   if ((is.null(model)) & (is.null(model_distances)){
#'     if(is.null(librarySdfPth){
#'       model <- system.file("extdata","models","hmdb_svm_bio_endo_model.RDS" package="msPurity")
#'       modelDistances <- system.file("extdata","models","hmdb_svm_bio_endo_distances.RDS" package="msPurity")
#'       librarySdfPth <- system.file("extdata","models","hmdb_bio_endo.sdf" package="msPurity")
#'     }else{
#'       model_out <- createBiologicalSVM(librarySdfPth)
#'       model <- model_out[[1]]
#'       modelDistances <- model_out[[1]]
#'     }
#'     # get population metrics
#'     popSd <- sd(modelDistances)*sqrt((length(modelDistances)-1)/(length(modelDistances)))
#'     popMean <- mean(modelDistances)
#'   }
#'
#'   librarySdfPth <- '/media/sf_DATA/structures_converted_small.sdf'
#'   librarySdfPth <- '/media/sf_DATA/structures_converted.sdf'
#'   querySdfPth <- '~/testhmdb.sdf'
#'
#'   # Get the query kernal (this also gives all comparisons to the library)
#'   # probably easier to split the query into chunks
#'   # need to remove snynomns as there is something which breaks code. Also need to remove very small
#'   # compounds (e.g. where there is one element) but potentially we just remove anything less than 50 Da as that
#'   # is the limit of our instrumentation
#'
#'   K1 <- Rchemcpp::sd2gram(querySdfPth,librarySdfPth)
#'
#'   # check with model
#'   queryKernal <- ::as.kernelMatrix(K1[,SVindex(model), drop=FALSE])
#'   queryDistance <- predict(model,queryKernal,type='decision')
#'
#'   # get pvalue for model
#'   pval <- pnorm(queryDistance, popMean, popSd)
#'
#'   maxSim <- max(K1)
#'   K1nm <- colnames(K1)
#'   K1nm <- K1nm[order(K1, decreasing = TRUE)]
#'   K1 <- K1[order(K1, decreasing = TRUE)]
#'
#'   numhits <- length(K1[K1>klim])
#'   tophits <- K1nm[1:topn]
#'   topscores <- K1[1:topn]
#'
#'   return(c('distance'=qeuryDistance, 'pval'=pval,maxSim=maxSim, 'numhits'= numhits,
#'            'tophits'= tophits,  'topscores' =topscores))
#'
#'
#'
#' }
#'
#' createBiologicalSVM <- function(librarySdfPth){
#'   #https://stackoverflow.com/questions/1753299/help-using-predict-for-kernlabs-svm-in-r
#'
#'   KCA <- sd2gramSpectrum(librarySdfPth,  detectArom=FALSE,depthMax=4,silentMode=TRUE)
#'   model <- ksvm(KCA, kernel="matrix",type="one-svc")
#'
#'   # get the distance between the kernal (will be negative if outside and positive inside)
#'   modelDistances <- predict(model,as.kernelMatrix(KCA[,SVindex(model),drop=FALSE]),type='decision')
#'
#'   return(list(model, modelDistances))
#' }
