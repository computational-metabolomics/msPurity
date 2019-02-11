#' #' @title Combine Annotations
#' #'
#' #' @description
#' #'
#' #' @examples
#' #'
#' #' @export
#' combineAnnotations <- function(sqlite_pth, metfrag_result, sirius_csi_result, probmetab_result){
#'   ########################################################
#'   # Check sqlite database has spectral matching results
#'   ########################################################
#'
#'   if (!is.data.frame(grp_peaklist)){
#'     if (is.null(xsa)){
#'       grp_peaklist <- xcms::peakTable(xset)
#'     }else{
#'       grp_peaklist <- CAMERA::getPeaklist(xsa)
#'     }
#'
#'     grp_peaklist <- data.frame(cbind('grpid'=1:nrow(grp_peaklist), grp_peaklist))
#'   }
#'
#'   message("Creating a database of fragmentation spectra and LC features")
#'   target_db_pth <- export_2_sqlite(pa, grp_peaklist, xset, xsa, out_dir, db_name)
#'
#'   return(target_db_pth)
#'
#' }
#'
#'
