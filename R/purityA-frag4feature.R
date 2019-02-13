#' @title Assign precursor purity scored fragmentation spectra to XCMS features
#'
#' @description
#'
#' Assign fragmentation spectra (MS/MS) scored via msPurity package to features
#' from an XCMS set object.
#'
#' Allows the user to filter out spectra below a certain threshold for purity.
#'
#' @aliases frag4feature
#'
#' @param pa object; purityA object
#' @param xset object; XCMS object derived from the same files as the puritydf
#' @param CSVfile character; file to be able to match MS files with their MSMS files
#' @param ppm numeric; ppm tolerance between precursor mz and feature mz
#' @param plim numeric; min purity of precursor to be included
#' @param intense boolean; If the most intense precursor or the centered precursor is used
#' @param use_group boolean; Ignore individual peaks and just find matching fragmentation spectra within the (full) rtmin rtmax of each grouped feature
#' @param convert2RawRT boolean; If retention time correction has been used in XCMS set this to TRUE
#' @param create_db boolean; SQLite database will be created of the results
#' @param db_name character; If create_db is TRUE, a custom database name can be used, default is a time stamp
#' @param out_dir character; Path where database will be created
#' @param grp_peaklist dataframe [optional]; Can use any peak dataframe to add to databse. Still needs to be derived from the xset object though
#' @return purityA object with slots for fragmentation-XCMS links
#'
#' @examples
#'
#' msmsPths <- list.files(system.file("extdata", "lcms", "mzML", package="msPurityData"), full.names = TRUE, pattern = "MSMS")
#' xset <- xcms::xcmsSet(msmsPths, nSlaves = 1)
#' xset <- xcms::group(xset)
#' xset <- xcms::retcor(xset)
#' xset <- xcms::group(xset)
#'
#' pa  <- purityA(msmsPths, interpol = "linear")
#' pa <- frag4feature(pa, xset, CSVfile)
#'
#' @export
setMethod(f="frag4feature", signature="purityA",
          definition = function(pa, xset, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE, create_db=FALSE,
                                out_dir='.', db_name=NA, grp_peaklist=NA){

  #Verify if there is an assess-purity input
  if(is.null(pa)){
    message("no pa files")
    return(NULL)

  }
  #Verify that there is a xset group input
  if(is.null(xdata)){
    message("no xdata files")
    return(NULL)
  }
  #Verify that there is a CSV file input to match files each one with each other
  if(is.null(CSVfile)){
    message("no CSV file")
    return(NULL)
  }else{
    #Have the data frame of the CSVfile
    fileMatch <- read.csv2(file=CSVfile, header=FALSE)
    names(fileMatch) <- c("MS1","MS2")
  }

  cat("\n===================================================================================================\n")
  cat("Processing",length(fileNames(xdata)),"xdata file(s) and",length(pa@fileList),"pa file(s)...\n")

  #Run process one couple of files by one couple
  grped_ms2 <- NULL
  grpall <- NULL
  use_group = FALSE
  for(i in 1:nrow(fileMatch)) {

    fileMS1 <- as.character(fileMatch[i,1])
    fileMS2 <- as.character(fileMatch[i,2])
    cat("\n====== Running",fileMS1,"with",fileMS2,"=====\n")

    #Matching the good MS/MS file
    numberofMS2file <- 0
    for(j in 1:length(pa@fileList)){
      if(names(pa@fileList[j]) == fileMS2){
        pathfileMS2 <- unname(pa@fileList[j])
        numberofMS2file <- j
        return
      }
    }

    #Matching the good MS file
    numberofMS1file <- 0
    for(l in 1:nrow(xdata@phenoData@data)){
      if(rownames(xdata@phenoData@data[l,]) == paste0("./",fileMS1)){
        numberofMS1file <- l
      }
    }

    #Verify that we have each file
    if(numberofMS1file == 0 || numberofMS2file == 0){
      error_message <- paste("/!\\ We can't find",fileMS1,"or",fileMS2,"/!\\")
      print(error_message)
      next
    }

    #Pass xdata into xset object
    if(class(xdata) == "XCMSnExp"){
      xset <- getxcmsSetObject(xdata)
    }else{
      xset <- xdata
    }

    #Verify if xset and pa comes from the same file or not
    if(fileMS1 != fileMS2){
      #Carreful if use_group = FALSE and if it is not the first couple of file !!
      if(i > 1 && use_group == FALSE) {
        error_message <- "/!\\ You already match with peaks and want now match group /!\\ \n\tEx not to do : LCMSMS1;LCMSMS1 first then now LCMS1;LCMSMS1"
        print(error_message)
        stop(error_message)
      } else {
        use_group = TRUE
        pa@f4f_link_type = 'group'
      }
    }else{
      #Carreful if use_group = TRUE already !!
      if(use_group == TRUE){
        error_message <- "/!\\ You already match with group peaks and want now match only on peaks /!\\ \n\tEx not to do : LCMS1;LCMSMS1 first and then now LCMSMS1;LCMSMS1 "
        print(error_message)
        stop(error_message)
      }
      use_group = FALSE
      pa@f4f_link_type = 'individual'
    }  

    # Get the purity data frame and the xcms peaks data frame of the good file
    puritydf <- pa@puritydf[which(pa@puritydf[,"fileid"] == numberofMS2file),]
    puritydf$fileid <- as.numeric(puritydf$fileid)
    cat(paste("Stock",nrow(puritydf),"MS/MS spectra from assess-purity\n"))

    # Check if is going to be multi-core
    if(pa@cores>1){
      cl <- parallel::makeCluster(pa@cores)
      doSNOW::registerDoSNOW(cl)
      para = TRUE
    }else{
      para = FALSE
    # perform multicore
    }
    para = FALSE

    if(use_group) {
      #Select groups which contain peaks in the same class as file class
      for(x in 1:nrow(xset@phenoData)) {
        if(basename(rownames(xset@phenoData)[x]) == fileMS1) {
          MSclass <- as.character(xset@phenoData[x,"class"])
        }
      }
      fullpeakw <- xset@groups[which(xset@groups[,MSclass] > 0),]
      fullpeakw <- data.frame(get_full_peak_width(fullpeakw, xset))
      fullpeakw$grpid <- seq(1, nrow(fullpeakw))
      cat("Stock",nrow(fullpeakw),"group peaks from class",MSclass,"from xset\n")

      # Map xcms features to the data frame (takes a while)
      matched <- plyr::ddply(puritydf, ~ pid, fsub2, allpeaks = fullpeakw, intense = intense, ppm = ppm, fullp = TRUE, use_grped = TRUE)

      if(nrow(matched) > 0) {
        grpedp <- matched
        savegrpid <- grpedp[,"grpid"]
        coldontwant <- "grpid"
        grpedp <- grpedp[, ! names(grpedp) %in% coldontwant, drop = F]
        grpedp <- cbind(savegrpid,grpedp)
        names(grpedp)[1] <- "grpid"
        cat(nrow(grpedp))
        cat(" peaks matched with parents ions\n")
      } else {
        error_message <- "/!\\ No peaks matched for these files /!\\\n"
        cat(error_message)
        next
      }
      
    } else {
      allpeaks <- data.frame(xset@peaks[which(xset@peaks[,"sample"] == numberofMS1file),])
      allpeaks$cid <- seq(1, nrow(allpeaks))
      allpeaks <- plyr::ddply(allpeaks, ~ sample, getname, xset = xset)

      if(convert2RawRT) {
        allpeaks$rtminCorrected <- allpeaks$rtmin
        allpeaks$rtmaxCorrected <- allpeaks$rtmax
        allpeaks <- plyr::ddply(allpeaks, ~ sample, convert2Raw, xset = xset)
      }

      cat(paste("Stock",nrow(allpeaks),"peaks from xset preprocessing\n"))
      # Map xcms features to the data frame (takes a while)
      matched <- plyr::ddply(puritydf, ~ fileid, .parallel = para, fsub1, allpeaks=allpeaks, ppm = ppm, intense = intense)

      #Group by the xcms groups
      if(nrow(matched) > 0) {
        grpedp <- plyr::llply(xset@groupidx, grpByXCMS, matched = matched)
        names(grpedp) <- seq(1, length(grpedp))
        grpedp <- plyr::ldply(grpedp, .id = TRUE)
        colnames(grpedp)[1] <- "grpid"
        cat(nrow(grpedp))
        cat(" peaks matched with parents ions\n")
      } else {
        error_message <- "/!\\ No peaks matched for these files /!\\\n"
        cat(error_message)
        next
      }    
    }

    # Filter out any precursor below purity threshold
    if (!is.na(plim) && plim>0){
      grpm <- grpm[grpm$inPurity>plim,]
    }

    #Output if no peaks matched
    if(nrow(grpedp) == 0) {
      message("/!\\ No peaks matched /!\\")
      next
    }

    # Add some extra info for filtering purposes
    #grpm <- merge(grpedp, shrt, by = c('pid', 'fileid', 'seqNum'))
    grpm <- grpedp
    # Make sure order is by grpid
    grpm <- grpm[order(grpm$grpid),]
    # Filter out any precursor below purity threshold
    if (!is.na(plim) && plim > 0) {
      grpm <- grpm[grpm$inPurity > plim,]
    }
    #Merge informations from each file
    grpall <- rbind(grpall,grpm)

    if(!is.list(grped_ms2)) {
      #For the first file processed
      grped_ms2 <- getMS2scans(grpm, pathfileMS2, mzRback = pa@mzRback)
    } else {
      #For the others files
      grped_ms2 <- c(grped_ms2,getMS2scans(grpm, pathfileMS2, mzRback = pa@mzRback))
    }

    if (create_db) {
      pa@db_path <- create_database(pa = pa, xset = xset, out_dir = out_dir, db_name = db_name, grp_peaklist = grp_peaklist, grped_df = grpall, fileMS1 = numberofMS1file, fileMS2 = numberofMS2file)
    }
  }

  #Save data in pa object
  if(is.null(grpall)){
    error_message <- "/!\\ Nothing found for grped_df /!\\"
    print(error_message)
    stop(error_message)
  } else {
    pa@grped_df <- grpall
  }
  if(is.null(grped_ms2)){
    error_message <- "/!\\ Nothing found for grped_ms2 /!\\"
    print(error_message)
    stop(error_message)
  } else {
    pa@grped_ms2 <- grped_ms2
  }

  cat("===================================================================================================\n")
  return(pa)
})

fsub1  <- function(prod, allpeaks, intense, ppm) {
  # go through all the MS/MS files from the each file
  allpeakfile <- allpeaks[allpeaks$filename==unique(prod$filename),]

  grpdFile <- plyr::ddply(prod, ~ seqNum,
                          fsub2, # FUNCTION
                          allpeaks = allpeakfile,
                          intense = intense,
                          ppm = ppm)
}

fsub2  <- function(pro, allpeaks, intense, ppm, fullp = FALSE, use_grped=FALSE) {

  # check for each MS/MS scan if there is an associated feature found in that region for that file
  if(intense) {
    mz1 <- pro$iMz
  } else {
    if (is.na(pro$aMz)) {
      mz1 <- pro$precursorMZ
    } else {
      mz1 <- pro$aMz
    }
  }

  if(is.na(mz1) | is.null(mz1)) {
    return(NULL)
  }

  prt <- pro$precursorRT

  if(is.na(prt)) {
    prt <- pro$retentionTime
  }

  if(fullp) {
    mtchRT <- allpeaks[prt >= allpeaks$rtmin_full & prt <= allpeaks$rtmax_full, ]
  } else {
    mtchRT <- allpeaks[prt >= allpeaks$rtmin & prt <= allpeaks$rtmax, ]
  }

  if(nrow(mtchRT)==0) {
    return(NULL)
  }

  if (use_grped){
    # can only use fullp when using the grouped peaklist
    mtchMZ <- plyr::ddply(mtchRT, ~ grpid, mzmatching, mz1=mz1, ppm=ppm, pro=pro)
  }else{
    mtchMZ <- plyr::ddply(mtchRT, ~ cid, mzmatching, mz1=mz1, ppm=ppm, pro=pro)
  }

  return(mtchMZ)
}



check_ppm <- function(mz1, mz2) { return(abs(1e6*(mz1-mz2)/mz2)) }

getMS2scans  <- function(grpm, filepths, mzRback) {

  # Get all MS2 scans
  scans <- getscans(filepths, mzRback)

  if(length(filepths)==1){
    scans = list(scans)
  }

  grpm$fid <- seq(1, nrow(grpm))

  ms2l <- plyr::dlply(grpm, ~ grpid, getScanLoop, scans=scans)

  return(ms2l)
}


mzmatching <- function(mtchRow, mz1=mz1, ppm=ppm, pro=pro){

  if ('mzmed' %in% colnames(mtchRow)){
    mz2 <- mtchRow$mzmed
  }else{
    mz2 <- mtchRow$mz
  }

  ppmerror <- check_ppm(mz1, mz2)

  if(ppmerror<ppm){

    mtchRow$inPurity <- pro$inPurity
    mtchRow$pid <- pro$pid
    mtchRow$precurMtchID <- pro$seqNum
    mtchRow$precurMtchScan <- pro$precursorScanNum
    mtchRow$precurMtchRT <- pro$precursorRT
    mtchRow$precurMtchMZ <- mz1
    mtchRow$precurMtchPPM <- ppmerror
    mtchRow$retentionTime <- pro$retentionTime
    mtchRow$fileid <- pro$fileid
    mtchRow$seqNum <- pro$seqNum
    return(mtchRow)
  }else{
    return(NULL)
  }
}

getScanLoop <- function(peaks, scans){
  grpl <-  list()

  if ('sample' %in% colnames(peaks)){
    idx_nm ='sample'
  }else{
    idx_nm = 'fileid'
  }

  for(i in 1:nrow(peaks)){
    x <- peaks[i,]
    grpl[[i]] <- scans[[1]][[x$precurMtchID]]
  }
  return(grpl)
}

getname <- function(x, xset){
  x$filename <- basename(xset@filepaths[x$sample])
  return(x)
}

grpByXCMS <- function(x, matched){
  matched[matched$cid %in% x,]
}

convert2Raw <- function(x, xset){
  sid <- unique(x$sample)
  # for each file get list of peaks
  x$rtmin <- xset@rt$raw[[sid]][match(x$rtmin, xset@rt$corrected[[sid]])]
  x$rtmax <- xset@rt$raw[[sid]][match(x$rtmax, xset@rt$corrected[[sid]])]
  return(x)
}

# This function retrieve a xset like object
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getxcmsSetObject <- function(xobject) {
    # XCMS 1.x
    if (class(xobject) == "xcmsSet")
        return (xobject)
    # XCMS 3.x
    if (class(xobject) == "XCMSnExp") {
        # Get the legacy xcmsSet object
        suppressWarnings(xset <- as(xobject, 'xcmsSet'))
        if (!is.null(xset@phenoData$sample_group))
            sampclass(xset) <- xset@phenoData$sample_group
        else
            sampclass(xset) <- "."
        return (xset)
    }
}