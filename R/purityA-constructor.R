NULL

#' @title Assess the purity of multiple LC-MS/MS or DI-MS/MS files (constructor)
#'
#' @description
#' Constructor for the purityA class.
#'
#' Given a vector of LC-MS/MS or DI-MS/MS mzML file paths calculate the
#' precursor purity of each MS/MS scan
#'
#' Will automatically determine the isolation widths offsets from the mzML file.
#' For some vendors though this is not recorded (Agilent).
#' In these cases the offsets should be given as a parameter.
#'
#' In the case of Agilent only the "narrow" isolation is supported.
#' This roughly equates to +/- 0.65 Da (depending on the instrument). If the
#' file is detected as originating from an Agilent instrument the isolation
#' widths will automatically be set as +/- 0.65 Da.
#'
#' @param fileList vector; mzML file paths for MS/MS spectra
#' @param cores numeric; Number of cores to use
#' @param mostIntense boolean; True if the most intense peak is used for calculation. False if the centered peak is used
#' @param nearest boolean; True if the peak selected is from either the preceding scan or the nearest.
#' @param offsets vector; Overide the isolation offsets found in the mzML filee.g. c(0.5, 0.5)
#' @param plotP boolean; If TRUE a plot of the purity is to be saved
#' @param plotdir vector; If plotP is TRUE plots will be saved to this directory
#' @param interpol character; type of interolation to be performed "linear" or "spline" (Spline option is only included for testing purposes,
#'                            linear should be used for all standard cases, isotope removal is also not available for the spline option)
#' @param iwNorm boolean; If TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5\% (0.05)
#' @param isotopes boolean; TRUE if isotopes are to be removed
#' @param im matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#' @param mzRback character; backend to use for mzR parsing
#'
#' @return a dataframe of the purity score of the ms/ms spectra
#'
#' @examples
#' filepths <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")
#' pa <- purityA(filepths)
#' @seealso \code{\link{assessPuritySingle}}
#' @export
purityA <- function(filepathsMS2, 
                    filepathsMS1 = NULL, 
                    CSVfile = NULL,
                    forcedMS1 = FALSE,
                    cores=1,
                    mostIntense=FALSE,
                    nearest=TRUE,
                    offsets=NA,
                    plotP=FALSE,
                    plotdir=NULL,
                    interpol="linear",
                    iwNorm=FALSE,
                    iwNormFun=NULL,
                    ilim=0.05,
                    mzRback='pwiz',
                    isotopes=TRUE,
                    im=NULL){

  if((length(filepathsMS2)>=1) && (filepathsMS2 == "" )){
    message("no MS2 file list")
    return(NULL)
  }
  
  requireNamespace('foreach')

  #Build the data frame of file matching
  #First case is when you only have filepathsMS2 (when MS and MSMS are in the same file)
  #Second case is when you have filepathsMS1 and CSVfile
  #Other cases are false and return NULL cause when you have filepathsMS1 and filepathsMS2 you always need to match them with a CSVfile
  if(is.null(filepathsMS1) && (is.null(CSVfile))){
    print("No MS1 and no CSV")
    filepathsMS1 <- filepathsMS2
    fileMatch <- data.frame(MS1 = filepathsMS1, MS2 = filepathsMS2)
    rownames(fileMatch) <- 1:nrow(fileMatch)
  }else{
    if(!is.null(filepathsMS1)){
      if(!is.null(CSVfile)){
        #Have the data frame of the CSVfile
        fileMatch <- read.csv2(file=CSVfile, header=FALSE)
        names(fileMatch) <- c("MS1","MS2")
        print(fileMatch)
      }else{
        message("no CSV file")
        return(NULL)
      }
    }else{
      message("You have no filepaths for MS1 files")
      return(NULL)
    }   
  }
  #Verify if the user set the extension in the fileMatch
  print("TODO : verify if the user set the extension in the fileMatch")

  #Verify that we have all the files containing in fileMatch and reorder them if needed
  if(NA %in% filepathsMS1[match(fileMatch[,1],filepathsMS1)]){
    cat("It misses some MS1 files to be able to match them !\n")
    #quit(save="no",status="It misses some MS1 files to be able to match them !")
  }else{
    filepathsMS1 <-filepathsMS1[match(fileMatch[,1],filepathsMS1)] 
  }
  if(NA %in% filepathsMS2[match(fileMatch[,2],filepathsMS2)]){
    cat("It misses some MS2 files to be able to match them !\n")
    #quit(save="no",status="It misses some MS1 files to be able to match them !")
  }else{
    filepathsMS2 <-filepathsMS2[match(fileMatch[,2],filepathsMS2)] 
  }
  print(filepathsMS2)
  #Build the purityA object
  pa <- new("purityA", fileList = filepathsMS2, fileListMS1 = filepathsMS1, fileMatch = fileMatch, cores = cores, mzRback=mzRback)

  # Check cores and choose if parallel or not (do or dopar)
  if(pa@cores<=1){
    operator <- foreach::'%do%'
  } else{
    cl1 <- parallel::makeCluster(pa@cores, type="SOCK")
    doSNOW::registerDoSNOW(cl1)
    operator <- foreach::'%dopar%'
  }

  #Run MS2 file(s) in parallel (or not) using foreach
  purityL <- operator(foreach::foreach(i = 1:nrow(fileMatch),
                                  .packages = 'mzR'),
                                  assessPuritySingle(filepathMS2 = pa@fileList[fileMatch[i,2]], 
                                                     filepathMS1 = filepathsMS1[fileMatch[i,1]], 
                                                     CSVfile = CSVfile,
                                                     forcedMS1 = forcedMS1,
                                                     mostIntense = mostIntense,
                                                     nearest=nearest,
                                                     offsets = offsets,
                                                     plotP = plotP,
                                                     plotdir = plotdir,
                                                     interpol = interpol,
                                                     iwNorm = iwNorm,
                                                     iwNormFun = iwNormFun,
                                                     ilim = ilim,
                                                     mzRback = mzRback,
                                                     isotopes = isotopes,
                                                     im = im
                                                     ))

  if(pa@cores>1){
      parallel::stopCluster(cl1)
  }

  names(purityL) <- seq(1, length(purityL))

  puritydf <- plyr::ldply(purityL, .id = TRUE)
  
  if(nrow(puritydf)>0){
    colnames(puritydf)[1] <- "fileid"
    pid <- seq(1, nrow(puritydf))
    puritydf <- cbind(pid, puritydf)
  }

  #Remove duplicated scans (maybe add columns filename and precursorFilename??????????????????????????????????????????????????)
  puritydf <- puritydf[!duplicated(puritydf[,-c(1,2)]),]

  pa@puritydf <- puritydf
  pa@cores = cores

  return(pa)
}


#' @title Assess the purity of a single LC-MS/MS or DI-MS/MS file
#'
#' @description
#'
#' Given a filepath to an mzML file the precursor purity for any MS/MS scans
#' will be outputed into a dataframe
#'
#' @param filepth character; mzML file path for MS/MS spectra
#' @param fileid numeric; adds a fileid column (primarily for internal use for msPurity)
#' @param mostIntense boolean; True if the most intense peak is used for calculation. False if the centered peak is used
#' @param nearest boolean; True if the peak selected is as the nearest MS1 scan. If False then the preceding scan is used
#' @param offsets vector; Overide the isolation offsets found in the mzML filee.g. c(0.5, 0.5)
#' @param cores numeric; Number of cores to use
#' @param plotP boolean; If TRUE a plot of the purity is to be saved
#' @param plotdir vector; If plotP is TRUE plots will be saved to this directory
#' @param interpol character; Type of interolation to be performed "linear", "spline" or "none"
#' @param iwNorm boolean; If TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function; A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric; All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5\% (0.05)
#' @param mzRback character; Backend to use for mzR parsing
#' @param isotopes boolean; TRUE if isotopes are to be removed
#' @param im matrix; Isotope matrix, default removes C13 isotopes (single, double and triple bonds)
#' @return a dataframe of the purity score of the ms/ms spectra
#'
#' @examples
#' filepth <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")
#'
#' puritydf <- assessPuritySingle(filepth)
#' @seealso \code{\link{purityA}}
#' @export
assessPuritySingle <- function(filepathMS2, 
                               filepathMS1 = NULL, 
                               CSVfile = NULL,
                               forcedMS1 = FALSE,
                               fileid=NA,
                               mostIntense=FALSE,
                               nearest=TRUE,
                               offsets=NA,
                               cores=1,
                               plotP=FALSE,
                               plotdir=NULL,
                               interpol="linear",
                               iwNorm=FALSE,
                               iwNormFun=NULL,
                               ilim=0,
                               mzRback='pwiz',
                               isotopes=TRUE,
                               im=NULL){

  cat("\n\n===================== ")
  cat(basename(filepathMS2))
  cat(" =====================\n")

  #=================================
  # Load in files and initial setup
  #=================================
  # Get the mzR dataframes
  mrdf <- getmrdf(filepathMS2, filepathsMS1 = filepathMS1, CSVfile = CSVfile, mzRback)
  if(is.null(mrdf)){
    message(paste("/!\\STOP/!\\ No MS/MS spectra for file: ", basename(filepathMS2)))
    return(NULL)
  }else if(length(mrdf) < 2){
    if(mrdf == FALSE){
      message(" /!\\STOP/!\\ You must have a MS1 fileList cause your file has only MS2")
      return(NULL)
    }
  }

  # having id column makes things easier to track
  mrdf$id <- seq(1, nrow(mrdf))

  # add the fileid info
  if(!is.na(fileid)){
    mrdf$fileid <- rep(fileid, nrow(mrdf))
  }

  # add filename (already done no?)
  mrdf$filename <- basename(filepathMS2)

  if(!("precursorFilename" %in% names(mrdf))){
    mrdf$precursorFilename <- mrdf$filename
  }
  # Get a shortened mzR dataframe
  mrdfshrt <- mrdf[mrdf$msLevel==2,][,c("seqNum","acquisitionNum","precursorIntensity",
                                        "precursorMZ", "precursorRT",
                                        "precursorScanNum", "id", "retentionTime", "filename", "precursorFilename")]

  #if((length(unique(mrdf$msLevel))<2) && (unique(mrdf$msLevel)==2)){
  #  message("only MS2 data, not possible to calculate purity")
  #  mrdfshrt[ , c("precursorNearest", "aMz", "aPurity", "apkNm", "iMz", "iPurity", "ipkNm", "inPkNm", "inPurity")] <- NA
  #  return(mrdfshrt)
  #}

  # get scans of MS2 (list of mz and i) from mzR
  scansMS2 <- getscans(filepathMS2, mzRback)
  # get scans of MS1 (list of mz and i) from mzR
  if(unique(mrdf$precursorFilename) == basename(filepathMS1)){
    print("precursorFilename et names = on prend le names")
    scansMS1 <- getscans(filepathMS1,mzRback)
  }else{
    print("precursorFilename et names != on regarde forcedMS1")
    if(forcedMS1){
      print("L'utilisateur force le MS1")
      scansMS1 <- getscans(filepathMS1,mzRback)
    }else{
      print("MS1 non forcé")
      if(unique(mrdf$precursorFilename) == basename(filepathMS1)){
        print("precursorFilename et unname = on prend unname (mais c bizarre ça ne devrait pas arriver !)")
        scansMS1 <- getscans(filepathMS1,mzRback)
      }else{
        print("Tout est différent du précurseur, faut faire quelquechose ! Genre chercher le file du prec ou bien faire quitter ou bien regarder si = filepathMS2")
        if(unique(mrdf$filename) == unique(mrdf$precursorFilename)){
          print("Les filename et precursorFilename sont les memes")
          scansMS1 <- getscans(filepathMS1,mzRback)
        }
      }
    }
  }
  #TODO verify if this function works well
  ########
  # Get offsets from mzML unless defined by user
  if(anyNA(offsets)){
    offsets <- get_isolation_offsets(filepathMS2)
  }
  minoff <- offsets[1]
  maxoff <- offsets[2]
  ########

  if(is.null(iwNormFun)){
    iwNormFun <- iwNormGauss(minOff = -minoff, maxOff = maxoff)
  }

  # For n MS1 scans before and after the MS2 scan. For linear interpolation
  # only two points needed. More needed for spline
  if((interpol=="linear") || (interpol=="none")){
    # if linear interpolation we only need the purity before and after.
    # If no interpolation required we still need to find the nearest ms1 scan
    # so the nump is 1 as well
    nump <- 1
  }else if (interpol=="spline"){
    # if performing spline we need multiple data points
    nump <- 10
  }
  # Get the precursor (ms1) scans
  prec_scans <- get_prec_scans(mrdf, nump)

  # add a column with the nearest precursor for all ms2 scans
  mrdfshrt$precursorNearest <- plyr::laply(prec_scans, function(x){ return(x$nearest)})

  #=====================================
  # Initial purity calculation
  #=====================================
  # Get the initial purity score and target mz associated with every ms2 scan.
  # Calculate the purity on the nearest ms1 scan. This function also determines
  # which is the target peak used for interpolations
  if(unique(mrdfshrt$filename) == unique(basename(mrdfshrt$precursorFilename))){
    print("We have the same filename and precursorFilename")
    initial <- plyr::ddply(mrdfshrt, ~ seqNum, get_init_purity,
                           scansMS2=scansMS2,
                           scansMS1=NULL,
                           minoff=minoff,
                           maxoff=maxoff,
                           nearest=nearest,
                           mostIntense=mostIntense,
                           iwNorm=iwNorm,
                           iwNormFun=iwNormFun,
                           ilim=ilim,
                           isotopes=isotopes,
                           im=im)
  }else{
    print("We have different filename and precursorFilename")
    initial <- plyr::ddply(mrdfshrt, ~ seqNum, get_init_purity,
                           scansMS2=scansMS2,
                           scansMS1=scansMS1,
                           minoff=minoff,
                           maxoff=maxoff,
                           nearest=nearest,
                           mostIntense=mostIntense,
                           iwNorm=iwNorm,
                           iwNormFun=iwNormFun,
                           ilim=ilim,
                           isotopes=isotopes,
                           im=im)
  }

  # Add results to mzR dataframe
  mrdfshrt <- cbind(mrdfshrt, initial)
  print(paste(nrow(mrdfshrt),"scans added to the MS2 scans"))

  # If no interpolation required then all purity calculations are completed
  if(interpol=="none"){
    return(mrdfshrt)
  }

  #=====================================
  # Interpolation
  #=====================================
  # The interpolation can be time consuming if using spline interpolation,
  # so can be run in parallel
  if(cores>1){
    requireNamespace(parallel)
    clust<-parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(clust)
    pBool = TRUE
  }else{
    pBool = FALSE
  }
  interDF <- plyr::ddply(mrdfshrt, ~ seqNum, .parallel = pBool,
                                   get_interp_purity, # FUNCTION
                                   scan_peaks=c(scansMS1,scansMS2),
                                   ppm=7, # hard coded at the moment
                                   ms2=mrdf[mrdf$msLevel==2,]$seqNum,
                                   prec_scans=prec_scans,
                                   minoff=minoff,
                                   maxoff=maxoff,
                                   mostIntense=mostIntense,
                                   plotP=plotP,
                                   plotdir=plotdir,
                                   interpol=interpol,
                                   nearest=nearest,
                                   iwNorm=iwNorm,
                                   iwNormFun=iwNormFun,
                                   ilim=ilim,
                                   isotopes=isotopes,
                                   im=im)
  if(cores>1){
    requireNamespace(parallel)
    parallel::stopCluster(clust)
  }

  if(interpol=="spline"){
    mrdfshrt$inPurity <- interDF$inPurity
  }else{
    mrdfshrt$inPkNm <- interDF$inPkNm
    mrdfshrt$inPurity <- interDF$inPurity
  }
  return(mrdfshrt)
}


get_init_purity <- function(ms2h, scansMS2, scansMS1, minoff, maxoff, nearest,
                            mostIntense, iwNorm, iwNormFun, ilim, isotopes, im){

  #========================================
  # Get the scan to perform calculations on
  #========================================
  if (nearest){
    precScn <- ms2h$precursorNearest
  }else{
    precScn <- ms2h$precursorScanNum
  }

  # What the mzML says the mz and intensity of the target is
  precMZ <- ms2h$precursorMZ
  precI <- ms2h$precursorIntensity

  mzmin <-  precMZ - minoff
  mzmax <-  precMZ + maxoff

  if(is.null(scansMS1)){

    allscans <- scansMS2[[precScn]]

    #=======================================
    # Get centered target peak
    #=======================================
    # Get target based on matching to be within isolation window of precursorMZ
    # and 1 decimal of intensity
    target_approx <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax)
                            & (round(allscans[,2],1) == round(precI,1)),]

    # no precusor within isolation window and intensity
    if(length(target_approx)==0){
      # get nearest match for just mz
      target_approx <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax),]
      if(length(target_approx)==0){
        return(c("aMz"= NA, "aPurity" = 0, "apkNm" = NA,
                 "iMz" = NA, "iPurity" = 0, "ipkNm" = NA ))
      }
    }

    if (length(target_approx)==2) {
      target <- target_approx
    }else{
      # Get the mz that is closest to the target
      target <- target_approx[which.min(abs(target_approx[,1]-precMZ)),]
    }

    # The centered peak
    aMz <- target[1]

    #==================================================
    # Get most intense peak within isolation window
    #==================================================
    subp <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax),]

    if(is.vector(subp)){
      iTarget <- subp
    }else if (nrow(subp)==0){
      # nothing within in the mz tolerance range
      return(c("aMz"=aMz, "aPurity" = 0, "iMz" = NA, "iPurity" = 0 ))
    }else{
      iTarget <- subp[order(subp[,2], decreasing = TRUE),][1,]
    }

    # the most intense peak
    iMz <- iTarget[1]

    #==================================================
    # Calculate initial purity
    #==================================================
    if(is.vector(subp)){
      aPurity = iPurity = apkNm = ipkNm =  1
    } else if (nrow(subp)==1){
      aPurity = iPurity = apkNm = ipkNm =  1
    } else {
      pouta <- pcalc(peaks=subp, mzmin=mzmin, mzmax=mzmax, mztarget=aMz, ppm=NA,
                     iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim, isotopes=isotopes, im=im)
      aPurity <- unname(pouta[1])
      apkNm <- unname(pouta[2])
      pouti <- pcalc(peaks=subp, mzmin=mzmin, mzmax=mzmax, mztarget=iMz, ppm=NA,
                     iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim, isotopes=isotopes, im=im)
      iPurity <- unname(pouti[1])
      ipkNm <- unname(pouti[2])

    }
    fileinfo <- c("aMz"=aMz, "aPurity" = aPurity, "apkNm" = apkNm,
                "iMz" = iMz, "iPurity" = iPurity, "ipkNm" = ipkNm )
  }else{
   
    scans<-c(scansMS1,scansMS2)
    allscans <- scans[[precScn]]
    
    #=======================================
    # Get centered target peak
    #=======================================
    # Get target based on matching to be within isolation window of precursorMZ
    # and 1 decimal of intensity
    target_approx <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax)
                            & (round(allscans[,2],1) == round(precI,1)),]
    # no precusor within isolation window and intensity
    if(length(target_approx)==0){
      # get nearest match for just mz
      target_approx <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax),]
      if(length(target_approx)==0){
        return(c("aMz" = NA, "aPurity" = 0, "apkNm" = NA,
                 "iMz" = NA, "iPurity" = 0, "ipkNm" = NA ))
      }
    }
    if (length(target_approx)==2) {
      target <- target_approx
    }else{
      # Get the mz that is closest to the target
      target <- target_approx[which.min(abs(target_approx[,1]-precMZ)),]
    }
    # The centered peak
    aMz <- target[1]
    #==================================================
    # Get most intense peak within isolation window
    #==================================================
    subp <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax),]
    if(is.vector(subp)){
      iTarget <- subp
    }else if (nrow(subp)==0){
      # nothing within in the mz tolerance range
      return(c("aMz"=aMz, "aPurity" = 0, "iMz" = NA, "iPurity" = 0 ))
    }else{
      iTarget <- subp[order(subp[,2], decreasing = TRUE),][1,]
    }
    # the most intense peak
    iMz <- iTarget[1]
    #==================================================
    # Calculate initial purity
    #==================================================
    if(is.vector(subp)){
      aPurity = iPurity = apkNm = ipkNm =  1
    } else if (nrow(subp)==1){
      aPurity = iPurity = apkNm = ipkNm =  1
    } else {
      pouta <- pcalc(peaks=subp, mzmin=mzmin, mzmax=mzmax, mztarget=aMz, ppm=NA,
                     iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim, isotopes=isotopes, im=im)
      aPurity <- unname(pouta[1])
      apkNm <- unname(pouta[2])
      pouti <- pcalc(peaks=subp, mzmin=mzmin, mzmax=mzmax, mztarget=iMz, ppm=NA,
                     iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim, isotopes=isotopes, im=im)
      iPurity <- unname(pouti[1])
      ipkNm <- unname(pouti[2])
    }
    fileinfo <- c("aMz"=aMz, "aPurity" = aPurity, "apkNm" = apkNm,
                "iMz" = iMz, "iPurity" = iPurity, "ipkNm" = ipkNm )
  }
  return(fileinfo)
}

get_interp_purity <- function(rowi, scan_peaks, prec_scans, ms2, ppm,
                              minoff, maxoff, mostIntense=FALSE, plotP=FALSE,
                              plotdir=NULL, interpol, nearest=TRUE,
                              iwNorm, iwNormFun, ilim, isotopes, im){
  
  # get the region of interest for scans
  scanids <- prec_scans[which(ms2==rowi$seqNum)][[1]]
  roi_scns <- scan_peaks[c(scanids$pre, scanids$post)]
  if(interpol=="linear"){
    purity <- linearPurity(rowi, scan_peaks, minoff, maxoff, ppm, scanids,
                           nearest, mostIntense, iwNorm, iwNormFun,
                           ilim, plotP, plotdir, isotopes, im)
  }else if (interpol=="spline"){
    purity <- splinePurity(rowi, roi_scns, minoff, maxoff, ppm,
                           mostIntense, scanids, plotP, plotdir)
  }
  return(purity)
}

linearPurity <- function(rowi, scan_peaks, minoff, maxoff, ppm, scanids,
                         nearest, mostIntense, iwNorm, iwNormFun, ilim,
                         plotP, plotdir, isotopes, im){
  if(nearest){
    scn1 <- rowi$precursorNearest
  }else{
    scn1 <- rowi$precursorScanNum
  }

  if(mostIntense){
    purity1 <- rowi$iPurity
    mztarget1 <- rowi$iMz
    pknm1 <- rowi$ipkNm
  }else{
    purity1 <- rowi$aPurity
    mztarget1 <- rowi$aMz
    pknm1 <- rowi$apkNm
  }

  precMZ <- rowi$precursorMZ

  mzmin <-  precMZ - minoff
  mzmax <-  precMZ + maxoff

  # Calculate ppm tolerance to look either side of the original target value
  targetMaxMZ <- mztarget1 + (mztarget1*0.000001)*ppm
  targetMinMZ <- mztarget1 - (mztarget1*0.000001)*ppm

  # Then calculate purity for the other nearest MS1 scan (purity2)
  if(scanids$pre==scn1){
    scn2 <- scanids$post
  }else{
    scn2 <- scanids$pre
  }

  # The last scan (nothing to interpolate from)
  if((length(scn2)==0) || (is.na(mztarget1))){
    return(c("inPkNm"=pknm1, "inPurity"=purity1))
  }
  # Calculate purity for the other nearest MS1 scan
  allp <- scan_peaks[[scn2]]
  pout <- pcalc(peaks=allp, mzmin=mzmin, mzmax=mzmax, mztarget=mztarget1,
                iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim,
                targetMinMZ = targetMinMZ, targetMaxMZ = targetMaxMZ,
                isotopes = isotopes, im = im)

  purity2 <- pout[1]
  pknm2 <- pout[2]
  ###########################################
  # Interpolate for purity
  ###########################################
  # Cant interpolate if only have 1 NA value
  if(is.na(purity2)){
    inPurity <- purity1
    fpurity <- NULL
  }else if (is.na(purity1)){
    inPurity <- purity2
    fpurity <- NULL
  }else{
    # Create linear functions for poth purity and the number of peaks
    fpurity <- approxfun(c(scn1, scn2), c(purity1, purity2))
    inPurity <- fpurity(scanids$scan)
  }
  ###########################################
  # Interpolate for peak number
  ###########################################
  if(is.na(pknm2)){
    inPkNm <- pknm1
  }else if (is.na(pknm1)){
    inPkNm <- pknm2
  }else{
    # Create linear functions for poth purity and the number of peaks
    fpk <- approxfun(c(scn1, scn2), c(pknm2, pknm2))
    inPkNm <- fpk(scanids$scan)
  }
  # save a plot of the interpolation
  if(plotP){
    df <- data.frame("idx"=c(scn1, scn2),"purity"=c(purity1, purity2))
    intername <- paste(rowi$file, "_", rowi$id, "_interpolate_plot_",
                       mztarget1,".png", sep="")
    name <- file.path(plotdir, intername)
    plotPurity(df, rowi$seqNum, name, fpurity)
  }

  return(c("inPkNm"=inPkNm, "inPurity"=inPurity))
}

get_prec_scans <- function(mrdf, num){
  # Verify if there is more than 1 level in the file
  if(1 %in% unique(mrdf$msLevel)){
    if(2 %in% unique(mrdf$msLevel)){
      print("MS1 and MS2 are in the same file")
      ms1 <- mrdf[mrdf$msLevel==1,]$seqNum
      ms2 <- mrdf[mrdf$msLevel==2,]$seqNum

      #Searching for the pre and post MS scan(s)
      prec_scans <- plyr::alply(ms2, 1, function(scan2){
        post <- ms1[ms1>scan2]
        pre <- ms1[ms1<scan2]

        if(length(pre)>num){
          pre <- pre[(length(pre)-(num-1)):length(pre)]
        }

        if(length(post)>num){
          post <- post[1:num]
        }

        #Searching for the nearest MS scan
        nearest <- ms1[which.min(abs(ms1 - scan2))]

        return(list("pre"=pre, "post"=post, "scan"=scan2, "nearest"= nearest))
      })
    }else{
      print("File with only MS1 data")
      specin <- 1
      next
    }
  }else{
    if(2 %in% unique(mrdf$msLevel)){
      print("MS1 and MS2 are in separated files")
      #Reading MS1 precursor file
      fileToLoad <- unique(mrdf$precursorFilename)
      s_groups <- sapply(fileToLoad, function(x) tail(unlist(strsplit(dirname(x),"/")), n=1))
      s_name <- tools::file_path_sans_ext(basename(fileToLoad))
      pd <- data.frame(sample_name=s_name, sample_group=s_groups, stringsAsFactors=FALSE)
      raw_data <- MSnbase::readMSData(files=fileToLoad, pdata = new("NAnnotatedDataFrame", pd), mode="onDisk")
      
      #ms1 = all scans MS find in the mzML file of MS
      ms1 <- raw_data@featureData@data[raw_data@featureData@data$msLevel==1,]$acquisitionNum
      ms1match <- mrdf[mrdf$msLevel==2,]$precursorScanNum
      ms2 <- mrdf[mrdf$msLevel==2,]$seqNum
      ms2list <- cbind(ms2,ms1match)

      #Searching for the pre and post MS scan(s)
      prec_scans <- plyr::alply(ms2list, 1, function(scan2){
        post <- ms1[ms1>=scan2[[2]]]
        pre <- ms1[ms1<=scan2[[2]]]
        
        if(length(post)>num){
          if(post[1:num] == pre){
            post <- post[1:num+1]
          }else{
            post <- post[1:num]
          }
        }
        if(length(pre)>num){
          #Problem : post and pre can be the same when we arrive at the last MS1 scan in a different file.
          #For example scan MS2 : 9620, precScan MS1 : 10902 => pre : 10902 and post : 10902
          if(pre[(length(pre)-(num-1)):length(pre)] != post){
            pre <- pre[(length(pre)-(num-1)):length(pre)]
          }else{
            prepre <- pre[(length(pre)-(num)):length(pre)]
            pre <- prepre[1]
          }
        }
        #Searching for the nearest MS scan
        nearest <- ms1[which.min(abs(ms1 - scan2[2]))]
        
        return(list("pre"=pre, "post"=post, "scan"=scan2[[1]], "nearest"= nearest))
      })
    }
  }
  return(prec_scans)
}

get_isolation_offsets <- function(inputfile){

  con  <- file(inputfile, open = "r")
  lowFound = FALSE
  highFound = FALSE

  while (TRUE) {
    oneLine <- readLines(con, n = 1)

    if (!lowFound){
      low <- as.numeric(stringr::str_match(oneLine, '^.*name=\"isolation window lower offset\" value=\"([0-9]+\\.[0-9]+).*$')[,2])
      if(!is.na(low)){lowFound=TRUE}
    }

    if (!highFound){
      high <- as.numeric(stringr::str_match(oneLine, '^.*name=\"isolation window upper offset\" value=\"([0-9]+\\.[0-9]+).*$')[,2])
      if(!is.na(high)){highFound=TRUE}
    }

    if (grepl('<cvParam cvRef="MS" accession="MS:1000490" name="Agilent instrument model" value=""/>', oneLine)){
      message("Agilent do not have isolation offset in mzML, default given as +- 0.65 (this is an approximate for the 'narrow' window type)")
      return(c(0.65, 0.65))
    }

    if(lowFound & highFound){
      break
    }
  }

  close(con)
  return(c(low, high))

}

# Get the Data from one file
getmrdf <- function(filepathMS2, filepathsMS1 = NULL, CSVfile = NULL, backend='pwiz'){
  print("Processing the file....")
  #requireNamespace('mzR') # problem with cpp libraries
  # need to be loaded here for parallel
  mrdf <- NULL
  all_file <- TRUE
  #does this "for" function useless? cause we do MS2 in each file maybe in parallel????????????????????????????????????????????????????????????????????????
  for(i in 1:length(filepathMS2)){
    mr <- mzR::openMSfile(filepathMS2[i], backend=backend)
    mrdfn <- mzR::header(mr)
    cat("============= STEP 1 : Find msLevel informations =============\n")
    # Verify if there is more than 1 level in the file
    if(1 %in% unique(mrdfn$msLevel)){
      if(2 %in% unique(mrdfn$msLevel)){
        print("File with MS1 and MS2 data")
        specin <- 12
      }else{
        print("File with only MS1 data, go to the next file")
        specin <- 1
      }
    }else{
      if(2 %in% unique(mrdfn$msLevel)){
        print("File with only MS2 data")
        specin <- 2
      }
    }
    if(specin==1){
      next
    }

    cat("============= STEP 2 : Find the precursors scans =============\n")
    if(length(unique(mrdfn$precursorScanNum))<2){
      # Note: will be of length 1 even if no scans associated because
      # the mrdf will be zero for not assigned
      print("MS2 data has no associated scan data, will use most recent full scan for information")
      if(specin == 12){
        print("MS1 and MS2 are in the same file")
        mrdfn <- missing_prec_scan(mrdfn)
      }else if(specin == 2){
        print("MS1 are in separated file(s)")
        #Vérification qu'on a bien donné une liste de MS1 et un CSV
        if(is.null(CSVfile) | is.null(filepathsMS1)){
          print("You have to complete the MS1 list and the CSV file to match files MS with MSMS")
          all_file = FALSE
        }else{
          print("C'est good on a tout")
          print(filepathMS2)
          mrdfn <- find_scanMS_for_MS2(mrdfn, filepathsMS1, filepathsMS2 = filepathMS2, CSVfile)
        }
        
      }
    }else{
      print("MS2 already have their precursorScanNum")
    }

    if(all_file == TRUE){
      cat("============= STEP 3 : Complete MS2s with precursors informations =============\n")
      if(specin == 12){
        print("MS1 are in the same file")
        #mrdfn$fileid <- rep(i,nrow(mrdfn))
        mrdfn$filename <- rep(basename(filepathMS2[i]),nrow(mrdfn))
        mrdfn$precursorRT <- NA
        # precursorScanNum matches to the acquisitionNum, get row matching row number and relevant retention time
        mrdfn[mrdfn$msLevel==2,]$precursorRT <- mrdfn[match(mrdfn[mrdfn$msLevel==2,]$precursorScanNum, mrdfn$acquisitionNum),]$retentionTime
      }else if(specin == 2){
        print("MS1 are from other files")
        #mrdfn$fileid <- rep(i,nrow(mrdfn))
        mrdfn$precursorRT <- NA
        fileToLoad <- unique(mrdfn$precursorFilename)
        s_groups <- sapply(fileToLoad, function(x) tail(unlist(strsplit(dirname(x),"/")), n=1))
        s_name <- tools::file_path_sans_ext(basename(fileToLoad))
        pd <- data.frame(sample_name=s_name, sample_group=s_groups, stringsAsFactors=FALSE)
        raw_data <- MSnbase::readMSData(files=fileToLoad, pdata = new("NAnnotatedDataFrame", pd), mode="onDisk")
        for(MSMSdata in 1:nrow(mrdfn)){
          mrdfn[MSMSdata,"filename"] <- basename(filepathMS2[i])
          datafind <- raw_data@featureData@data[raw_data@featureData@data[,"acquisitionNum"] == mrdfn[MSMSdata,"precursorScanNum"],]
          mrdfn[MSMSdata,"precursorRT"] <- datafind$retentionTime
        }
      }

      cat("============= STEP 4 : Completing output =============\n")
      if(!is.data.frame(mrdf)){
        mrdf <- mrdfn
      }else{
        mrdf <- rbind(mrdf,mrdfn)
      }
    }else{
      return(all_file)
    }
    
  }
  return(mrdf)
}

find_scanMS_for_MS2 <- function(mrdfn, filepathsMS1 = NULL, filepathsMS2 = NULL, CSVfile = NULL){
  print("Find the MS1 scans for each MS2scan")
 
  print(paste("MS1 files :",filepathsMS1))
  print(paste("MS2 file :",filepathsMS2,sep=""))
  #Find the good MS file corresponding to the MSMS file
  print(paste("CSV matching files filename :",CSVfile))
  #Explore the CSV file to find the MSMS filename
  readCSV<-basename(CSVfile)
  fileToMatch <- read.csv2(file=readCSV, header=FALSE)
  names(fileToMatch) <- c("MS1","MS2")
  fileToBeMatch <- names(filepathsMS2)
  for(row in 1:nrow(fileToMatch)){
    if(fileToMatch[row,2] == basename(fileToBeMatch)){
      MS1matching <- fileToMatch[row,1]
    }
  }

  #Find the complete filename of the MS file matched
  for(i in 1:length(filepathsMS1)){
    if(names(filepathsMS1) == MS1matching){
      fileToLoad <- paste("./",names(filepathsMS1),sep="")
      print(fileToLoad)
    }
  }

  #Read the MS1 file
  s_groups <- sapply(fileToLoad, function(x) tail(unlist(strsplit(dirname(x),"/")), n=1))
  s_name <- tools::file_path_sans_ext(basename(fileToLoad))
  pd <- data.frame(sample_name=s_name, sample_group=s_groups, stringsAsFactors=FALSE)

  print(paste("Reading the MS file ",fileToLoad,"...",sep=""))

  raw_data <- MSnbase::readMSData(files=fileToLoad, pdata = new("NAnnotatedDataFrame", pd), mode="onDisk")
  # Transform the files absolute pathways into relative pathways
  raw_data@processingData@files <- sub(paste(getwd(), "/", sep="") , "", raw_data@processingData@files)

  #Filter msLevel1 if the file has also some MS2
  raw_data <- MSnbase::filterMsLevel(raw_data,msLevel=1)

  #Stock RTs of each MS scan
  vectorMS1RT <- c()
  for(i in 1:nrow(raw_data@featureData@data)){
    vectorMS1RT[i]<-raw_data@featureData@data[i,"retentionTime"]
  }

  for(i in 1:nrow(mrdfn)){
    if(mrdfn[i,]$msLevel>1){
      diff <- abs(vectorMS1RT-mrdfn[i,"retentionTime"])
      #Obtain the scan numbers for each (MS2 and peak)
      numscanMS2 <- which(diff==min(diff)) #findInterval(MS2data["retentionTime"],vectorMS1RT)

      if(numscanMS2 != mrdfn[i,"precursorScanNum"]){
        mrdfn[i,"precursorScanNum"] = raw_data@featureData@data[numscanMS2,"acquisitionNum"]
      }
      
      mrdfn[i,"precursorFilename"] = basename(fileToLoad)
    }
  }
  return(mrdfn)
}

#This function find all MS1 scans (x) until they arrive at the scan number of the MS2 scans that we are using (i). The last MS1 scan is saved.
#We break if we have the same scan number for i and x (means that it is a MS1 scan).
missing_prec_scan <- function(mrdfn){
  for(i in 1:nrow(mrdfn)){
    if(mrdfn[i,]$msLevel>1){
      # Find most recent ms1 level scan
      for (x in i:1){
        if(mrdfn[x,]$msLevel==1){
          mrdfn[i,]$precursorScanNum = mrdfn[x,]$acquisitionNum
          break
        }
      }
    }
  }
  return(mrdfn)
}

getscans <- function(files, backend='pwiz'){
  if(length(files)==1){
    mr <- mzR::openMSfile(files, backend=backend)
    scan_peaks <- mzR::peaks(mr)
    return(scan_peaks)
  }else{
    scan_peaks <- plyr::alply(files, 1 ,function(x){
      print(x)
      mr <- mzR::openMSfile(x, backend=backend)
      scan_peaks <- mzR::peaks(mr)
      return(scan_peaks)
    })
    return(scan_peaks)
  }
}

# MSMSperMS <- function(filepths){
#   l <- lapply(filepths, MSMSperMSsingle)
#   return(plyr::ldply(l, summary))
# }
#
#
# MSMSperMSsingle <- function(filepth){
#   mr <- mzR::openMSfile(filepth)
#   mrdf <- mzR::header(mr)
#
#   m <- mrdf$msLevel
#
#   # First remove repeating ms1 elements
#   c <- 1
#   idx <- c()
#   for (i in 1:length(m)){
#
#     if(i==length(m)){
#       break
#     }
#     if(m[i]==1 && m[i+1]==1) {
#       idx[c] <- i
#       c = c+1
#     }
#   }
#   if(!is.null(idx)){
#     m <- m[-idx]
#   }
#
#
#   # then count how many ms2 are present at each ms1---ms2ms2...ms2---ms1 interval
#   ms2all <- c()
#   c <- 0
#   cidx <- 1
#   for(i in 1:length(m)){
#
#     if(m[i]==2) {
#       c = c+1
#     }
#
#     if(m[i]==2 && (m[i+1]==1 || is.na(m[i+1]))) {
#       ms2all[cidx] <- c
#       cidx = cidx+1
#       c = 0
#     }
#
#   }
#
#   return(ms2all)
#
# }

plotPurity <- function(df, inter, name, f){

  dfm <- reshape2::melt(df, id.vars = c("purity", "idx"))
  png(filename = name)
  p1 <- ggplot2::ggplot(dfm, ggplot2::aes(x=idx, y=purity)) +
    ggplot2::geom_point(ggplot2::aes(colour=variable),shape=20, colour="red") +
    ggplot2::geom_vline(xintercept = inter, colour="red", linetype = "longdash")+
    ggplot2::theme_bw(base_size = 18)+
    ggplot2::theme(text = ggplot2::element_text(size=23))+
    ggplot2::labs(x="scan" , y="precursor purity score") +
    ggplot2::stat_function(fun=f, colour="blue")+
    ggplot2::ylim(0, 1)+
    ggplot2::scale_color_manual(
            "Legend Title\n",labels = c("EIC feature", "+/- 0.5 Contamination"),
             values = c("red", "black") )
  print(p1)
  dev.off()
  dfm$pick <- inter
  saveRDS(f, paste(name, ".rds", sep=""))
  write.csv(dfm, paste(name, ".csv", sep=""), row.names = FALSE)
}