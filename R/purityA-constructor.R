#' @include iw-norm.R pcalc.R
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
#' @param fileList vector = mzML file paths for MS/MS spectra
#' @param cores numeric = number of cores to use
#' @param mostIntense boolean = True if the most intense peak is used for calculation. False if the centered peak is used
#' @param nearest boolean = True if the peak selected is from either the preceding scan or the nearest.
#' @param offsets vector = overide the isolation offsets found in the mzML filee.g. c(0.5, 0.5)
#' @param plotP boolean = if TRUE a plot of the purity is to be saved
#' @param plotdir vector = if plotP is TRUE plots will be saved to this directory
#' @param interpol character = type of interolation to be performed "linear" or "spline"
#' @param iwNorm boolean = if TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function = A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric = All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5\% (0.05)
#' @param mzRback character = backend to use for mzR parsing
#'
#' @return a dataframe of the purity score of the ms/ms spectra
#'
#' @examples
#' filepths <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")
#' pa <- purityA(filepths)
#' @seealso \code{\link{assessPuritySingle}}
#' @export
purityA <- function(fileList,
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
                    mzRback='pwiz'){

  if((length(fileList)>=1) && (fileList == "" )){
    message("no file list")
    return(NULL)
  }

  requireNamespace('foreach')
  pa <- new("purityA", fileList = fileList , cores = cores, mzRback=mzRback)
  
  # Check cores and choose if parallel or not (do or dopar)
  if(pa@cores<=1){
    operator <- foreach::'%do%'
  } else{
    cl1 <- parallel::makeCluster(pa@cores)
    doSNOW::registerDoSNOW(cl1)
    operator <- foreach::'%dopar%'
  }

  # if iwNorm is TRUE and iwNormFun is NULL
  if(is.null(iwNormFun)){
    iwNormFun <- iwNormGauss()
  }

  # run parallel (or not) using foreach
  purityL <- operator(foreach::foreach(i = 1:length(pa@fileList),
                                  .packages = 'mzR'),
                                  assessPuritySingle(filepth = pa@fileList[[i]],
                                  mostIntense = mostIntense,
                                  nearest=nearest,
                                  offsets = offsets,
                                  plotP = plotP,
                                  plotdir = plotdir,
                                  interpol = interpol,
                                  iwNorm = iwNorm,
                                  iwNormFun = iwNormFun,
                                  ilim = ilim,
                                  mzRback = mzRback
                                  ))

  if(pa@cores>1){
      parallel::stopCluster(cl1)
  }

  names(purityL) <- seq(1, length(purityL))

  puritydf <- plyr::ldply(purityL, .id = TRUE)

  if(nrow(puritydf)>0){
    colnames(puritydf)[1] <- "fileid"
  }
  pid <- seq(1, nrow(puritydf))
  puritydf <- cbind(pid, puritydf)

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
#' @param filepth character = mzML file path for MS/MS spectra
#' @param fileid numeric = adds a fileid column (primarily for internal use for msPurity)
#' @param mostIntense boolean = True if the most intense peak is used for calculation. False if the centered peak is used
#' @param nearest boolean = True if the peak selected is as the nearest MS1 scan. If False then the preceding scan is used
#' @param offsets vector = overide the isolation offsets found in the mzML filee.g. c(0.5, 0.5)
#' @param cores numeric = number of cores to use
#' @param plotP boolean = if TRUE a plot of the purity is to be saved
#' @param plotdir vector = if plotP is TRUE plots will be saved to this directory
#' @param interpol character = type of interolation to be performed "linear", "spline" or "none"
#' @param iwNorm boolean = if TRUE then the intensity of the isolation window will be normalised based on the iwNormFun function
#' @param iwNormFun function = A function to normalise the isolation window intensity. The default function is very generalised and just accounts for edge effects
#' @param ilim numeric = All peaks less than this percentage of the target peak will be removed from the purity calculation, default is 5\% (0.05)
#' @param mzRback character = backend to use for mzR parsing
#' @return a dataframe of the purity score of the ms/ms spectra
#'
#' @examples
#' filepth <- system.file("extdata", "lcms", "mzML", "LCMSMS_1.mzML", package="msPurityData")
#'
#' puritydf <- assessPuritySingle(filepth)
#' @seealso \code{\link{purityA}}
#' @export
assessPuritySingle <- function(filepth,
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
                               mzRback='pwiz'){
  #=================================
  # Load in files and initial setup
  #=================================
  # Get the mzR dataframes
  mrdf <- getmrdf(filepth, mzRback)
  if(is.null(mrdf)){
    message(paste("No MS/MS spectra for file: ", filepth))
    return(NULL)
  }

  # having id column makes things easier to track
  mrdf$id <- seq(1, nrow(mrdf))

  # add the fileid info
  if(!is.na(fileid)){
    mrdf$fileid <- rep(fileid, nrow(mrdf))
  }

  # add filename
  mrdf$filename <- basename(filepth)

  # get scans (list of mz and i) from mzR
  scans <- getscans(filepth, mzRback)

  # Get offsets from mzML unless defined by user
  if(anyNA(offsets)){
    offsets <- get_isolation_offsets(filepth)
  }
  minoff <- offsets[1]
  maxoff <- offsets[2]

  # Get a shortened mzR dataframe
  mrdfshrt <- mrdf[mrdf$msLevel==2,][,c("seqNum","precursorIntensity",
                                        "precursorMZ", "precursorRT",
                                        "precursorScanNum", "id", "filename")]

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


  # if iwNorm is TRUE and iwNormFun is NULL
  # then a gaussian model of the isolation window will be used to normalise
  # intensity
  if(is.null(iwNormFun)){
    # Using a gaussian curve 3 SD either side
    iwNormFun <- iwNormGauss(3)
  }

  #=====================================
  # Initial purity calculation
  #=====================================
  # Get the initial purity score and target mz associated with every ms2 scan.
  # Calculate the purity on the nearest ms1 scan. This function also determines
  # which is the target peak used for interpolations
  initial <- plyr::ddply(mrdfshrt, ~ seqNum, get_init_purity,
                           scans=scans,
                           minoff=minoff,
                           maxoff=maxoff,
                           nearest=nearest,
                           mostIntense=mostIntense,
                           iwNorm=iwNorm,
                           iwNormFun=iwNormFun,
                           ilim=ilim)

  # Add results to mzR dataframe
  mrdfshrt <- cbind(mrdfshrt, initial)

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
                                   scan_peaks=scans,
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
                                   ilim=ilim)

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


get_init_purity <- function(ms2h, scans, minoff, maxoff, nearest,
                            mostIntense, iwNorm, iwNormFun, ilim){

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

  allscans <- scans[[precScn]]

  #=======================================
  # Get centered target peak
  #=======================================
  # Get target based on matching to be within isolation window of precursorMZ
  # and 1 decimal of intensity
  target_approx <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax)
                            & (round(allscans[,2],1) == round(precI,1)),]

  # no precusor within isolation window and itensity
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
    target <- target_approx[which.min(abs(target_approx[,1]- precMZ)),]
  }

  # The centered peak
  aMz <- target[1]

  #==================================================
  # Get most intense peak within isolation window
  #==================================================
  subp <- allscans[(allscans[,1]>=mzmin) & (allscans[,1]<=mzmax),  ]

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
                   iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim)
    aPurity <- pouta[1]
    apkNm <- pouta[2]
    pouti <- pcalc(peaks=subp, mzmin=mzmin, mzmax=mzmax, mztarget=iMz, ppm=NA,
                   iwNorm=iwNorm, iwNormFun=iwNormFun, ilim=ilim)
    iPurity <- pouti[1]
    ipkNm <- pouti[2]

  }
  fileinfo <- c("aMz"=aMz, "aPurity" = aPurity, "apkNm" = apkNm,
                "iMz" = iMz, "iPurity" = iPurity, "ipkNm" = ipkNm )
  return(fileinfo)

}

get_interp_purity <- function(rowi, scan_peaks, prec_scans, ms2, ppm,
                              minoff, maxoff, mostIntense=FALSE, plotP=FALSE,
                              plotdir=NULL, interpol, nearest=TRUE,
                              iwNorm, iwNormFun, ilim){

  # get the region of interest for scans
  scanids <- prec_scans[which(ms2==rowi$seqNum)][[1]]
  roi_scns <- scan_peaks[c(scanids$pre, scanids$post)]

  if(interpol=="linear"){
    purity <- linearPurity(rowi, scan_peaks, minoff, maxoff, ppm, scanids,
                           nearest, mostIntense, iwNorm, iwNormFun,
                           ilim, plotP, plotdir)
  }else if (interpol=="spline"){
    purity <- splinePurity(rowi, roi_scns, minoff, maxoff, ppm,
                           mostIntense, scanids, plotP, plotdir)
  }


  return(purity)

}

linearPurity <- function(rowi, scan_peaks, minoff, maxoff, ppm, scanids,
                         nearest, mostIntense, iwNorm, iwNormFun, ilim,
                         plotP, plotdir){
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
                targetMinMZ = targetMinMZ, targetMaxMZ = targetMaxMZ)

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


# Get closest precursor MZ value (or most intense)
get_prec_scans <- function(mrdf, num){

  ms1 <- mrdf[mrdf$msLevel==1,]$seqNum
  ms2 <- mrdf[mrdf$msLevel==2,]$seqNum

  prec_scans <- plyr::alply(ms2, 1, function(scan2){
    post <- ms1[ms1>scan2]
    pre <- ms1[ms1<scan2]

    if(length(pre)>num){
      pre <- pre[(length(pre)-(num-1)):length(pre)]
    }

    if(length(post)>num){
      post <- post[1:num]
    }

    nearest <- ms1[which.min(abs(ms1 - scan2))]

    return(list("pre"=pre, "post"=post, "scan"=scan2, "nearest"= nearest))
  })

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


# Get the Data
getmrdf <- function(files, backend='pwiz'){
  #requireNamespace('mzR') # problem with cpp libraries
  # need to be loaded here for parallel
  mrdf <- NULL

  for(i in 1:length(files)){
    #message(paste("processing file:" ,i))
    mr <- mzR::openMSfile(files[i], backend=backend)
    mrdfn <- mzR::header(mr)
    if(length(unique(mrdfn$msLevel))<2){
      #message("only MS1 data")
      next
    }

    if(length(unique(mrdfn$precursorScanNum))<2){
      # Note: will be of length 1 even if no scans associated because
      # the mrdf will be zero for not assigned
      message("MS2 data has no associated scan data, will use most recent full scan for information")
      mrdfn  <- missing_prec_scan(mrdfn)
    }

    # mrdfn  <- missing_prec_scan(mrdfn)
    if(length(files)==1){

    }
    #mrdfn$fileid <- rep(i,nrow(mrdfn))
    mrdfn$filename <- rep(basename(files[1]),nrow(mrdfn))
    mrdfn$precursorRT <- NA
    mrdfn[mrdfn$msLevel==2,]$precursorRT <- mrdfn[mrdfn$precursorScanNum,]$retentionTime

    if(!is.data.frame(mrdf)){
      mrdf <- mrdfn
    }else{
      mrdf <- rbind(mrdf,mrdfn)
    }
  }
  return(mrdf)
}

missing_prec_scan <- function(mrdfn){
  for(i in 1:nrow(mrdfn)){
    if(mrdfn[i,]$msLevel>1){
      # Find most recent ms1 level scan
      for (x in i:1){
        if(mrdfn[x,]$msLevel==1){
          mrdfn[i,]$precursorScanNum = mrdfn[x,]$seqNum
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
