#' @title Get files for DI-MS processing
#'
#' @description
#' Takes in a folder path and outputs the a data frame structure for purityPD.
#' Function modified from mzmatch.
#'
#' @param projectFolder character: directory path
#' @param recursive boolean: recursively check for files
#' @param pattern character file suffix to check for
#' @param check boolean check with a GUI the files
#' @param cStrt boolean use the first word as the class name for files
#' @param raw (REDUNDANT)
#' @param peakout (REDUNDANT)
#' @param mzml_out (REDUNDANT)
#' @return dataframe of files
#' @examples
#'
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#' @export
Getfiles = function(projectFolder=NULL, recursive = FALSE, pattern = ".csv",
                    check = TRUE, raw = FALSE, peakout = NA, cStrt = TRUE,
                    mzml_out=FALSE){


  #If no list provided, search Path for files.
  if (is.null(projectFolder)) {
    #Pop up box will appear to prompt selection of directory by user.
    projectFolder <- tcltk::tclvalue(tcltk::tkchooseDirectory())
  }
  files.fullnames <- NULL

  ##Get file information from selected directory
  HIT <- dir(projectFolder, full.names=TRUE, pattern=pattern,
             ignore.case = TRUE, recursive=recursive,include.dirs = FALSE)

  if(raw){
    sampleList <- data.frame(RAW=HIT)
  }else{
    sampleList <- data.frame(filepth=HIT)
  }

  for(i in 1:length(HIT)){

    # GET NAME ID
    name <- strsplit(basename(HIT[i]), pattern)[[1]]
    sampleList[i, 'name'] = name

    # GET POLARITY ID
    if (grepl('pos', HIT[i]) == TRUE){
      sampleList[i, 'polarity'] = 'pos'
    } else if (grepl('neg', HIT[i]) == TRUE){
      sampleList[i, 'polarity'] = 'neg'
    }

    # GET SAMPLE TYPE
    if (grepl('*Blank*', HIT[i], ignore.case = TRUE) == TRUE){
      sampleList[i, 'sampleType'] = "blank"
    } else {
      sampleList[i, 'sampleType'] = "sample"
    }

    # GET CLASS
    ns <- strsplit(name, "_")[[1]]
    if(cStrt){
      sampleList[i, 'class'] = ns[1]
    }else{
      sampleList[i, 'class'] = ns[length(ns)]
    }


    if(raw){
      if(mzml_out){
        sampleList[i, 'filepth'] = file.path(peakout, paste(name, ".mzML", sep=""))
      }else{
        sampleList[i, 'filepth'] = file.path(peakout, paste(name, ".csv", sep=""))
      }

    }


  }

  sampleList <- sampleList[order(sampleList$sampleType), ]
  sampleList <- sampleList[order(sampleList$class), ]

  rownames(sampleList) <- seq(1, nrow(sampleList))

  if(check){
    fix(sampleList)
    ## fix writes edited object to the .GlobalEnv,
    # we have to read it back into function environmnet.
    sampleList <- get("sampleList", envir=.GlobalEnv)
  }


  return(sampleList)
}
