#@authors ABiMS TEAM, J. Saint-Vanne
# lib.r for MSMS study

#@author G. Le Corguille
# This function will
# - load the packages
# - display the sessionInfo
loadAndDisplayPackages <- function(pkgs) {
    for(pkg in pkgs) suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))

    sessioninfo = sessionInfo()
    cat(sessioninfo$R.version$version.string,"\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) { cat(paste(pkg,packageVersion(pkg)),"\t") }; cat("\n")
}

# This function get the raw file path from the arguments
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
getRawfilePathFromArguments <- function(singlefile, zipfile, args, prefix="") {
  if (!(prefix %in% c("","Positive","Negative","MS1","MS2"))) stop("prefix must be either '', 'Positive', 'Negative', 'MS1' or 'MS2'")

  if (!is.null(args[[paste0("zipfile",prefix)]])) zipfile <- args[[paste0("zipfile",prefix)]]

  if (!is.null(args[[paste0("singlefile_galaxyPath",prefix)]])) {
    singlefile_galaxyPaths <- args[[paste0("singlefile_galaxyPath",prefix)]]
    singlefile_sampleNames <- args[[paste0("singlefile_sampleName",prefix)]]
  }
  if (exists("singlefile_galaxyPaths")){
    singlefile_galaxyPaths <- unlist(strsplit(singlefile_galaxyPaths,"\\|"))
    singlefile_sampleNames <- unlist(strsplit(singlefile_sampleNames,"\\|"))

    singlefile <- NULL
    for (singlefile_galaxyPath_i in seq(1:length(singlefile_galaxyPaths))) {
      singlefile_galaxyPath <- singlefile_galaxyPaths[singlefile_galaxyPath_i]
      singlefile_sampleName <- singlefile_sampleNames[singlefile_galaxyPath_i]
      singlefile[[singlefile_sampleName]] <- singlefile_galaxyPath
    }
  }
  return(list(zipfile=zipfile, singlefile=singlefile))
}

# This function retrieve the raw file in the working directory
#   - if zipfile: unzip the file with its directory tree
#   - if singlefiles: set symlink with the good filename
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr
retrieveRawfileInTheWorkingDirectory <- function(singlefile, zipfile) {
    if(!is.null(singlefile) && (length("singlefile")>0)) {
        for (singlefile_sampleName in names(singlefile)) {
            singlefile_galaxyPath <- singlefile[[singlefile_sampleName]]
            if(!file.exists(singlefile_galaxyPath)){
                error_message <- paste("Cannot access the sample:",singlefile_sampleName,"located:",singlefile_galaxyPath,". Please, contact your administrator ... if you have one!")
                print(error_message); stop(error_message)
            }

            if (!suppressWarnings( try (file.link(singlefile_galaxyPath, singlefile_sampleName), silent=T)))
                file.copy(singlefile_galaxyPath, singlefile_sampleName)

        }
        directory <- "."

    }
    if(!is.null(zipfile) && (zipfile != "")) {
        if(!file.exists(zipfile)){
            error_message <- paste("Cannot access the Zip file:",zipfile,". Please, contact your administrator ... if you have one!")
            print(error_message)
            stop(error_message)
        }

        #list all file in the zip file
        #zip_files <- unzip(zipfile,list=T)[,"Name"]

        #unzip
        suppressWarnings(unzip(zipfile, unzip="unzip"))

        #get the directory name
        suppressWarnings(filesInZip <- unzip(zipfile, list=T))
        directories <- unique(unlist(lapply(strsplit(filesInZip$Name,"/"), function(x) x[1])))
        directories <- directories[!(directories %in% c("__MACOSX")) & file.info(directories)$isdir]
        directory <- "."
        if (length(directories) == 1) directory <- directories

        cat("files_root_directory\t",directory,"\n")

    }
    return (directory)
}

#This function list the compatible files within the directory as xcms did
#@author Gildas Le Corguille lecorguille@sb-roscoff.fr ABiMS TEAM
getMSFiles <- function (directory) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]","[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep=""),collapse="|")
    info <- file.info(directory)
    listed <- list.files(directory[info$isdir], pattern=filepattern,recursive=TRUE, full.names=TRUE)
    files <- c(directory[!info$isdir], listed)
    exists <- file.exists(files)
    files <- files[exists]
    return(files)
}

# Just to verify the msPurity version
use <- function(package, version=0, ...) {
  package <- as.character(substitute(package))
  library(package, ..., character.only=TRUE)
  pver <- packageVersion(package)
  if (compareVersion(as.character(pver), as.character(version)) < 0)
    stop("Version ", version, " of '", package, 
         "' required, but only ", pver, " is available")
  invisible(pver)
}

loadRData <- function(rdata_path, name){
#loads an RData file, and returns the named xset object if it is there
    load(rdata_path)
    return(get(ls()[ls() %in% name]))
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