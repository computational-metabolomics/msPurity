#' @title Using Subtract MZ values based on ppm tolerance and noise ratio
#'
#' @description
#' Uses a purityD object with references to multiple MS files. Subtract blank peaks from the sample peaks
#' see subtractMZ for more information
#'
#' @aliases subtract
#'
#' @param Object object; purityD object
#' @param mapping parameter not functional (TODO)
#' @param byClass boolean; subtract within each class
#' @inheritParams subtractMZ
#' @examples
#' datapth <- system.file("extdata", "dims", "mzML", package="msPurityData")
#' inDF <- Getfiles(datapth, pattern=".mzML", check = FALSE, cStrt = FALSE)
#'
#' ppDIMS <- purityD(inDF, cores=1)
#' ppDIMS <- averageSpectra(ppDIMS)
#' ppDIMS <- filterp(ppDIMS, thr = 5000)
#' ppDIMS <- subtract(ppDIMS)
#' @return  purityD object with averaged spectra
#' @seealso \code{\link{subtractMZ}}
#' @export
setMethod(f="subtract", signature = "purityD",
          definition = function(Object, byClass = TRUE, mapping=c("sample", "blank"), ppm = 5, s2bthres=10){
  fileList <- Object@fileList

  requireNamespace('foreach')

  if (Object@cores>1){
    operator <- foreach::'%dopar%'
    cl<-parallel::makeCluster(Object@cores, type = "SOCK") #change the 2 to your number of CPU cores
    doSNOW::registerDoSNOW(cl)
  }else{
    operator <- foreach::'%do%'
  }

  processedP <- operator(foreach::foreach(i=1:nrow(fileList)), subtractSingle(Object, i, ppm, s2bthres))

  # Update the object
  Object@avPeaks$processed <- processedP
  names(Object@avPeaks$processed) <- fileList$name

  return(Object)
})

subtractSingle <- function(Object, i, ppm, s2bthres){

  fileList <- Object@fileList

  file <- fileList[i, ]

  if(file$sampleType=="blank"){
    # we dont change the blank peaks
    return(Object@avPeaks$processed[[i]])
  }

  # get the blank and sample idx for the peaklist
  files <- fileList[fileList$class==file$class,]
  blankIdx <- as.numeric(rownames(files[files$sampleType=='blank',]))
  sampleIdx <- as.numeric(rownames((files[files$sampleType=='sample',])))

  # Get the peaks (both processed and original)
  avPeaks <- Object@avPeaks
  # Extract the blank and sample mz and i values
  blanMz <- avPeaks$processed[[blankIdx]]$mz
  blani <- avPeaks$processed[[blankIdx]]$i

  samplePeaks <- avPeaks$processed[[sampleIdx]]
  sampMz <- samplePeaks$mz
  sampi <- samplePeaks$i

  if(nrow(samplePeaks)>0){
    #print(head(blanMz))
    #print(head(sampMz))

    # subtract the mz values within a ppm tolerance
    mzsub <- subtractMZ(mz1 = sampMz, mz2 = blanMz, i1 = sampi, i2 = blani, s2bthres = s2bthres, ppm = ppm)

    # Filter out any peaks that have been subtracted
    if(!is.null(mzsub)){
      subtractedP <- samplePeaks[samplePeaks$mz %in% mzsub,]
    }else{
      subtractedP <- samplePeaks
    }
  }else{
    subtractedP <- samplePeaks
  }

  return(subtractedP)

}



#' @title Subtract MZ values based on ppm tolerance and noise ratio
#'
#' @description
#' This function is intended for blank subtraction of mz values from two peaklists. It takes in 2 vectors of mz values and 2
#' coresponding vectors of Intensity values.
#'
#' The second mz values are subtracted from the first set within an MZ tolerance.
#'
#' However, if the mz match but the intensity is above a defined threshold then they are not subtracted
#'
#' @param mz1 vector = mz values to start with
#' @param mz2 vector = mz values to subtract
#' @param i1 vector = i values for mz1
#' @param i2 vector = i values for mz2
#' @param s2bthres numeric = threshold for the samp2blank (i1/i2)
#' @param ppm numeric = ppm tolerance
#' @return a vector of the remaining mz values
#'
#' @examples
#' mz1 <- c(100.001, 200.002, 300.302)
#' mz2 <- c(100.004, 200.003, 500.101)
#' i1 <- c(100, 100, 100)
#' i2 <- c(100, 10000, 100)
#'
#' subtractMZ(mz1, mz2, i1, i2, ppm=5, s2bthres =10)
#' @export
subtractMZ <- function(mz1, mz2, i1, i2, ppm = 5, s2bthres = 10){

  # Get ranges of the "sample"
  mz1Lo = round(mz1 - ((mz1*0.000001)*ppm), 10)
  mz1Up = round(mz1 + ((mz1*0.000001)*ppm), 10)

  # get ranges of the "blank" to remove
  mz2Lo = round(mz2 - ((mz2*0.000001)*ppm), 10)
  mz2Up = round(mz2 + ((mz2*0.000001)*ppm), 10)

  removemz <-  vector()

  for (i in 1:length(mz1)){
    for(j in 1:length(mz2)){
      if(overlap(mz1Lo[i], mz1Up[i], mz2Lo[j], mz2Up[j])){
        # if overlap is true then remove this mz value
        s2b <- i1[i]/i2[j]

        if(s2b<s2bthres){
          removemz[i] <- i
          break
        }

      }
    }
  }

  removemz <- removemz[!is.na(removemz)]
  newmz <- mz1[-removemz]

  if(length(newmz)==0){
    newmz <- NULL
  }

  return(newmz)

}

overlap <- function(start1, end1, start2, end2){
  #Does the range (start1, end1) overlap with (start2, end2)? (returns boolean)
  #based on "De Morgan's laws", see http://nedbatchelder.com/blog/201310/range_overlap_in_two_compares.html
  return((end1 >= start2) & (end2 >= start1))

}
