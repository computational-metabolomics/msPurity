# peaksPmzSingle <- function(row, mzw=1, ilimit=0){
#   mzmlPth <- as.character(row$filepth)
#   mr <- mzR::openMSfile(as.character(mzmlPth))
#
#   # get ms1 peaks
#   scanPeaks <- mzR::peaks(mr)
#   h <- header(mr)
#   hms1 <- h[h$msLevel==1,]
#   scanList <- scanPeaks[hms1$seqNum]
#
#   # bin into windows for each scan
#
#   names(scanList) <- seq(1, length(scanList))
#
#   counts <- ldply(scanList, .id='id', cutup4rawp, mzw=mzw, ilimit=ilimit)
#   count_hm <- dcast(counts, .id ~ id, value.var="V1")
#
#   avC <- ddply(count_hm, .(.id), function(x){
#     median(as.numeric(x[-1]))
#     })
#
#   return(median(avC$V1))
#
# }
#
#
#
# spectralComplexityExp <- function(exP, mzw=10, ilimit=1000){
#
#   sp <- ddply(exP@fileList[exP@sampleIdx,], .(name), peaksPmzSingle, mzw=mzw,ilimit=ilimit)
#
#   colnames(sp) <- c("name", "spectralComplexity")
#
#   return(sp)
#
# }
#
#
# cutup4rawp <- function(scani, mzw, ilimit){
#   c1 <- cut(scani[,1], breaks = seq(100, 1000, by = mzw))
#   df.split <- split(scani, c1)
#   count <- ldply(df.split, function(x){length(x[x>ilimit])})
# }
#
#
#
# peaksPscan <- function(files){
#
#   count <- sapply(files, function(x){
#     print(x)
#     mr <- mzR::openMSfile(x)
#     h <- mzR::header(mr)
#     ms1 <- h$seqNum[h$msLevel<2]
#     pc <- mzR::peaksCount(mr, ms1)
#     return(sum(pc)/length(pc))
#   })
#
#   return(summary(count))
#
# }
#
#
# scansPfile <- function(files){
#
#   count <- sapply(files, function(x){
#     mr <- mzR::openMSfile(x)
#     h <- header(mr)
#     ms2 <- nrow(h[h$msLevel>1,])
#     return(ms2)
#   })
#
#
#   s <- summary(count)
#   spf <- sum(count)/length(count)
#   return(c(s, "MS2 count per file"=spf, "file count"=length(count)))
# }
#
#
#
# imsms <- function(files){
#
#   i <- alply(files,1, function(x){
#     mr <- mzR::openMSfile(x)
#     h <- header(mr)
#     ms2 <- h[h$msLevel>1,]$precursorIntensity
#     return(ms2)
#   })
#
#   i <- unlist(i)
#
#   return(summary(i))
# }
#
#
#
# purityPosition <- function(puritydf){
#
#   seqNum <- puritydf$seqNum
#   c=0
#   id <- c()
#
#   for(i in 1:nrow(puritydf)){
#     if(i==1){
#       seqim <- 0
#       id[i] <- 1
#
#       next
#
#     }else{
#       seqim <- seqNum[i-1]
#
#     }
#
#
#     if (i==nrow(puritydf)){
#       seqip <- 0
#     }else{
#       seqip <- seqNum[i+1]
#     }
#
#
#     seqi <- seqNum[i]
#
#     if((seqi+1)==seqip){
#       print(seqi)
#       c = c+1
#       id[i] <- c
#
#     }else if((seqi-1)==seqim){
#       print(seqi)
#       c = c+1
#       id[i] <- c
#       c = 0
#
#
#     }else{
#       print("just1")
#       print(seqi)
#       id[i] <- 1
#       c=0
#     }
#
#   }
#
#
#   puritydf$spa <- id
#
#   summarySE
#   ggplot(puritydf, aes(spa, inPurity)) +geom_line()
#
#   library(ggplot2)
#
# }
