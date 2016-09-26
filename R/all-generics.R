setGeneric("getP", function(x) standardGeneric("getP"))

setGeneric("updatePeaks", function(x, newlist) standardGeneric("updatePeaks"))

setGeneric(name="averageSpectra",
           def=function(Object,rtscn = "all", scanRange=NA, timeRange = NA,
                        clustType="hc", ppm=1.5, snthr = 3, av="median",
                        missingV="zero", minfrac=0.6667,
                        normTIC=FALSE, snMeth="median", mzRback="pwiz") {
             standardGeneric("averageSpectra")
           }
)


setGeneric(name="groupPeaks",
           def=function(Object,ppm=3, sampleOnly=FALSE, clustType='hc') {
             standardGeneric("groupPeaks")
           }
)


setGeneric(name="filterp",
           def=function(Object, thr = 5000, rsd=20, sampleOnly = TRUE) {
             standardGeneric("filterp")
           }
)

setGeneric(name="subtract",
           def=function(Object, byClass = TRUE, mapping=c("sample", "blank"), ppm =5,  s2bthres=10){
             standardGeneric("subtract")
           }
)

setGeneric(name="dimsPredictPurity",
           def=function(Object, ppm = 1.5, minOffset=0.5, maxOffset=0.5,
                        iwNorm=FALSE, iwNormFun=NULL, ilim=0.05, sampleOnly=TRUE) {
             standardGeneric("dimsPredictPurity")
           }
)

setGeneric(name="writeOut",
           def=function(Object, outDir, original) {
             standardGeneric("writeOut")
           }
)

setGeneric(name="frag4feature",
           def=function(pa, xset, ppm = 5, plim = 0, intense=TRUE, convert2RawRT=TRUE){
             standardGeneric("frag4feature")
           }
)

setGeneric(name="validate",
           def=function(pa, ppLCMS){
             standardGeneric("validate")
           }
)


