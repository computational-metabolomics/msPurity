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
                        iwNorm=FALSE, iwNormFun=NULL, ilim=0.05, sampleOnly=TRUE,
                        isotopes=TRUE, im=NULL) {
             standardGeneric("dimsPredictPurity")
           }
)

setGeneric(name="writeOut",
           def=function(Object, outDir, original) {
             standardGeneric("writeOut")
           }
)

setGeneric(name="frag4feature",
           def=function(pa, xset, CSVfile, ppm=5, plim=NA, intense=TRUE, convert2RawRT=TRUE,
                        create_db=FALSE, out_dir='.', db_name=NA, grp_peaklist=NA){
             standardGeneric("frag4feature")
           }
)

setGeneric(name="averageIntraFragSpectra",
           def=function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                        av='median', sumi=TRUE, rmp=FALSE, cores=1) {
             standardGeneric("averageIntraFragSpectra")
           }
)

setGeneric(name="averageInterFragSpectra",
           def=function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                        av='median', sumi=TRUE, rmp=FALSE, cores=1) {
             standardGeneric("averageInterFragSpectra")
           }
)

setGeneric(name="averageAllFragSpectra",
           def=function(pa, minfrac=0.5, minnum=1, ppm=5, snr=0.0, ra=0.0,
                        av='median', sumi=TRUE,  rmp=FALSE, cores=1) {
             standardGeneric("averageAllFragSpectra")
           }
)


setGeneric(name="filterFragSpectra",
           def=function(pa, ilim=0, plim=0.8, ra=0, snr=3, cores=1, rmp=FALSE, snmeth='median') {
             standardGeneric("filterFragSpectra")
           }
)


setGeneric(name="createMSP",
           def=function(pa, msp_file_pth=NULL, metadata=NULL, metadata_cols=NULL,
                        xcms_groupids=NULL, method="all", adduct_split=TRUE, filter=TRUE,
                        msp_schema='massbank', intensity_ra='intensity_ra') {
             standardGeneric("createMSP")
           }
)


setGeneric(name="validate",
           def=function(pa, ppLCMS){
             standardGeneric("validate")
           }
)


