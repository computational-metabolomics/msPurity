with_ms_pth <- '../../Dropbox/msms_and_ms.mzML'
xset <- xcms::xcmsSet(with_ms_pth, nSlaves = 1)
xset <- xcms::group(xset)
xset <- xcms::retcor(xset)
xset <- xcms::group(xset)


with_just_msms_pth <- '../../Dropbox/just_msms.mzML'
pa  <- purityA(with_just_msms_pth, interpol = "linear")

f4f_pa <- frag4feature(pa, xset, use_group=T)
