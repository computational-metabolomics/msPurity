with_just_msms_pth <- '../../Dropbox/julein_france/just_msms.mzML'
with_ms_pth <- '../../Dropbox/julein_france/msms_and_ms.mzML'

# perform xcms on the data with that has MS1 data
xset <- xcms::xcmsSet(with_ms_pth, nSlaves = 1)
xset <- xcms::group(xset)

# process any files with purityA that have ms/ms
pa  <- purityA(c(with_just_msms_pth, with_ms_pth), interpol = "linear")
pa_linked <- frag4feature(pa, xset, use_group=T, create_db = T)  # make sure use_group set to TRUE
