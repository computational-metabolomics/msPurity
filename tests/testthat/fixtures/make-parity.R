# Regenerate the parity references for the move to Spectra containers.
#
# Run from the package root with the code whose output is the reference:
#
#   Rscript tests/testthat/fixtures/make-parity.R
#
# The references were first generated before any change to how msPurity
# reads or stores spectra (feature/113/refactor at 8a87160). Regenerating them accepts the
# current output as correct, so only do it for an intended change in results.
#
# One entry has changed since: two_features_msp$max_metadata, after the fix
# to createMSP(method = "max"), which used to pick a scan other than the one
# with the highest precursor intensity when a scan was linked to a feature
# more than once.

suppressPackageStartupMessages({
    pkgload::load_all(".", quiet = TRUE)
    library(xcms)
    library(testthat)
})
for (h in c("helper-sqlite-golden.R", "helper-parity.R"))
    source(file.path("tests", "testthat", h))

out <- file.path("tests", "testthat", "fixtures", "parity")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

message("purityA workflows")
saveRDS(.parity_purityA(), file.path(out, "purityA.rds"), compress = "xz")
message("purityD workflow")
saveRDS(.parity_purityD(), file.path(out, "purityD.rds"), compress = "xz")
