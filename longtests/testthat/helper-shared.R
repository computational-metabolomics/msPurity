helpers <- system.file("test-helpers", package = "msPurity", mustWork = TRUE)
for (helper in list.files(helpers, "[.]R$", full.names = TRUE))
    sys.source(helper, envir = environment())
