test_that(
    "PurityXParam retains complete settings and validates their meaning",
    {
        defaults <- as.list(PurityXParam())
        expect_identical(names(defaults), c("purityType", "offsets", "iwNorm",
                "iwNormFun", "ilim", "isotopes", "im"))
        expect_null(defaults$iwNormFun)
        expect_null(defaults$im)
        expect_equal(defaults$offsets, c(0.5, 0.5))
        expect_s4_class(PurityXParam(), "Param")

        weighting <- function(offset) exp(-offset^2)
        parameter <- PurityXParam(iwNorm = TRUE, iwNormFun = weighting)
        path <- tempfile(fileext = ".rds")
        on.exit(unlink(path))
        saveRDS(parameter, path)
        expect_equal(as.list(readRDS(path)), as.list(parameter))
        expect_error(PurityXParam(offsets = c(0, 0)), "positive total width")
        expect_error(PurityXParam(ilim = 1.1), "between zero and one")
        expect_error(PurityXParam(purityType = "median"), "prediction mode")
        expect_error(PurityXParam(iwNormFun = 1), "function or NULL")
        expect_error(PurityXParam(im = matrix(1, 1, 2)), "eight columns")
    }
)

test_that(
    "prediction parameters replace defaults but reject explicit conflicts",
    {
        scalar <- as.list(PurityXParam())
        special <- PurityXParam(ilim = 0)
        expect_equal(.px_resolve_parameters(special, scalar)$ilim, 0)
        expect_error(.px_resolve_parameters(special, scalar, "ilim"),
            "Explicit ilim conflicts")
        scalar$ilim <- 0
        expect_equal(.px_resolve_parameters(special, scalar, "ilim"),
            as.list(special))
        expect_error(.px_resolve_parameters(list(default = special), scalar),
            "one PurityXParam")

        model <- matrix(seq_len(8L), nrow = 1L)
        parameter <- PurityXParam(im = model)
        scalar$im <- matrix(as.double(model), nrow = 1L)
        expect_identical(.px_resolve_parameters(parameter, scalar, "im")$im,
            model)
        dim(scalar$im) <- c(4L, 2L)
        expect_error(.px_resolve_parameters(parameter, scalar, "im"),
            "Explicit im conflicts")
    }
)

test_that(
    "prediction history keeps actual references and stable input fingerprints",
    {
        parameter <- PurityXParam()
        prepared <- list(
            inputInfo = list(class = "xcmsSet", source_files = "sample.mzML"),
            targets = data.frame(sample = 1L, c_peak_id = 3L, mz = 205.0972),
            groups = data.frame(grpid = c(2L, 3L),
                feature_id = c("feature-two", "feature-three")),
            exclusions = data.frame(), scanpeaks = list(matrix(c(1, 2), 1)))
        calculations <- list(
            predictions = data.frame(output_id = 2L, sample = 1L,
                c_peak_id = 3L, status = "available"),
            observations = data.frame(output_id = 2L, sample = 1L,
                scan = c(4L, 6L), msLevel = 1L),
            excluded_peaks = data.frame())
        object <- methods::new("purityX",
            predictions = data.frame(grpid = c(2L, 3L),
                median = c(0.8, NA_real_)))
        controls <- list(cores = 2L, effective_cores = 1L)
        first <- .px_record(object, prepared, calculations, as.list(parameter),
            controls, Sys.time())
        second <- .px_record(object, prepared, calculations, as.list(parameter),
            controls, Sys.time())
        record <- processingHistory(first, step = "purityX")[[1L]]
        expect_identical(record$input$observations, calculations$observations)
        expect_identical(record$output$predictions, calculations$predictions)
        expect_identical(record$output$prediction_rows$feature_id,
            c("feature-two", "feature-three"))
        expect_identical(record$parameters$resolved$isotope_model_source,
            "not_used")
        expect_false("scanpeaks" %in% names(record$input))
        expect_identical(record$input$fingerprint,
            processingHistory(second)[[1L]]$input$fingerprint)
        expect_false(identical(names(processingHistory(first)),
                names(processingHistory(second))))
        parameter@ilim <- 0
        expect_equal(processingParameters(first)[[1L]]$ilim, 0.05)
        expect_length(processingHistory(first, "purityA"), 0L)

        path <- tempfile(fileext = ".rds")
        on.exit(unlink(path))
        saveRDS(first, path)
        expect_identical(processingHistory(readRDS(path)),
            processingHistory(first))
        prepared$targets$mz <- 205.1000
        changed <- .px_record(object, prepared, calculations,
            as.list(parameter), controls, Sys.time())
        expect_false(identical(record$input$fingerprint,
                processingHistory(changed)[[1L]]$input$fingerprint))
    }
)

test_that(
    "saved purityX results upgrade without inventing execution settings",
    {
        path <- system.file("extdata", "tests", "purityX",
            "1_purityX_grouped_px.rds", package = "msPurity")
        historical <- readRDS(path)
        current <- updateObject(historical)
        expect_identical(current@predictions, historical@predictions)
        expect_true(methods::validObject(current))
        record <- processingHistory(current)[[1L]]
        expect_identical(record$status, "imported")
        expect_false(record$complete)
        expect_true(is.na(record$time$started))
        expect_true(is.na(record$software$version))
        expect_false("ilim" %in% names(record$parameters))
        expect_identical(processingHistory(historical),
            processingHistory(current))
        expect_identical(processingHistory(updateObject(current)),
            processingHistory(current))
    }
)
