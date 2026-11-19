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
