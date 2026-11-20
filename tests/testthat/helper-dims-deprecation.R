# Compatibility regressions still exercise the retired DIMS interfaces.
# Only their retirement warnings are expected; other warnings remain visible.
.with_dims_deprecation <- function(expr) {
    withCallingHandlers(expr, deprecatedWarning = function(w) {
        if (identical(w$package, "msPurity") &&
            grepl("is deprecated as part of the DIMS workflow retirement.",
                  conditionMessage(w), fixed = TRUE))
            invokeRestart("muffleWarning")
    })
}
