.qvq_args <- list(q_spectraTypes = "av_all", q_spectraFilter = TRUE,
                  q_pol = NA, l_spectraTypes = "av_all", l_pol = NA,
                  l_spectraFilter = TRUE, cores = 1)

.example <- function(...) system.file("extdata", "tests", ...,
                                      package = "msPurity")

# Tests append results and edit manifests. Convert once and give each test
# its own writable copy of the original generation.
.query_dataset <- local({
    template <- NULL
    function() {
        skip_if_no_arrow()
        if (is.null(template) || !dir.exists(template)) {
            path <- tempfile("q-")
            .without_identity_warning(convertSqliteToParquet(
                .example("db", "createDatabase_example.sqlite"), path))
            template <<- path
        }
        .copy_dataset(template)
    }
})

.sm <- function(q, l, ...) suppressMessages(do.call(
    spectralMatching, c(list(q_dbPth = q, l_dbPth = l, format = "parquet"),
                        utils::modifyList(.qvq_args, list(...)))))
