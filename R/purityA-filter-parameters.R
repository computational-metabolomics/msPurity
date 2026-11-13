# Acquisition metadata shared by processing and averaging summaries.

.PA_FILTER_AXES <- list(
    msLevel = c("msLevel", "ms_level"),
    massAnalyzer = c("massAnalyzer", "analyzer", "analyzerType",
        "massAnalyzerType", "analyzerName"),
    scanFilterAnalyzer = "scanFilterAnalyzer",
    activationType = c("activationType", "activationMethod"),
    collisionEnergy = c("collisionEnergy", "collision_energy"),
    collisionEnergyUnit = c("collisionEnergyUnit", "collisionEnergyUnits",
        "collisionEnergyUnitName"),
    polarity = "polarity"
)

.PA_ACTIVATION_TYPES <- c(
    CID = "CID", "collision-induced dissociation" = "CID",
    "MS:1000133" = "CID",
    HCD = "HCD", "beam-type collision-induced dissociation" = "HCD",
    "MS:1000422" = "HCD",
    ETD = "ETD", "electron transfer dissociation" = "ETD",
    "MS:1000598" = "ETD")

.pa_activation_type <- function(values) {
    values <- as.character(values)
    index <- match(tolower(trimws(values)),
        tolower(names(.PA_ACTIVATION_TYPES)))
    known <- !is.na(index)
    values[known] <- unname(.PA_ACTIVATION_TYPES[index[known]])
    values
}

.pa_filter_string_metadata <- function(strings, ms_level) {
    result <- data.frame(
        scanFilterAnalyzer = rep(NA_character_, length(strings)),
        activationType = rep(NA_character_, length(strings)))
    number <- "(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
    event <- paste0("^", number, "@([[:alpha:]]+)", number, "$")
    range <- paste0("\\[\\s*", number, "\\s*-\\s*", number,
        "(?:\\s*,\\s*", number, "\\s*-\\s*", number, ")*\\s*\\]")
    for (i in seq_along(strings)) {
        text <- strings[i]
        if (is.na(text)) next
        prefix <- regmatches(text, regexpr("^(FTMS|ITMS)\\b", text,
                ignore.case = TRUE, perl = TRUE))
        if (length(prefix)) result$scanFilterAnalyzer[i] <- toupper(prefix)
        marker <- regexpr("\\bms[0-9]*\\b", text, ignore.case = TRUE,
            perl = TRUE)
        if (marker < 0L || is.na(ms_level[i])) next
        level <- sub("^ms", "", tolower(regmatches(text, marker)))
        level <- if (nzchar(level)) as.numeric(level) else 1
        if (level < 2L || level != ms_level[i]) next
        tail <- substring(text, marker + attr(marker, "match.length"))
        tail <- trimws(sub(paste0("\\s*", range, "\\s*$"), "", tail,
                perl = TRUE))
        events <- strsplit(tail, "[[:space:]]+")[[1L]]
        if (length(events) != level - 1L ||
                !all(grepl(event, events, ignore.case = TRUE, perl = TRUE)))
            next
        method <- toupper(sub(event, "\\1", events[length(events)],
                ignore.case = TRUE, perl = TRUE))
        if (method %in% c("CID", "HCD"))
            result$activationType[i] <- method
    }
    result
}

.pa_filter_metadata <- function(source) {
    n <- if (is.data.frame(source)) nrow(source) else length(source)
    available <- if (is.data.frame(source)) names(source) else
        Spectra::spectraVariables(source)
    values <- lapply(names(.PA_FILTER_AXES), function(axis) {
        numeric_axis <- axis %in% c("msLevel", "collisionEnergy", "polarity")
        out <- if (numeric_axis) rep(NA_real_, n) else rep(NA_character_, n)
        for (name in intersect(.PA_FILTER_AXES[[axis]], available)) {
            value <- source[[name]]
            if (numeric_axis) {
                if (!is.numeric(value) && !all(is.na(value)))
                    stop("Acquisition metadata ", name, " must be numeric.")
                value <- as.numeric(value)
                known <- value[!is.na(value)]
                if (any(!is.finite(known)))
                    stop("Acquisition metadata ", name,
                        " must be finite or NA.")
                if (axis == "msLevel" && any(known < 1 |
                            known != trunc(known)))
                    stop("Acquisition msLevel must be positive integer or NA.")
                if (axis == "collisionEnergy" && any(known < 0))
                    stop("Acquisition collisionEnergy must be nonnegative.")
                if (axis == "polarity" && any(!known %in% c(0, 1)))
                    stop("Acquisition polarity must be 0, 1 or NA.")
            } else {
                value <- as.character(value)
                value[!is.na(value) & (!nzchar(value) |
                            tolower(value) == "unknown")] <- NA_character_
            }
            compared_out <- if (axis == "activationType")
                .pa_activation_type(out) else out
            compared_value <- if (axis == "activationType")
                .pa_activation_type(value) else value
            if (any(!is.na(out) & !is.na(value) &
                        compared_out != compared_value))
                stop("Conflicting known acquisition aliases for ", axis, ".")
            take <- is.na(out) & !is.na(value)
            out[take] <- value[take]
        }
        if (axis %in% c("msLevel", "polarity")) as.integer(out) else out
    })
    names(values) <- names(.PA_FILTER_AXES)
    values <- as.data.frame(values, stringsAsFactors = FALSE)
    structured_activation <- rep(FALSE, n)
    if ("activationTypeAccession" %in% available) {
        accession <- as.character(source$activationTypeAccession)
        structured_activation <- !is.na(accession) & nzchar(trimws(accession))
        mapped <- .pa_activation_type(accession)
        label <- .pa_activation_type(values$activationType)
        methods <- unique(unname(.PA_ACTIVATION_TYPES))
        conflict <- label %in% methods & mapped %in% methods & label != mapped
        if (any(conflict))
            stop("Conflicting activation label and activationTypeAccession ",
                "at spectrum ", which(conflict)[1L], ".")
        take <- is.na(values$activationType) & mapped %in% methods
        values$activationType[take] <- mapped[take]
    }
    if ("filterString" %in% available) {
        derived <- .pa_filter_string_metadata(as.character(source$filterString),
            values$msLevel)
        for (field in names(derived)) {
            take <- is.na(values[[field]]) & !is.na(derived[[field]])
            if (field == "activationType")
                take <- take & !structured_activation
            values[[field]][take] <- derived[[field]][take]
        }
    }
    values
}

