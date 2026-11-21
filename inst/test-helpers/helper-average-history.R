# Original measured occurrences selected through each output's input graph.
.average_selected_links <- function(averages) {
    lapply(averages$spectrumId, function(id)
        averageMembers(averages, id, depth = "original"))
}

# Immediate spectrum IDs of each output.
.average_input_ids <- function(averages) {
    lapply(averages$spectrumId, function(id)
        averageMembers(averages, id, depth = "immediate")$member)
}
