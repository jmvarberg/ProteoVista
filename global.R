suppressPackageStartupMessages({
    library(shiny)
    #library(msdap)
    library(data.table)
    library(waiter)
    library(shinyalert)
    library(shinyBS)
    library(shinyWidgets)
    library(rmarkdown)
    library(cowplot)
    library(plotly)
    library(ggsci)
    library(tidyverse)
})

#specify maximum file upload size (set at 500MB)
options(shiny.maxRequestSize=500*1024^2)

#HTML tag for data upload wait window
data_ingest_waiting <- tagList(
    spin_dots(),
    br(),
    br(),
    h4("Input data uploading, please wait...")
)

#HTML tag for data upload wait window
data_processing_waiting <- tagList(
    spin_dots(),
    br(),
    br(),
    h4("MS-DAP processing, please wait...")
)

#function to convert to contrast pairs to output filtered list pretty
print_list <- function(x) {
    cont <- paste0(x[1], " vs. ", x[2])
}

# # Function to filter out items that don't have reference group as the first object
# filter_ref <- function(lst, reference_group) {
#     lst[sapply(lst, function(x) x[1] == reference_group)] #convert to make this ref to the selected group.
# }

# Function to generate all pairs with a specified reference from a character vector
generate_contrasts_to_ref <- function(elements, reference) {

    # Check if the reference is in the elements vector
    if (!(reference %in% elements)) {
        stop("Reference element not found in the input vector.")
    }

    # Create pairs where each element is paired with the reference except itself
    pairs <- setdiff(elements, reference) # Exclude the reference from the elements
    result <- paste(reference, pairs, sep = "__") # Create pairs

    return(result)
}


#' Format character vector into HTML bulleted list - stolen from Kyle B. comment here: https://stackoverflow.com/questions/22923784/how-to-add-bullet-points-in-r-shinys-rendertext
#'
#' @param char a character vector. Each element will be a bullet
#' @param ordered logical (T/F). If `TRUE`, return numbered list.
#'
#' @keywords internal
#'
format_html_list <- function(char, ordered = FALSE){

    seps <- c("<li>", "</li>")
    html_wrapper <-  if(ordered) c("<ol>", "</ol>") else c("<ul>", "</ul>")

    bullets <- paste0(seps[1], char, seps[2], collapse = "")

    html_list <- paste0(html_wrapper[1], bullets, html_wrapper[2])

    return(html_list)
}
