suppressPackageStartupMessages({
    library(shiny)
    library(msdap)
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
    library(bslib)
})

#specify maximum file upload size (set at 500MB)
options(shiny.maxRequestSize=700*1024^2)

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


#' Generate a data.frame from list of vectors of uneven length
#'
#' @note For situations where you want to generate a data frame with columns that have different lengths.
#' This works by calculating the maximum length of objects in the list, then padding the remaining list objects
#' to fill with empty spaces so that a data frame can be created.
#'
#' @param list Named list of character vectors to be used to generate data frame, with each vector ending up as a unique column.
#'
#' @return Data frame object with each vector from the input list as a column.
#' @export
#'
#' @examples
#' List <- list(A = sample(letters, 15), B = sample(letters, 20), C = sample(letters, 5))
#' output <- jmv_mixedLengthDF(List)
#'
jmv_mixedLengthDF <- function(list) {

    #get maxlength of up clusters
    maxlength <- max(unlist(lapply(list, function(x) length(x))))

    #add "" to list objects to make equal length
    fill_vector <- function(x, max) {

        filled <- c(x, rep("", max - length(x)))

    }

    #get the max length in list
    list_filled <- lapply(list, fill_vector, max = maxlength)

    #make data frame from the list
    df <- do.call(cbind.data.frame, list_filled)

}


#summary stats: Total/Avg Peptides, Total/Avg Proteins, Avg DE proteins as valueBoxes
summary_stats_dashboard <- function(dataset) {

    #total/average number of peptides
    peptide_data <- dataset$peptides |>
        dplyr::filter(detect = TRUE)

    total_peptides <- length(unique(peptide_data$peptide_id))
    avg_peptides <- peptide_data |> dplyr::group_by(sample_id) |> dplyr::count() |> dplyr::ungroup() |> dplyr::summarise(across(n, mean, na.rm=TRUE)) |> dplyr::pull(n) |> round(digits = 0)

    #get protein information
    total_proteins <- length(unique(peptide_data$protein_id))
    avg_proteins <- peptide_data |> dplyr::group_by(sample_id, protein_id) |> dplyr::count() |> dplyr::ungroup() |> dplyr::select(-n) |>  dplyr::group_by(sample_id) |> dplyr::count() |> dplyr::ungroup() |> dplyr::summarise(across(n, mean, na.rm=TRUE)) |> dplyr::pull(n) |> round(digits = 0)

    #get avg number of DE proteins
    avg_de_proteins <- dataset$de_proteins |>
        dplyr::filter(signif==TRUE) |>
        dplyr::group_by(contrast) |>
        dplyr::count() |>
        dplyr::ungroup() |>
        dplyr::summarise(across(n, mean, na.rm=TRUE)) |>
        dplyr::pull(n) |>
        round(digits = 0)

    return(list(total_peptides = total_peptides, avg_peptides = avg_peptides, total_proteins = total_proteins, avg_proteins = avg_proteins, avg_de_proteins = avg_de_proteins))


}
