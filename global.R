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

#specify maximum file upload size (set at 1000MB)
options(shiny.maxRequestSize=1000*1024^2)

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

#function to add differential detection results to DEA results
add_dd_to_de <- function(de_results, dd_results) {

    #get min pval and adj. p val for dea results
    min_pval <- min(de_results$P.Value, na.rm=TRUE)
    print(min_pval)
    min_adjpval <- min(de_results$Adj.P.Value, na.rm=TRUE)
    print(min_adjpval)

    #replace zero's with min values from DE results to get properly displayed in Quickomics plots.
    dd_results <- dd_results |>
        dplyr::mutate(P.Value = dplyr::if_else(dea_algorithm=="Differential.Detection", min_pval, P.Value),
                      Adj.P.Value = dplyr::if_else(dea_algorithm=="Diferential.Detection", min_adjpval, Adj.P.Value))

    #combine the results by binding rows together.
    combined <- de_results |> dplyr::bind_rows(dd_results)
    return(combined)
}

#quickomics expression set subfolders
quickomics_expression_sets <- function(msdap_output_directory, msdap_data) {

    #msdap_output_directory <- "./ProteoVista_output/PROT-1027 Mouse Amyloid Fractionation_20240918140518/msdap_output/"

    # #create quickomics output folder in parent folder of msdap_output folder
    # dir.create(paste0(dirname(msdap_output_directory), "/quickomics_files/"))

    #load("./ProteoVista_output/PROT-1027 Mouse Amyloid Fractionation_20240918140518/msdap_output/2024-09-18_14-06-23/dataset.RData")
    #load("./ProteoVista_output/PROT-1027 Mouse Amyloid Fractionation_20240918140518/msdap_output/ExpressionSet_proteins_contrast wt_ev vs wt_mv.RData")
    dataset = msdap_data

    #need to load in the dd_proteins.tsv.gz results to filter and use later on.
    diff_detect_file <- list.files(msdap_output_directory, pattern = "dd_proteins.tsv.gz", recursive = TRUE, full.names=TRUE)
    diff_detect_results <- data.table::fread(diff_detect_file)

    #get the list of protein abundance data files for each contrast
    prot_abund_by_contrast <- list.files(path=msdap_output_directory, pattern = "protein_abundance__filter by contrast", recursive=TRUE, full.names = TRUE)

    #get the contrast names
    contrasts <- basename(stringr::str_remove_all(prot_abund_by_contrast, pattern = "protein_abundance__filter by contrast; ")) |>
        stringr::str_remove_all(pattern = ".tsv") |>
        stringr::str_replace_all(pattern = " ", "_")

    #for each contrast, create output subfolder in quickomics output folder
    lapply(contrasts, FUN = function(x) {dir.create(paste0(dirname(msdap_output_directory), "/quickomics_files/", x, "/"))})

    #now, process the files for each contrast and format/save output into the corresponding directory
    process_contrast_for_quickomics <- function(contrast_name, protein_abundance_results) {

        #contrast_name <- contrasts[1]
        #contrast_name = "contrast: bio_lysate vs bio_pellet"

        #Step 1: find the protein abundance data for the contrast and format for quickomics, save out
        prot_quan <- data.table::fread(prot_abund_by_contrast[stringr::str_detect(stringr::str_replace_all(prot_abund_by_contrast, pattern = " ", "_"), pattern = contrast_name)])

        quickomics_data <- prot_quan |>
            dplyr::select(-fasta_headers, -gene_symbols_or_id) |>
            dplyr::rename("UniqueID" = "protein_id")
        write.csv(quickomics_data, paste0(dirname(msdap_output_directory), "/quickomics_files/", contrast_name, "/quickomics_expression_data.csv"), row.names = FALSE)

        #Step 2: get the DEA results, filter for contrast, export for each algorithm used.
        de_results <- data.table::fread(list.files(path = msdap_output_directory, pattern = "de_proteins.tsv.gz", recursive=TRUE, full.names=TRUE))

        #select protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm. Filter for dea_algorithm of choice. Modify column names to match expected inputs for quickomics.
        #Required column names are UniqueId (protein group), test (comparison), Adj.P.Value, P.Value, and logFC
        quickomics_de <- de_results |>
            dplyr::select(protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm) |>
            dplyr::rename("UniqueID" = "protein_id",
                          "test" = "contrast",
                          "Adj.P.Value" = "qvalue",
                          "P.Value" = "pvalue",
                          "logFC" = "foldchange.log2") |>
            dplyr::mutate(test = stringr::str_remove_all(test, pattern = "contrast: "),
                          test = stringr::str_replace(test, pattern = " vs ", "-"),
                          logFC = -logFC) |>  #this reverses the standard notation for FC that MS-DAP uses (FC = A/B instead of conventional FC = B/A)
            dplyr::filter(stringr::str_detect(string = test, pattern = stringr::str_replace(contrast_name, pattern = "_vs_", replacement = "-")))

        #Now, split out by test and save csv's for each test method
        quickomics_de_list <- quickomics_de |> dplyr::group_split(dea_algorithm)
        de_names <- lapply(quickomics_de_list, function(x) unique(x$dea_algorithm))
        names(quickomics_de_list) <- as.character(unlist(de_names))

        #before saving out, bind rows for differentially detected proteins for this contrast to all de outputs.
        #step 1: filter the dd results for the contrast of interest and only retaining the differentiall detected with |z-score| >= 4, select necessary columns and rename to match the DE results columns
        dd_contrast <- diff_detect_results |>
            dplyr::select(protein_id, log2fc, contrast, zscore) |>
            dplyr::mutate(contrast = stringr::str_remove_all(contrast, pattern = "contrast: "),
                          contrast = stringr::str_replace_all(contrast, pattern = " ", "_"),
                          log2fc = -log2fc,
                          Adj.P.Value = 0,
                          P.Value = 0,
                          dea_algorithm = "Differential.Detection") |>
            dplyr::filter(abs(zscore) >= 4, contrast == contrast_name) |>
            dplyr::select(-zscore) |>
            dplyr::mutate(contrast = stringr::str_replace(contrast, pattern = "_vs_", "-")) |>
            dplyr::rename("UniqueID" = "protein_id",
                          "test" = "contrast",
                          "logFC" = "log2fc") |>
            dplyr::select(UniqueID, test, Adj.P.Value, P.Value, logFC, dea_algorithm)

        quickomics_de_list <- lapply(quickomics_de_list, add_dd_to_de, dd_results = dd_contrast)

        #save out each csv file
        paths <- lapply(de_names, function(x) paste0(dirname(msdap_output_directory), "/quickomics_files/", contrast_name, "/", x, "_dea_results.csv")) |> unlist()
        purrr::map2(quickomics_de_list, paths, .f = function(x, y) write.csv(x, file = y, row.names = F))

        #Step 3: Get the Sample Metadata Table - load the sample table then filter for the samples in the contrast

        #sample_table <- data.table::fread(list.files(msdap_output_directory, pattern = "samples.tsv.gz", recursive=TRUE, full.names=TRUE))
        #msdap_output_directory <- "./ProteoVista_output/PROT-1007 Smed X1 Subsort Populations OnePot Actin vs Tubulin_20240923165754/msdap_output/"
        eset_protein_files <- list.files(msdap_output_directory, pattern = "ExpressionSet_proteins_contrast", recursive=TRUE, full.names=TRUE)
        contrast_eset <- eset_protein_files[stringr::str_detect(stringr::str_replace_all(eset_protein_files, pattern = " ", replacement = "_"), pattern = contrast_name)]
        #load(eset_protein_files[1])
        load(contrast_eset)
        eset_proteins

        #Make sure that the values match the column names in the protein expression data
        sample_table <- eset_proteins@phenoData@data

        #Filter to remove any excluded samples
        sample_table_filtered <- sample_table |> dplyr::filter(!exclude)

        columns_present <- colnames(quickomics_data)
        test <- length(intersect(sample_table_filtered$sample_id, columns_present)) == length(unique(sample_table_filtered$sample_id))

        #if not true, then throw an error and stop
        if(!test) {
            stop("Columns in sample metadata not matching to names in the expression data. Check dataset$samples$sample_id values and compare with values in column names for protein abundance global filter dataset.")
        } else (
            print("All samples in eset object are present in the expression data. Proceeding with extraction.")
        )

        #Modify columns to match expected for quickomics upload
        #Sample MetaData file. Must have columns "sampleid" and "group", "Order" and "ComparePairs" columns
        quickomics_md <- jmv_mixedLengthDF(list(sampleid = sample_table_filtered$sample_id,
                                                group = sample_table_filtered$group,
                                                Order = unique(sample_table_filtered$group),
                                                ComparePairs = unique(quickomics_de$test)))
        quickomics_md
        print('Now saving metadata file: ')
        print(paste0(dirname(msdap_output_directory), "/quickomics_files/", contrast_name, "/quickomics_sample_metadata.csv"))

        write.csv(quickomics_md, paste0(dirname(msdap_output_directory), "/quickomics_files/", contrast_name, "/quickomics_sample_metadata.csv"), row.names=FALSE)

        #Step 4: Gene name table
        gene_table <- dataset$proteins

        #Gene/Protein Name File - must have four columns: id (sequential numbers), UniqueID (must match the rownames in expression data and comparison data), Gene.Name, Protein.ID. Can have additional

         quickomics_gene_table <- gene_table |>
             dplyr::mutate(id = dplyr::row_number()) |>
             dplyr::rename("UniqueID" = "protein_id",
                           "Gene.Name" = "gene_symbols_or_id",
                           "Protein.ID" = "accessions")
        write.csv(quickomics_gene_table, paste0(dirname(msdap_output_directory), "/quickomics_files/", contrast_name, "/quickomics_gene_protein_table.csv"), row.names = F)


    }

    #now run the processing steps
    #lapply(contrasts, FUN = process_contrast_for_quickomics, protein_abundance_results = prot_abund_by_contrast) #COMMENTED OUT 10/7 due to issue with FBC expression values not matching global filtered.

    #Now, repeat above but for global proteins for all contrasts.
    #Modified 9/18/24 to use the protein abundances filtered by groups instead of globally.
    process_global_for_quickomics <- function(msdap_output_directory) {

        #create global quickomics output folder
        global_output <- paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/")
        dir.create(global_output)

        #msdap_output_directory <- "./ProteoVista_output/PROT-1007 Smed X1 Subsort OnePot Tubulin and Actin, Non-redundant Smed DB 2024 Eric Ross_20240918085352/msdap_output/"

        #Step 1: load the global proteomics results and format/save out
        prot_quan_file <- list.files(msdap_output_directory, pattern = "protein_abundance__filter by group independently.tsv", recursive = TRUE, full.names=TRUE)
        #prot_quan_file <- data.table::fread("./ProteoVista_output/PROT-1007 Smed X1 Subsort Cytometry OnePot, Tubulin and Actin_20240917144856/msdap_output/2024-09-17_14-52-44/protein_abundance__global data filter.tsv")
        prot_quan <- data.table::fread(prot_quan_file)

        quickomics_data <- prot_quan |>
            dplyr::select(-fasta_headers, -gene_symbols_or_id) |>
            dplyr::rename("UniqueID" = "protein_id")
        write.csv(quickomics_data, paste0(global_output, "/global_quickomics_expression_data.csv"), row.names = FALSE)

        #Step 2:
        #Step 2: get the DEA results, filter for contrast, export for each algorithm used.
        de_results <- data.table::fread(list.files(path = msdap_output_directory, pattern = "de_proteins.tsv.gz", recursive=TRUE, full.names=TRUE))

        #select protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm. Filter for dea_algorithm of choice. Modify column names to match expected inputs for quickomics.
        #Required column names are UniqueId (protein group), test (comparison), Adj.P.Value, P.Value, and logFC
        quickomics_de <- de_results |>
            dplyr::select(protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm) |>
            dplyr::rename("UniqueID" = "protein_id",
                          "test" = "contrast",
                          "Adj.P.Value" = "qvalue",
                          "P.Value" = "pvalue",
                          "logFC" = "foldchange.log2") |>
            dplyr::mutate(test = stringr::str_remove_all(test, pattern = "contrast: "),
                          test = stringr::str_replace(test, pattern = " vs ", "-"),
                          logFC = -logFC)   #this reverses the standard notation for FC that MS-DAP uses (FC = A/B instead of conventional FC = B/A)


        #Now, split out by test and save csv's for each test method
        quickomics_de_list <- quickomics_de |> dplyr::group_split(dea_algorithm)
        de_names <- lapply(quickomics_de_list, function(x) unique(x$dea_algorithm))
        names(quickomics_de_list) <- as.character(unlist(de_names))

        #before saving out, bind rows for differentially detected proteins for this contrast to all de outputs.
        #step 1: filter the dd results and only retain the differential detected with |z-score| >= 4, select necessary columns and rename to match the DE results columns
        dd_contrast <- diff_detect_results |>
            dplyr::select(protein_id, log2fc, contrast, zscore) |>
            dplyr::mutate(contrast = stringr::str_remove_all(contrast, pattern = "contrast: "),
                          contrast = stringr::str_replace_all(contrast, pattern = " ", "_"),
                          log2fc = -log2fc,
                          Adj.P.Value = 0,
                          P.Value = 0,
                          dea_algorithm = "Differential.Detection") |>
            dplyr::filter(abs(zscore) >= 4) |>
            dplyr::select(-zscore) |>
            dplyr::mutate(contrast = stringr::str_replace(contrast, pattern = "_vs_", "-")) |>
            dplyr::rename("UniqueID" = "protein_id",
                          "test" = "contrast",
                          "logFC" = "log2fc") |>
            dplyr::select(UniqueID, test, Adj.P.Value, P.Value, logFC, dea_algorithm)

        quickomics_de_list <- lapply(quickomics_de_list, add_dd_to_de, dd_results = dd_contrast)

        #save out each csv file
        paths <- lapply(de_names, function(x) paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/", x, "_global_dea_results.csv")) |> unlist()
        purrr::map2(quickomics_de_list, paths, .f = function(x, y) write.csv(x, file = y, row.names = F))

        #Step 3: Get the Sample Metadata Table - load the sample table then filter for the samples in the contrast

        sample_table <- data.table::fread(list.files(msdap_output_directory, pattern = "samples.tsv.gz", recursive=TRUE, full.names=TRUE))

        #remove any samples that were excluded
        #samples <- data.table::fread("./ProteoVista_output/PROT-1007 Smed X1 Subsort Cytometry OnePot, Tubulin and Actin_20240917144856/msdap_output/2024-09-17_14-52-44/samples.tsv.gz")

        sample_table_filtered <- sample_table |> dplyr::filter(!exclude)

        columns_present <- colnames(quickomics_data)
        test <- length(intersect(sample_table_filtered$sample_id, columns_present)) == length(unique(sample_table_filtered$sample_id))

        #if not true, then throw an error and stop
        if(!test) {
            stop("Columns in sample metadata not matching to names in the expression data. Check dataset$samples$sample_id values and compare with values in column names for protein abundance global filter dataset.")
        } else (
            print("All samples in are present in the expression data. Proceeding with extraction.")
        )

        #Modify columns to match expected for quickomics upload
        #Sample MetaData file. Must have columns "sampleid" and "group", "Order" and "ComparePairs" columns
        quickomics_md <- jmv_mixedLengthDF(list(sampleid = sample_table_filtered$sample_id,
                                                group = sample_table_filtered$group,
                                                Order = unique(sample_table_filtered$group),
                                                ComparePairs = unique(quickomics_de$test)))
        quickomics_md
        print('Now saving metadata file: ')
        print(paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/global_quickomics_sample_metadata.csv"))

        write.csv(quickomics_md, paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/global_quickomics_sample_metadata.csv"), row.names=FALSE)

        #Step 4: Gene name table
        gene_table <- dataset$proteins

        #Gene/Protein Name File - must have four columns: id (sequential numbers), UniqueID (must match the rownames in expression data and comparison data), Gene.Name, Protein.ID. Can have additional

        quickomics_gene_table <- gene_table |>
            dplyr::mutate(id = dplyr::row_number()) |>
            dplyr::rename("UniqueID" = "protein_id",
                          "Gene.Name" = "gene_symbols_or_id",
                          "Protein.ID" = "accessions")
        write.csv(quickomics_gene_table, paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/global_quickomics_gene_protein_table.csv"), row.names = F)


    }

    process_global_for_quickomics(msdap_output_directory)


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

