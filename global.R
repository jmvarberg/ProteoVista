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
    library(openxlsx)
    library(bslib)
    library(yaml)
    library(sessioninfo)
    library(shinymeta)
})

#specify maximum file upload size (set at 2000MB)
options(shiny.maxRequestSize=2000*1024^2)

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

#quickomics expression set subfolders
quickomics_expression_sets <- function(msdap_output_directory, msdap_data) {

    #DEBUGGING 12/16/25
    #msdap_output_directory <- "./ProteoVista_output/SIMR_ShortCourse_Jan2025_All_Combined_20250114211908/msdap_output/"
    #load("./ProteoVista_output/SIMR_ShortCourse_Jan2025_All_Combined_20250114211908/msdap_output/2025-01-14_21-24-59/dataset.RData")
    # #create quickomics output folder in parent folder of msdap_output folder
    # dir.create(paste0(dirname(msdap_output_directory), "/quickomics_files/"))

    #load("./ProteoVista_output/PROT-1027 Mouse Amyloid Fractionation_20240918140518/msdap_output/2024-09-18_14-06-23/dataset.RData")
    #load("./ProteoVista_output/PROT-1027 Mouse Amyloid Fractionation_20240918140518/msdap_output/ExpressionSet_proteins_contrast wt_ev vs wt_mv.RData")
    dataset = msdap_data

    #need to load in the dd_proteins.tsv.gz results to filter and use later on.

    #Need to modify to only do this if DD was performed. Fix after update to UI and specification of that input parameter
    diff_detect_file <- list.files(msdap_output_directory, pattern = "dd_proteins.tsv.gz", recursive = TRUE, full.names=TRUE)
    diff_detect_results <- data.table::fread(diff_detect_file)

    #Now, repeat above but for global proteins for all contrasts.
    #Modified 9/18/24 to use the protein abundances filtered by groups instead of globally.
    process_for_quickomics <- function(msdap_output_directory) {

        #create global quickomics output folder
        quickomics_output <- paste0(dirname(msdap_output_directory), "/quickomics_files/")
        dir.create(quickomics_output)

        #msdap_output_directory <- "./ProteoVista_output/PROT-1007 Smed X1 Subsort OnePot Tubulin and Actin, Non-redundant Smed DB 2024 Eric Ross_20240918085352/msdap_output/"

        #Step 1: load the MS-DAP FBG proteomics results and reformat/save out
        prot_quan_file <- list.files(msdap_output_directory, pattern = "protein_abundance__filter by group independently.tsv", recursive = TRUE, full.names=TRUE)
        #prot_quan_file <- data.table::fread("./ProteoVista_output/PROT-1007 Smed X1 Subsort Cytometry OnePot, Tubulin and Actin_20240917144856/msdap_output/2024-09-17_14-52-44/protein_abundance__global data filter.tsv")
        prot_quan <- data.table::fread(prot_quan_file)

        quickomics_data <- prot_quan |>
            dplyr::select(-fasta_headers, -gene_symbols_or_id) |>
            dplyr::rename("UniqueID" = "protein_id")
        write.csv(quickomics_data, paste0(quickomics_output, "/global_quickomics_expression_data.csv"), row.names = FALSE)

        #Step 2:
        #Step 2: get the DEA results export for each DE algorithm used.
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
            dplyr::mutate(test = stringr::str_remove_all(test, pattern = paste(c(" # condition.*", "contrast: "), collapse = "|")),
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
            dplyr::mutate(contrast = stringr::str_remove_all(contrast, pattern = paste(c(" # condition.*", "contrast: "), collapse = "|")),
                          log2fc = -log2fc,
                          Adj.P.Value = 0,
                          P.Value = 0,
                          dea_algorithm = "Differential.Detection") |>
            dplyr::filter(abs(zscore) >= 4) |>
            dplyr::select(-zscore) |>
            dplyr::mutate(contrast = stringr::str_replace(contrast, pattern = " vs ", "-")) |>
            dplyr::rename("UniqueID" = "protein_id",
                          "test" = "contrast",
                          "logFC" = "log2fc") |>
            dplyr::select(UniqueID, test, Adj.P.Value, P.Value, logFC, dea_algorithm)

        quickomics_de_list <- lapply(quickomics_de_list, add_dd_to_de, dd_results = dd_contrast)

        #save out each csv file
        paths <- lapply(de_names, function(x) paste0(dirname(msdap_output_directory), "/quickomics_files/", x, "_dea_results.csv")) |> unlist()
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
        print(paste0(dirname(msdap_output_directory), "/quickomics_files/quickomics_sample_metadata.csv"))

        write.csv(quickomics_md, paste0(dirname(msdap_output_directory), "/quickomics_files/quickomics_sample_metadata.csv"), row.names=FALSE)

        #Step 4: Gene name table
        gene_table <- dataset$proteins

        #Gene/Protein Name File - must have four columns: id (sequential numbers), UniqueID (must match the rownames in expression data and comparison data), Gene.Name, Protein.ID. Can have additional

        quickomics_gene_table <- gene_table |>
            dplyr::mutate(id = dplyr::row_number()) |>
            dplyr::rename("UniqueID" = "protein_id",
                          "Gene.Name" = "gene_symbols_or_id",
                          "Protein.ID" = "accessions")
        write.csv(quickomics_gene_table, paste0(dirname(msdap_output_directory), "/quickomics_files/quickomics_gene_protein_table.csv"), row.names = F)


    }

    process_global_for_quickomics(msdap_output_directory)


}

#summary stats: Total/Avg Peptides, Total/Avg Proteins, Avg DE proteins as valueBoxes
# summary_stats_dashboard <- function(dataset) {
#
#     #total/average number of peptides
#     peptide_data <- dataset$peptides |>
#         dplyr::filter(detect = TRUE)
#
#     total_peptides <- length(unique(peptide_data$peptide_id))
#     avg_peptides <- peptide_data |> dplyr::group_by(sample_id) |> dplyr::count() |> dplyr::ungroup() |> dplyr::summarise(across(n, mean, na.rm=TRUE)) |> dplyr::pull(n) |> round(digits = 0)
#
#     #get protein information
#     total_proteins <- length(unique(peptide_data$protein_id))
#     avg_proteins <- peptide_data |> dplyr::group_by(sample_id, protein_id) |> dplyr::count() |> dplyr::ungroup() |> dplyr::select(-n) |>  dplyr::group_by(sample_id) |> dplyr::count() |> dplyr::ungroup() |> dplyr::summarise(across(n, mean, na.rm=TRUE)) |> dplyr::pull(n) |> round(digits = 0)
#
#     #get avg number of DE proteins
#     avg_de_proteins <- dataset$de_proteins |>
#         dplyr::filter(signif==TRUE) |>
#         dplyr::group_by(contrast) |>
#         dplyr::count() |>
#         dplyr::ungroup() |>
#         dplyr::summarise(across(n, mean, na.rm=TRUE)) |>
#         dplyr::pull(n) |>
#         round(digits = 0)
#
#     return(list(total_peptides = total_peptides, avg_peptides = avg_peptides, total_proteins = total_proteins, avg_proteins = avg_proteins, avg_de_proteins = avg_de_proteins))
#
#
# }

#SIMR Excel output
simr_excel_report <- function(msdap_output_directory, msdap_dataset) {

    # #testing
    #load("./ProteoVista_output/PROT-1090 Smed Hoechst Neoblasts OnePot_20250307115425/msdap_output/dataset.qs")
    #dataset <- qs::qread("./ProteoVista_output/PROT-1090 Smed Hoechst Neoblasts OnePot_20250307115425/msdap_output/dataset.qs")
    #msdap_dataset <- dataset

    #read in the MSDAP report TSV directly to troubleshoot
    #msdap_report <- data.table::fread("./ProteoVista_output/PROT-1090 Smed Hoechst Neoblasts OnePot_20250307115425/input_data/PROT-1090_Smed_hoechst_neoblasts_g1_s_fed_unfed_with_bovine_contaminants_sme_Report_MS-DAP Format (Normal).tsv")

    #testing the yaml import
    #msdap_output_directory = "./ProteoVista_output/PROT-1090 Smed Hoechst Neoblasts OnePot_20250307115425/msdap_output/"
    #find, import and format yaml params
    yaml_file <- list.files(msdap_output_directory, pattern = ".yaml", full.names = T, recursive=T)

    if(file.exists(yaml_file)) {

        yaml_list <- yaml::yaml.load_file(yaml_file)

        yaml_df <- data.frame(
            Key = names(unlist(yaml_list, recursive = TRUE)),
            Value = unlist(yaml_list, recursive = TRUE),
            stringsAsFactors = FALSE
        )
        print(yaml_df)

    }

    #get sample md
    md <- msdap_dataset$samples

    md_out <- md |>
        dplyr::select(sample_index, sample_id, shortname, exclude, group)

    #get peptide info
    peptides <- msdap_dataset$peptides

    #get proteins info
    proteins <- msdap_dataset$proteins

    #create summary of number of peptides that pass detection, are not decoys, and are not contaminants per sample, and passed the by_group filtering criteria used
    pep_per_sample <- peptides |>
        dplyr::filter(detect, !stringr::str_detect(protein_id, "Cont_"), !isdecoy) |>
        dplyr::select(peptide_id, sample_id, intensity_by_group) |>
        na.omit() |>
        dplyr::group_by(sample_id) |>
        dplyr::count(name = "Num.Peptides") |>
        dplyr::left_join(select(md, sample_id, group))

    pep_per_group_summary <- pep_per_sample |>
        dplyr::group_by(group) |>
        rstatix::get_summary_stats(Num.Peptides, type = "common")

    #Now, need to get number of proteins. need to do the roll-up using the filter by group pep intensities.
    pep_data <- peptides |>
        dplyr::select(peptide_id, protein_id, sample_id, intensity_by_group) |>
        dplyr::rename(intensity = intensity_by_group)

    prot_data <- msdap::rollup_pep2prot_maxlfq(pep_data, intensity_is_log2 = T, implementation = "iq", return_as_matrix=F)

    #get wide matrix format with individual protein * sample expression values
    prot_data_wide <- msdap::rollup_pep2prot_maxlfq(pep_data, intensity_is_log2 = T, implementation = "iq", return_as_matrix=T)

    #map back the additional protein information, i.e., fasta header, gene symbols etc.
    prot_data_wide_out <- as.data.frame(prot_data_wide) |>
        tibble::rownames_to_column(var = "protein_id") |>
        dplyr::left_join(proteins) |>
        dplyr::select(-accessions) |>
        dplyr::select(protein_id, gene_symbols_or_id, fasta_headers, everything())

    #create summary of number of proteins per sample/group
    prot_per_sample <- prot_data |>
        dplyr::group_by(sample_id) |>
        dplyr::count(name = "Num.ProtGroups") |>
        dplyr::left_join(select(md, sample_id, group))

    pg_per_group_summary <- prot_per_sample |>
        dplyr::group_by(group) |>
        rstatix::get_summary_stats(Num.ProtGroups, type = "common")

    #now, get protein mean and SD expression value within group
    prot_group_exp_summary <- prot_data |>
        dplyr::left_join(select(md, sample_id, group)) |>
        dplyr::group_by(protein_id, group) |>
        dplyr::summarise(meanInt = mean(intensity, na.rm=T),
                         sdInt = sd(intensity, na.rm=T)) |>
        tidyr::pivot_wider(id_cols = c("protein_id"),
                           names_from = group,
                           names_glue = "{group}.{.value}",
                           values_from = c(meanInt, sdInt)) |>
        dplyr::ungroup() |>
        dplyr::left_join(proteins) |>
        dplyr::select(-accessions) |>
        dplyr::select(protein_id, gene_symbols_or_id, fasta_headers, everything())

    #get differential expression/detection results
    #To add: Filter this to a specific q-value for output in Shiny UI? for now, default to 0.05 for testing
    de <- msdap_dataset$de_proteins |>
        dplyr::filter(qvalue <= 0.05) |>
        dplyr::select(protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm, peptides_used_for_dea) |>
        dplyr::rename("UniqueID" = "protein_id",
                      "test" = "contrast",
                      "Adj.P.Value" = "qvalue",
                      "P.Value" = "pvalue",
                      "logFC" = "foldchange.log2") |>
        dplyr::mutate(test = stringr::str_remove_all(test, pattern = paste(c(" # condition.*", "contrast: "), collapse = "|")),
                      test = stringr::str_replace(test, pattern = " vs ", "-"),
                      logFC = -logFC) |>   #this reverses the standard notation for FC that MS-DAP uses (FC = A/B instead of conventional FC = B/A)
        dplyr::left_join(proteins, by = c("UniqueID" = "protein_id")) |>
        dplyr::select(-accessions) |>
        dplyr::select(UniqueID, gene_symbols_or_id, fasta_headers, everything())

    #now, want to make this a list of separate dfs for each DE algorithm used. Will then flatten this so that each one gets its own sheet in the Excel report
    de_list <- de |> dplyr::group_split(dea_algorithm)
    de_names <- lapply(de_list, function(x) unique(x$dea_algorithm))
    names(de_list) <- as.character(unlist(de_names))

    # Function to expand only lists that contain data frames
    expand_list <- function(lst, parent_name = "") {
        expanded_list <- list()

        for (name in names(lst)) {
            item <- lst[[name]]

            if (is.data.frame(item)) {
                # Keep data frames as they are
                expanded_list[[name]] <- item
            } else if (is.list(item)) {
                # Recursively expand sublists with parent name prefix
                sub_expanded <- expand_list(item, paste0(parent_name, name, "_"))
                expanded_list <- c(expanded_list, sub_expanded)
            }
        }

        return(expanded_list)
    }

    #To add: Filter this to a specific z-score for output in Shiny UI? for now, default to +/- 4 for testing
    dd <- msdap_dataset$dd_proteins |>
        dplyr::filter(abs(zscore) >= 4) |>
        dplyr::rename("UniqueID" = "protein_id",
                      "test" = "contrast") |>
        dplyr::mutate(test = stringr::str_remove_all(test, pattern = paste(c(" # condition.*", "contrast: "), collapse = "|")),
                      test = stringr::str_replace(test, pattern = " vs ", "-"),
                      log2fc = -log2fc) |> #this reverses the standard notation for FC that MS-DAP uses (FC = A/B instead of conventional FC = B/A)
        dplyr::left_join(proteins, by = c("UniqueID" = "protein_id")) |>
        dplyr::select(-accessions) |>
        dplyr::select(UniqueID, gene_symbols_or_id, fasta_headers, everything())

    #summarize number of up/down genes per contrast, per dea algorithm
    de_summary <- de |>
        dplyr::mutate(direction = dplyr::if_else(logFC > 0, "Up", "Down")) |>
        dplyr::group_by(test, dea_algorithm, direction) |>
        dplyr::count()

    dd_summary <- dd |>
        dplyr::mutate(direction = dplyr::if_else(log2fc > 0, "Up", "Down"),
                      dea_algorithm = "Differential.Detection") |>
        dplyr::group_by(test, dea_algorithm, direction) |>
        dplyr::count()

    #combine output
    de_dd_comb_summary <- de_summary |>
        dplyr::bind_rows(dd_summary)

    #create object for saving out as excel file
    excel_output <- list(
        params = yaml_df,
        metadata = md_out,
        pep_per_sample = pep_per_sample,
        pep_per_group = pep_per_group_summary,
        prot_per_sample = prot_per_sample,
        prot_per_group = pg_per_group_summary,
        prot_expression = prot_data_wide_out,
        prot_expression_group_summary = prot_group_exp_summary,
        differential_expression = de_list,
        differetial_detection = dd,
        differential_analysis_summary = de_dd_comb_summary
    )

    # Expand nested lists into separate sheets
    excel_output <- expand_list(excel_output)

    #save excel file to output directory
    print('Now saving SIMR Excel report file: ')
    openxlsx::write.xlsx(excel_output, paste0(dirname(msdap_output_directory), "/ProteoVista_summary_report.xlsx"))

    #return the list object
    return(excel_output)
}

#function to generate PAR - built from the generate_pdf_report() function in MS-DAP, modified to suit our needs.

generate_par_report <- function(dataset, output_dir) {

    ################ render Rmarkdown ################

    f = "./simr_report.Rmd"
    if(!file.exists(f)) {
        append_log(paste("cannot find report template file:", f), type = "error")
    }

    append_log("report: rendering report (this may take a while depending on dataset size)", type = "progress")

    ### prepare report files for RMarkdown rendering
    # rmarkdown known bugs: do not use output_dir  @  https://github.com/rstudio/rmarkdown/issues/861
    # a recommended work-around is to set the base/root directories inside the rmarkdown report using knitr::opts_chunk$set(...)
    # ref; https://github.com/yihui/knitr/issues/277#issuecomment-6528846
    # our robust work-around: since the render function uses the report document and it's dir as a base, we simply copy the Rmarkdown file to output dir and run render() with all default settings
    output_dir__temp = paste0(output_dir, "/temp", floor(as.numeric(Sys.time())) ) # to ensure unique dirname, simply add unix timestamp with seconds as precision
    dir.create(output_dir__temp, showWarnings = F) # don't care if directory already exists
    if(!dir.exists(output_dir__temp)) {
        append_log(paste("failed to create temp directory at;", output_dir__temp), type = "error")
    }

    # copy .Rmd file into temp directory, nested in chosen output dir
    f_newlocation = paste0(output_dir__temp, "/", basename(f))
    if(!file.copy(from = f, to = f_newlocation)) {
        append_log(paste("failed to copy report template from", f, "into to the temp directory:", f_newlocation), type = "error")
    }

    ### create the actual report
    rmarkdown::render(input = f_newlocation, output_file = "proteovista_analysis_report.docx", quiet = T, clean = T)

    # sanity check; was a report PDF created ?
    fpdf_templocation = paste0(output_dir__temp, "/proteovista_analysis_report.docx")
    if(!file.exists(fpdf_templocation)) {
        append_log(paste("failed to create ProteoVista report at:", fpdf_templocation), type = "error")
    }

    ### move report to output dir and remove temp dir
    fpdf_finallocation = paste0(output_dir, "/proteovista_analysis_report.docx")
    file.rename(fpdf_templocation, fpdf_finallocation)
    if(!file.exists(fpdf_finallocation)) {
        append_log(paste("failed to move the ProteoVista report from", fpdf_templocation, "to", fpdf_finallocation), type = "error")
    }

    # try to remove entire temp dir; may fail if user opened one of the files or is inside the dir in windows explorer
    # should be safe because we use a unique name in a dir we created previously AND we checked that this is an existing path where we have write access (otherwise above code would have failed)
    unlink(output_dir__temp, recursive = T, force = T) # use recursive=T, because unlink() documentation states: "If recursive = FALSE directories are not deleted, not even empty ones"

    append_log_timestamp("report:", start_time)
}

#testing capture of R commands to log file for non-shiny reprocessing
log_eval_with_args <- function(expr, log_file = "run_log.R") {
    # Evaluate first to capture values
    result <- eval(substitute(expr), envir = parent.frame())

    # Deparse the evaluated expression with evaluated arguments
    expr_str <- paste0("# Command issued at ", Sys.time(), "\n",
                       capture.output(dput(result)), "\n")

    writeLines(expr_str, con = log_file, sep = "\n", useBytes = TRUE, append = TRUE)
    return(result)
}
