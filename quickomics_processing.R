#Functions to process MS-DAP output directory and generate files for Quickomics upload

#For testing
# msdap_output_directory = "./PROT-975_Smed_EV_MV_OnePot_ProteoVista_for_publication_20241028140925/msdap_output"
# quickomics_output_directory = "./PROT-975_Smed_EV_MV_OnePot_ProteoVista_for_publication_20241028140925/quickomics_files/all_contrasts_filtered_by_group"
# dd_zscore_thresh = 4
# remove_dd_in_de = T

#reformat protein-level expression data from MS-DAP - log2 transformed, normalized. Export for Quickomics upload.
quickomics_expression <- function(msdap_output_directory, quickomics_output_directory) {

    # Check that the ms-dap directory exists
    if (!dir.exists(msdap_output_directory)) {
        stop("The specified MS-DAP output directory does not exist: ", msdap_output_directory)
    }

    #check that the quickomics output directory exists
    if (!dir.exists(quickomics_output_directory)) {
        stop("The specified output directory does not exist: ", quickomics_output_directory)
    }

    # Look for the expected msd-dap protein-level quantitation file
    prot_quan_file <- list.files(
        msdap_output_directory,
        pattern = "protein_abundance__filter by group independently\\.tsv$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Handle case: file not found
    if (length(prot_quan_file) == 0) {
        stop("No matching protein abundance file found in: ", msdap_output_directory)
    }

    # Handle case: multiple matches
    if (length(prot_quan_file) > 1) {
        warning("Multiple matching files found; using the first one:\n", paste(prot_quan_file, collapse = "\n"))
        prot_quan_file <- prot_quan_file[1]
    }

    # Attempt to read the file
    prot_quan <- tryCatch({
        data.table::fread(prot_quan_file)
    }, error = function(e) {
        stop("Failed to read the protein abundance file: ", e$message)
    })

    #now modify and save out the quickomics formatted file
    quickomics_data <- prot_quan |>
        dplyr::select(-fasta_headers, -gene_symbols_or_id) |>
        dplyr::rename("UniqueID" = "protein_id")
    output_file <- file.path(quickomics_output_directory, "quickomics_expression_data.csv")
    write.csv(quickomics_data, output_file, row.names = FALSE)
}

#reformat MS-DAP Differential Expression and Differential Detection results. Filter DD hits as specified, based on |z-score| values and
#presences/absence in DE results. Export for Quickomics upload.
quickomics_de_dd_results <- function(msdap_output_directory, quickomics_output_directory, remove_dd_in_de = F, dd_zscore_thresh = 4) {

    # Check that the ms-dap directory exists
    if (!dir.exists(msdap_output_directory)) {
        stop("The specified MS-DAP output directory does not exist: ", msdap_output_directory)
    }

    #check that the quickomics output directory exists
    if (!dir.exists(quickomics_output_directory)) {
        stop("The specified output directory does not exist: ", quickomics_output_directory)
    }

    # Look for the expected ms-dap de results file
    de_file <- list.files(
        msdap_output_directory,
        pattern = "de_proteins\\.tsv.gz$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Handle case: file not found
    if (length(de_file) == 0) {
        stop("No matching DE results file found in: ", msdap_output_directory)
    }

    # Handle case: multiple matches
    if (length(de_file) > 1) {
        warning("Multiple matching DE results files found; using the first one:\n", paste(de_file, collapse = "\n"))
        de_file <- de_file[1]
    }

    # Attempt to read the file
    de_results <- tryCatch({
        data.table::fread(de_file)
    }, error = function(e) {
        stop("Failed to read the DE results file: ", e$message)
    })

    # Look for the expected ms-dap DD results file
    dd_file <- list.files(
        msdap_output_directory,
        pattern = "dd_proteins\\.tsv.gz$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Handle case: DD file not found
    if (length(dd_file) == 0) {
        stop("No matching DD results file found in: ", msdap_output_directory)
    }

    # Handle case: multiple DD matches
    if (length(de_file) > 1) {
        warning("Multiple matching DD results files found; using the first one:\n", paste(dd_file, collapse = "\n"))
        dd_file <- dd_file[1]
    }

    # Attempt to read the DD file
    dd_results <- tryCatch({
        data.table::fread(dd_file)
    }, error = function(e) {
        stop("Failed to read the DD results file: ", e$message)
    })

    #modify as needed for Quickomics
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

    #Now, split out by DEA algorithm and to allow saving for each method
    quickomics_de_list <- quickomics_de |> dplyr::group_split(dea_algorithm)
    de_names <- lapply(quickomics_de_list, function(x) unique(x$dea_algorithm))
    names(quickomics_de_list) <- as.character(unlist(de_names))

    #process/reformat DD results
    dd_results_full <- dd_results |>
        dplyr::select(protein_id, log2fc, contrast, zscore) |>
        dplyr::mutate(contrast = stringr::str_remove_all(contrast, pattern = paste(c(" # condition.*", "contrast: "), collapse = "|")),
                      log2fc = -log2fc,
                      Adj.P.Value = 0,
                      P.Value = 0,
                      dea_algorithm = "Differential.Detection") |>
        dplyr::filter(abs(zscore) >= dd_zscore_thresh) |>
        dplyr::select(-zscore) |>
        dplyr::mutate(contrast = stringr::str_replace(contrast, pattern = " vs ", "-")) |>
        dplyr::rename("UniqueID" = "protein_id",
                      "test" = "contrast",
                      "logFC" = "log2fc") |>
        dplyr::select(UniqueID, test, Adj.P.Value, P.Value, logFC, dea_algorithm)

    #now, for each contrast in the DE list object, add in the DD hits

    add_dd_to_de <- function(de_results, dd_results) {

        # #only keep differential detection results for proteins that weren't used for DE testing
        # de_prots <- de_results$UniqueID
        #
        # dd_results <- dd_results |> dplyr::filter(!UniqueID %in% de_prots)

        #get min pval and adj. p val for dea results
        min_pval <- min(de_results$P.Value, na.rm=TRUE)
        #print(min_pval)
        min_adjpval <- min(de_results$Adj.P.Value, na.rm=TRUE)
        #print(min_adjpval)

        #replace zero's with min values from DE results to get properly displayed in Quickomics plots.
        dd_results <- dd_results |>
            dplyr::mutate(P.Value = dplyr::if_else(dea_algorithm=="Differential.Detection", min_pval, P.Value),
                          Adj.P.Value = dplyr::if_else(dea_algorithm=="Diferential.Detection", min_adjpval, Adj.P.Value))

        #combine the results by binding rows together.
        combined <- de_results |> dplyr::bind_rows(dd_results)
        return(combined)
    }

    quickomics_de_list <- lapply(quickomics_de_list, add_dd_to_de, dd_results = dd_results_full)

    # Now, the list object contains all DE hits from each DE algorithm, plus each contains the full DD results.
    # Objective now is to go into each object in the list, group by each test/contrast
    # and if within that contrast the protein has DE testing results, then remove the Differential Detection results. This should be rare.

    if(remove_dd_in_de) {
        print("For each contrast, removing any differential detection results if the protein had differential expression results.")
        quickomics_de_list <- lapply(quickomics_de_list, remove_dd_if_de_exists)
    }

    #save out each csv file
    paths <- lapply(de_names, function(x) paste0(dirname(msdap_output_directory), "/quickomics_files/", x, "_dea_results.csv")) |> unlist()
    purrr::map2(quickomics_de_list, paths, .f = function(x, y) write.csv(x, file = y, row.names = F))

}

#reformat MS-DAP samples file for quickomics metadata table
quickomics_metadata <- function(msdap_output_directory, quickomics_output_directory) {

    #Get the Sample Metadata Table - load the sample table then filter for the samples in the contrast

    # Check that the ms-dap directory exists
    if (!dir.exists(msdap_output_directory)) {
        stop("The specified MS-DAP output directory does not exist: ", msdap_output_directory)
    }

    #check that the quickomics output directory exists
    if (!dir.exists(quickomics_output_directory)) {
        stop("The specified output directory does not exist: ", quickomics_output_directory)
    }

    # Look for the expected ms-dap samples file
    metadata_file <- list.files(
        msdap_output_directory,
        pattern = "^samples.tsv.gz$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Handle case: file not found
    if (length(metadata_file) == 0) {
        stop("No matching metadata file found in: ", msdap_output_directory)
    }

    # Handle case: multiple matches
    if (length(metadata_file) > 1) {
        warning("Multiple matching metadata files found; using the first one:\n", paste(metadata_file, collapse = "\n"))
        metadata_file <- metadata_file[1]
    }

    # Attempt to read the file
    sample_table <- tryCatch({
        data.table::fread(metadata_file)
    }, error = function(e) {
        stop("Failed to read the metadata file: ", e$message)
    })

    #remove any samples that were excluded
    sample_table_filtered <- sample_table |> dplyr::filter(!exclude)

    #now, load the quickomics protein expression data file and make sure that the column names match the metadata
    # Look for the expected ms-dap samples file
    prot_quan_file <- list.files(
        quickomics_output_directory,
        pattern = "quickomics_expression_data.csv$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Attempt to read the file
    quickomics_data <- tryCatch({
        data.table::fread(prot_quan_file)
    }, error = function(e) {
        stop("Failed to read the Quickomics expression data file: ", e$message)
    })

    columns_present <- colnames(quickomics_data)
    test <- length(intersect(sample_table_filtered$sample_id, columns_present)) == length(unique(sample_table_filtered$sample_id))

    #if not true, then throw an error and stop
    if(!test) {
        stop("Columns in sample metadata not matching to names in the expression data. Check dataset$samples$sample_id values and compare with values in column names for protein abundance global filter dataset.")
    } else (
        print("All samples in are present in the expression data. Proceeding with extraction.")
    )

    #also need to have the de/dd results to extract all unique contrast names after formatting
    quickomics_de_files <- list.files(
        quickomics_output_directory,
        pattern = "dea_results.csv$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Attempt to read the files and combine into one
    quickomics_de <- tryCatch({
        lapply(quickomics_de_files, data.table::fread) |> dplyr::bind_rows()
    }, error = function(e) {
        stop("Failed to read the Quickomics differential expression files: ", e$message)
    })

    #Modify columns to match expected for quickomics upload
    #Sample MetaData file. Must have columns "sampleid" and "group", "Order" and "ComparePairs" columns
    quickomics_md <- jmv_mixedLengthDF(list(sampleid = sample_table_filtered$sample_id,
                                            group = sample_table_filtered$group,
                                            Order = unique(sample_table_filtered$group),
                                            ComparePairs = unique(quickomics_de$test)))

    print('Now saving metadata file: ')
    output_file = file.path(quickomics_output_directory, "quickomics_sample_metadata.csv")
    print(output_file)

    write.csv(quickomics_md, output_file, row.names=FALSE)
}

#generate Gene/Protein table for Quickomics
quickomics_gene_table <- function(msdap_output_directory, quickomics_output_directory) {

    # Check that the ms-dap directory exists
    if (!dir.exists(msdap_output_directory)) {
        stop("The specified MS-DAP output directory does not exist: ", msdap_output_directory)
    }

    #check that the quickomics output directory exists
    if (!dir.exists(quickomics_output_directory)) {
        stop("The specified output directory does not exist: ", quickomics_output_directory)
    }

    # Look for the expected msd-dap protein-level quantitation file
    proteins_file <- list.files(
        msdap_output_directory,
        pattern = "^proteins.tsv.gz$",
        recursive = TRUE,
        full.names = TRUE
    )

    # Handle case: file not found
    if (length(proteins_file) == 0) {
        stop("No proteins.tsv.gz file found in: ", msdap_output_directory)
    }

    # Handle case: multiple matches
    if (length(proteins_file) > 1) {
        warning("Multiple proteins.tsv.gz files found; using the first one:\n", paste(proteins_file, collapse = "\n"))
        proteins_file <- proteins_file[1]
    }

    # Attempt to read the file
    proteins_list <- tryCatch({
        data.table::fread(proteins_file)
    }, error = function(e) {
        stop("Failed to read the proteins file: ", e$message)
    })

    #Gene/Protein Name File - must have four columns: id (sequential numbers),
    #UniqueID (must match the rownames in expression data and comparison data), Gene.Name, Protein.ID. Can have additional

    quickomics_gene_table <- proteins_list |>
        dplyr::mutate(id = dplyr::row_number()) |>
        dplyr::rename("UniqueID" = "protein_id",
                      "Gene.Name" = "gene_symbols_or_id",
                      "Protein.ID" = "accessions")

    output_file <- file.path(quickomics_output_directory, "quickomics_gene_protein_table.csv")
    write.csv(quickomics_gene_table, output_file, row.names = F)
}

#from ProteoVista_cerebro - same as ProteoVista original
add_dd_to_de <- function(de_results, dd_results) {

    # #only keep differential detection results for proteins that weren't used for DE testing
    # de_prots <- de_results$UniqueID
    #
    # dd_results <- dd_results |> dplyr::filter(!UniqueID %in% de_prots)

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

#From ProteoVista original
#quickomics expression set subfolders
# quickomics_expression_sets <- function(msdap_output_directory, msdap_data) {
#
#     #DEBUGGING 12/16/25
#     #msdap_output_directory <- "./ProteoVista_output/SIMR_ShortCourse_Jan2025_All_Combined_20250114211908/msdap_output/"
#     #load("./ProteoVista_output/SIMR_ShortCourse_Jan2025_All_Combined_20250114211908/msdap_output/2025-01-14_21-24-59/dataset.RData")
#     # #create quickomics output folder in parent folder of msdap_output folder
#     # dir.create(paste0(dirname(msdap_output_directory), "/quickomics_files/"))
#
#     #DEBUGGING - are all DD hits being applied to all tests?
#     # msdap_output_directory <- "./PROT-975_Smed_EV_MV_OnePot_ProteoVista_for_publication_20241028140925/msdap_output/"
#     # load("./PROT-975_Smed_EV_MV_OnePot_ProteoVista_for_publication_20241028140925/msdap_output/2024-10-28_14-11-19/dataset.RData")
#     #
#     #load("./ProteoVista_output/PROT-1027 Mouse Amyloid Fractionation_20240918140518/msdap_output/2024-09-18_14-06-23/dataset.RData")
#     #load("./ProteoVista_output/PROT-1027 Mouse Amyloid Fractionation_20240918140518/msdap_output/ExpressionSet_proteins_contrast wt_ev vs wt_mv.RData")
#     dataset = msdap_data
#
#     #need to load in the dd_proteins.tsv.gz results to filter and use later on.
#     diff_detect_file <- list.files(msdap_output_directory, pattern = "dd_proteins.tsv.gz", recursive = TRUE, full.names=TRUE)
#     diff_detect_results <- data.table::fread(diff_detect_file)
#
#     #example dd results for debugging
#     #diff_detect_results <- data.table::fread("./ProteoVista_output/SIMR_ShortCourse_Jan2025_All_Combined_20250114211908/msdap_output/2025-01-14_21-24-59/dd_proteins.tsv.gz")
#
#     #get the list of protein abundance data files for each contrast
#     prot_abund_by_contrast <- list.files(path=msdap_output_directory, pattern = "protein_abundance__filter by contrast", recursive=TRUE, full.names = TRUE)
#
#     #get the contrast names
#     contrasts <- basename(stringr::str_remove_all(prot_abund_by_contrast, pattern = "protein_abundance__filter by contrast; ")) |>
#         stringr::str_remove_all(pattern = ".tsv") |>
#         stringr::str_replace_all(pattern = " ", "_")
#
#
#     #Now, repeat above but for global proteins for all contrasts.
#     #Modified 9/18/24 to use the protein abundances filtered by groups instead of globally.
#     process_global_for_quickomics <- function(msdap_output_directory) {
#
#         #create global quickomics output folder
#         global_output <- paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/")
#         dir.create(global_output)
#
#         #msdap_output_directory <- "./ProteoVista_output/PROT-1007 Smed X1 Subsort OnePot Tubulin and Actin, Non-redundant Smed DB 2024 Eric Ross_20240918085352/msdap_output/"
#
#         #Step 1: load the global proteomics results and format/save out
#         prot_quan_file <- list.files(msdap_output_directory, pattern = "protein_abundance__filter by group independently.tsv", recursive = TRUE, full.names=TRUE)
#         #prot_quan_file <- data.table::fread("./ProteoVista_output/PROT-1007 Smed X1 Subsort Cytometry OnePot, Tubulin and Actin_20240917144856/msdap_output/2024-09-17_14-52-44/protein_abundance__global data filter.tsv")
#         prot_quan <- data.table::fread(prot_quan_file)
#
#         quickomics_data <- prot_quan |>
#             dplyr::select(-fasta_headers, -gene_symbols_or_id) |>
#             dplyr::rename("UniqueID" = "protein_id")
#         write.csv(quickomics_data, paste0(global_output, "/global_quickomics_expression_data.csv"), row.names = FALSE)
#
#         #Step 2:
#         #Step 2: get the DEA results, filter for contrast, export for each algorithm used.
#         de_results <- data.table::fread(list.files(path = msdap_output_directory, pattern = "de_proteins.tsv.gz", recursive=TRUE, full.names=TRUE))
#
#         #select protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm. Filter for dea_algorithm of choice. Modify column names to match expected inputs for quickomics.
#         #Required column names are UniqueId (protein group), test (comparison), Adj.P.Value, P.Value, and logFC
#         quickomics_de <- de_results |>
#             dplyr::select(protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm) |>
#             dplyr::rename("UniqueID" = "protein_id",
#                           "test" = "contrast",
#                           "Adj.P.Value" = "qvalue",
#                           "P.Value" = "pvalue",
#                           "logFC" = "foldchange.log2") |>
#             dplyr::mutate(test = stringr::str_remove_all(test, pattern = paste(c(" # condition.*", "contrast: "), collapse = "|")),
#                           test = stringr::str_replace(test, pattern = " vs ", "-"),
#                           logFC = -logFC)   #this reverses the standard notation for FC that MS-DAP uses (FC = A/B instead of conventional FC = B/A)
#
#
#         #Now, split out by test and save csv's for each test method
#         quickomics_de_list <- quickomics_de |> dplyr::group_split(dea_algorithm)
#         de_names <- lapply(quickomics_de_list, function(x) unique(x$dea_algorithm))
#         names(quickomics_de_list) <- as.character(unlist(de_names))
#
#         #before saving out, bind rows for differentially detected proteins for this contrast to all de outputs.
#         #step 1: filter the dd results and only retain the differential detected with |z-score| >= 4, select necessary columns and rename to match the DE results columns
#         dd_contrast <- diff_detect_results |>
#             dplyr::select(protein_id, log2fc, contrast, zscore) |>
#             dplyr::mutate(contrast = stringr::str_remove_all(contrast, pattern = paste(c(" # condition.*", "contrast: "), collapse = "|")),
#                           log2fc = -log2fc,
#                           Adj.P.Value = 0,
#                           P.Value = 0,
#                           dea_algorithm = "Differential.Detection") |>
#             dplyr::filter(abs(zscore) >= 4) |>
#             dplyr::select(-zscore) |>
#             dplyr::mutate(contrast = stringr::str_replace(contrast, pattern = " vs ", "-")) |>
#             dplyr::rename("UniqueID" = "protein_id",
#                           "test" = "contrast",
#                           "logFC" = "log2fc") |>
#             dplyr::select(UniqueID, test, Adj.P.Value, P.Value, logFC, dea_algorithm)
#
#         quickomics_de_list <- lapply(quickomics_de_list, add_dd_to_de, dd_results = dd_contrast)
#
#         #now, if specified, apply the dd filter to remove any dd hits that were detected in de results - applied per test/contrast
#         if()
#
#         #save out each csv file
#         paths <- lapply(de_names, function(x) paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/", x, "test_global_dea_results.csv")) |> unlist()
#         purrr::map2(quickomics_de_list, paths, .f = function(x, y) write.csv(x, file = y, row.names = F))
#
#         #Step 3: Get the Sample Metadata Table - load the sample table then filter for the samples in the contrast
#
#         sample_table <- data.table::fread(list.files(msdap_output_directory, pattern = "samples.tsv.gz", recursive=TRUE, full.names=TRUE))
#
#         #remove any samples that were excluded
#         #samples <- data.table::fread("./ProteoVista_output/PROT-1007 Smed X1 Subsort Cytometry OnePot, Tubulin and Actin_20240917144856/msdap_output/2024-09-17_14-52-44/samples.tsv.gz")
#
#         sample_table_filtered <- sample_table |> dplyr::filter(!exclude)
#
#         columns_present <- colnames(quickomics_data)
#         test <- length(intersect(sample_table_filtered$sample_id, columns_present)) == length(unique(sample_table_filtered$sample_id))
#
#         #if not true, then throw an error and stop
#         if(!test) {
#             stop("Columns in sample metadata not matching to names in the expression data. Check dataset$samples$sample_id values and compare with values in column names for protein abundance global filter dataset.")
#         } else (
#             print("All samples in are present in the expression data. Proceeding with extraction.")
#         )
#
#         #Modify columns to match expected for quickomics upload
#         #Sample MetaData file. Must have columns "sampleid" and "group", "Order" and "ComparePairs" columns
#         quickomics_md <- jmv_mixedLengthDF(list(sampleid = sample_table_filtered$sample_id,
#                                                 group = sample_table_filtered$group,
#                                                 Order = unique(sample_table_filtered$group),
#                                                 ComparePairs = unique(quickomics_de$test)))
#         quickomics_md
#         print('Now saving metadata file: ')
#         print(paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/global_quickomics_sample_metadata.csv"))
#
#         write.csv(quickomics_md, paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/global_quickomics_sample_metadata.csv"), row.names=FALSE)
#
#         #Step 4: Gene name table
#         gene_table <- dataset$proteins
#
#         #Gene/Protein Name File - must have four columns: id (sequential numbers), UniqueID (must match the rownames in expression data and comparison data), Gene.Name, Protein.ID. Can have additional
#
#         quickomics_gene_table <- gene_table |>
#             dplyr::mutate(id = dplyr::row_number()) |>
#             dplyr::rename("UniqueID" = "protein_id",
#                           "Gene.Name" = "gene_symbols_or_id",
#                           "Protein.ID" = "accessions")
#         write.csv(quickomics_gene_table, paste0(dirname(msdap_output_directory), "/quickomics_files/all_contrasts_filtered_by_group/global_quickomics_gene_protein_table.csv"), row.names = F)
#
#
#     }
#
#     process_global_for_quickomics(msdap_output_directory)
#
#
# }


#function to apply to a data frame - group by 'test' (contrast) and remove DD hits if there are DE results
remove_dd_if_de_exists <- function(x) {

    #pseudocode: group by contrast (test value) and UniqueID (Protein/Protein Group).
    #Capture number of unique dea_algorithm values - if it was a hit in multiple (both DE and DD) then set TRUE
    #Capture string of all algorithms it was a hit for. This is just an intermediate column for sanity check during debugging and is not returned
    #Filter to only keep rows that were a single hit (detected de and dd == F), or if TRUE, then only remove the DD value.
    nrows_start <- nrow(x)

    filtered_df <- x |>
        dplyr::group_by(test, UniqueID) |>
        dplyr::mutate(detected_de_and_dd = length(unique(dea_algorithm)) > 1,
                      algorithms_detected = toString(unique(dea_algorithm))) |>
        dplyr::ungroup() |>
        dplyr::filter(!(detected_de_and_dd & dea_algorithm == "Differential.Detection"))

    nrows_end <- nrow(filtered_df)
    nrow_removed = nrows_start - nrows_end
    pct_double_detected = round(100*(nrows_start-nrows_end)/nrows_start, digits=2)
    print(paste0("Removed ", nrow_removed, " rows, (", pct_double_detected, "% of input), with proteins/protein groups with values in both DE and DD results."))
    print("In these rows, the DE results were retained, and the additional Differential Detection rows were removed.")
    return(filtered_df)
}




