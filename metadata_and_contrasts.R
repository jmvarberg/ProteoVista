#Step 2: Metadata and Contrast Specifications
# Specification of Sidebar UI elements --------------------------------------------

#UI element for users to upload metadata file after modifying the template that was saved out.
output$inputMetadata <- renderUI({

    fileInput(inputId = "modified_metadata",
              label = "Select User-Modified Sample Metadata File (.xlsx)",
              multiple = F,
              accept = c(".xlsx")
    )
})

#UI element to allow input for msdap::analysis_quickstart() function for processing dataset.
output$mdsap_run_params <- renderUI({

    tagList(
        numericInput(
            inputId = "msdap_min_pep_per_prot",
            label = "Protein Filter: Minimum number of peptides per protein",
            value = 2,
            min = 1,
            width = "100%"
        ),
        checkboxInput(
            inputId = "msdap_filter_by_contrast",
            label = "Filter By Contrast (If unchecked, filtering and normalizations applied globally; if checked, filtering and normalization applied on subset of relevant samples within each specified contrast pair)",
            value = TRUE,
            width = "100%"
        ),
        selectInput(
            inputId = "msdap_norm_method",
            label = "Normalization Method",
            width = "100%",
            choices = list(
                "Variance Stabilizing Normalization (VSN)" = "vsn",
                "Variation Within, Mode Between" = "vwmb",
                "Mode Within, Mode Between" = "mwmb",
                "Median" = "median",
                "Loess" = "loess",
                "Robust Linear Regression (RLR)" = "rlr",
                "MS-EmpiRe" = "msempire",
                "Mode Between (Peptide)" = "mb_pep",
                "Mode Between (Protein)" = "mb_prot"
            ),
            selected = "vsn"
        ),
        checkboxInput("add_mbprot", "Add Mode Between (Protein) after Peptide-level normalization? (Recommended)", TRUE),
        numericInput(
            inputId = "dea_qval_thresh",
            label = "Q-value Threshold (Differential Expression Significance)",
            value = 0.05,
            width = "100%",
            min = 0,
            max = 1,
            step = 0.01
        ),
        numericInput(
            inputId = "dea_log2sig_thresh",
            label = "Absolute log2-Fold Change Significance Threshold",
            value = 0,
            width = "100%",
            min = 0
        )

    ) #close tagList()

}) #close renderUI()

#UI element to specify peptide filtering criteria for DEA
output$mdsap_dea_filt_params <- renderUI({

    tagList(
        numericInput(
            inputId = "msdap_min_detect_dea",
            label = "DEA Peptide Filter (Confidence Score): Minimum number of samples per group with high confidence ID",
            value = 3,
            width = "100%"
        ),
        numericInput(
            inputId = "msdap_min_quant_dea",
            label = "DEA Peptide Filter (Quantitation): Minimum number of samples per group with quantitative value",
            value = 3,
            width = "100%"
        )
    ) #close tagList

        ##Commenting out for now 9/17/24 JMV - just force people to specify number of samples/condition instead of option for fraction.
        # numericInput(
        #     inputId = "msdap_fraction_detect",
        #     label = "Peptide Filter (Confidence Score): Minimum fraction of samples per group with high confidence ID",
        #     value = 0,
        #     width = "100%",
        #     step = 0.1,
        #     min = 0,
        #     max = 1
        # ),
        # numericInput(
        #     inputId = "msdap_fraction_quant",
        #     label = "Peptide Filter (Quantitation): Minimum fraction of samples per group with quantitative value",
        #     value = 0,
        #     min = 0,
        #     max = 1,
        #     step = 0.1,
        #     width = "100%"
        # )

}) #close renderUI()

#UI element to specify peptide filtering criteria for Differential Detection analysis (dd)
output$mdsap_dd_filt_params <- renderUI({

    #diffdetect_min_peptides_observed = 2, #minimum number of peptides for a protein to be included in diff. detection analyses (protein level filter)
    #diffdetect_min_samples_observed = 3, #minimum number of samples within a condition that a peptide must be detected in to be kept for diff. detection (peptide level filter)


    tagList(
        numericInput(
            inputId = "msdap_min_pept_dd",
            label = "DD Protein Filter: Minimum number of peptides per protein required to keep protein for Differential Detection:",
            value = 2,
            width = "100%"
        ),
        numericInput(
            inputId = "msdap_min_samples_dd",
            label = "DD Peptide Filter: Minimum number of samples per group required to use peptide for Differential Detection:",
            value = 3,
            width = "100%"
        )
    ) #close tagList

    ##Commenting out for now 9/17/24 JMV - just force people to specify number of samples/condition instead of option for fraction.
    # numericInput(
    #     inputId = "msdap_fraction_detect",
    #     label = "Peptide Filter (Confidence Score): Minimum fraction of samples per group with high confidence ID",
    #     value = 0,
    #     width = "100%",
    #     step = 0.1,
    #     min = 0,
    #     max = 1
    # ),
    # numericInput(
    #     inputId = "msdap_fraction_quant",
    #     label = "Peptide Filter (Quantitation): Minimum fraction of samples per group with quantitative value",
    #     value = 0,
    #     min = 0,
    #     max = 1,
    #     step = 0.1,
    #     width = "100%"
    # )

}) #close renderUI()

#UI element for setting up contrasts
output$contrastMethod <- renderUI({

    selectInput(
        inputId = "contrast_method",
        label = "Contrast Design",
        choices = list("All vs. Reference" = "contrast_ref",
                       "All Pairwise Comparisons" = "contrast_pairwise",
                       "User Defined Comparisons" = "constrast_custom"),
        selected = "contrast_pairwise",
        width = "100%")

})

#UI element to specify DE algorithms to use
output$dea_selection <- renderUI({

    selectInput(
        inputId = "de_algorithm",
        label = "DE Algorithm(s)",
        choices = list("DEqMS" = "deqms",
                       "msEmpiRe" = "msempire",
                       "MSqRob" = "msqrob"),
        multiple=TRUE,
        selected = c("deqms", "msempire", "msqrob"),
        width = "100%")


})

#TO DO: add drop down list of column names from uploaded metadata that allows user to specify random variables/batch effects to account for during regression analysis

# Processing and Reactivity -----------------------------------------------

#read in the metadata into an object to display in the main panel as a data table for visual confirmation
final_metadata <- reactive({
    req(input$modified_metadata)
    readxl::read_xlsx(input$modified_metadata$datapath)
})

#generate output Data Table for visualizing the User-modified metadata table.
output$finalMeta <- DT::renderDataTable(final_metadata())

#generate dynamic UI element to select reference group
output$referenceSelect <- renderUI({

    #from the final_metadata, get all unique group labels and coerrce to character if needed
    user_md <- final_metadata()
    groups <- unique(user_md$group)

    selectInput(inputId = "reference_group",
                label = "Select Reference Group for Comparisons",
                choices = groups,
                multiple = F)


})

#if custom is selected, then generate UI element to allow user to select from all pairwise comparisons.
output$contrastSelection <- renderUI({

    groups <- metadata_groups()
    contrasts <- contrast_pairs <- apply(combn(groups,2), 2, paste, collapse='_vs_')

    req(input$contrast_method)
    pickerInput(inputId = "user_contrasts",
                label = "Contrasts List: ",
                choices = contrasts,
                multiple=TRUE,
                selected=NULL,
                width = '100%',
                options = list(
                    title = "Select Contrasts To Use",
                    `actions-box` = TRUE))

})

#reactive object of all unique groups in user uploaded metadata with checks for validity
metadata_groups <- reactive({

    #generate based on selection of UI input for contrast method
    req(input$contrast_method)

    #from the final_metadata, get all unique group labels and coerrce to character if needed
    user_md <- final_metadata()
    groups <- unique(user_md$group)

    if(!is.character(groups)) {
        print("Coercing group labels from metadata to character type.")
        groups <- as.character(groups)
    }

    #check input: make sure that all samples have group labels
    missing_groups <- sum(is.na(groups)) + sum(groups == "" | groups == "NA", na.rm =TRUE)

    if(missing_groups > 0){

        #alert user if not all samples have group labels.
        shinyalert::shinyalert(title = "Error: Detected Samples Without Group Assignment",
                               text = "Please correct and re-upload your metadata file. All samples must be assigned to a group for processing.",
                               type = "error")
    }

    #check input: make sure that at least two groups are specified
    if(length(groups) < 2) {

        #alert user if there are not at least two groups detected
        shinyalert::shinyalert(title = "Error: No Groups To Compare",
                               text = "Please correct and re-upload your metadata file. Only one group was specified, unable to do contrasts.",
                               type = "error")

    }

    #check input: do we need/want a filter to require at least N number of samples per group? Is there a requirement for MS-DAP?

    print("Groups:")
    print(groups)
    return(groups)
})

#make output UI object with info about the number and identity of groups in metadata
output$groupsCheck <- renderText({
    groups <- metadata_groups()
    format_html_list(groups)
})

#now, handle how to select and filter
selected_contrast_pairs <- reactive({

    #if user specifies contrasting all groups to a reference, generate ui element to select what the reference group is.
    if(input$contrast_method == "contrast_ref") {

        req(input$reference_group)

        #generate list of all comparisons to the reference group
        comparisons <- generate_contrasts_to_ref(metadata_groups(), input$reference_group)

        list_of_vectors <- lapply(comparisons, function(x) strsplit(x, "__")[[1]])
        #contrasts_out <- c(unlist(lapply(list_of_vectors, print_list)))
        #print(contrasts_out)
        return(list_of_vectors)

    } else if(input$contrast_method == "contrast_pairwise") {
        comparisons <- apply(combn(metadata_groups(),2), 2, paste, collapse='__')
        list_of_vectors <- lapply(comparisons, function(x) strsplit(x, "__")[[1]])
        #contrasts_out <- c(unlist(lapply(list_of_vectors, print_list)))
        #print(contrasts_out)
        return(list_of_vectors)

    } else if(input$contrast_method == "constrast_custom") {

        user_selections <- input$user_contrasts
        list_of_vectors <- lapply(user_selections, function(x) strsplit(x, "_vs_")[[1]])
        return(list_of_vectors)
    }
})

#generate output text object listing the selected contrast pairs
output$finalContrasts <- renderText({
    filtered_list <- selected_contrast_pairs()
    format_html_list(c(unlist(lapply(filtered_list, print_list))))
})

#once user specifies, run MS-DAP processing using the uploaded dataset and metadata
#browser()

observeEvent(input$submit_msdap, {

    req(input$modified_metadata, input$data, input$database)

    dataset <- msdap_dataset$dataset
    #import the final metadata and add to the dataset object
    dataset <- msdap::import_sample_metadata(dataset, filename = input$modified_metadata$datapath)

    #define contrasts
    dataset <- msdap::setup_contrasts(dataset, contrast_list = selected_contrast_pairs(), random_variables = )

    #get the project_dir path
    projectDir <- project_dir()
    print(projectDir)

    waiter_show(html = data_processing_waiting)

    #check that the metadata has been added to the msdap_dataset object
    check_samples <- length(dataset$samples)

    if(check_samples < 1) {
        shinyalert::shinyalert(title = "Metadata Missing?",
                               text = "Looks like you have tried running MS-DAP before you have submitted your sample metadata. Please upload and submit the final metadata table and then try again.",
                               type = "error")
    }

    #prep input for normalization methods
    norm_to_use <- input$msdap_norm_method

    #if mbprot selected, then add to this
    if(input$add_mbprot) {
        norm_to_use <- c(norm_to_use, "modebetween_protein")
    }

    #make ms-dap output directory inside of project_dir()
    msdap_dir <- paste0(projectDir, "/msdap_output/")
    dir.create(msdap_dir)

    #sink processing into txt file
    sink(paste0(msdap_dir, "sink_console_output.txt"), append = TRUE)

    #check if genes provided for filtering or not
    genes_to_filter <- input$geneFilter
    #genes_to_filter <- stringr::str_remove_all(as.character(unlist(stringr::str_split(genes_to_filter, ","))), pattern = " ")
    print(genes_to_filter)

    if(genes_to_filter == "") {
        #filter dataset according to provided regex
        dataset <- msdap::remove_proteins_by_name(dataset=dataset, remove_irt_peptides = F, regular_expression = input$regexFilter)

    } else {

        #filter dataset according to user provided regex and/or gene/symbol filters
        dataset <- msdap::remove_proteins_by_name(dataset=dataset, remove_irt_peptides = F, regular_expression = input$regexFilter, gene_symbols = input$geneFilter)

    }




    #run quickstart analysis for MS-DAP
    dataset <- msdap::analysis_quickstart(
        dataset,
        filter_min_detect = input$msdap_min_detect_dea, #inputId = msdap_min_detect_dea
        #filter_fraction_detect = input$msdap_fraction_detect, #inputId = msdap_fraction_detect
        filter_min_quant = input$msdap_min_quant_dea, #inputId = msdap_min_quant
        #filter_fraction_quant = input$msdap_fraction_quant, #inputId = msdap_fraction_quant
        filter_min_peptide_per_prot = input$msdap_min_pep_per_prot, #inputId = msdap_min_pep_per_prot
        filter_topn_peptides = 0, #set to zero to disable topN filtering - not used with MaxLFQ roll-up
        filter_by_contrast = input$msdap_filter_by_contrast, #inputId = msdap_filter_by_contrast
        norm_algorithm = norm_to_use,
        rollup_algorithm = "maxlfq", #other options are sum, maxlfq-diann, and tukey median polish.
        dea_algorithm = input$de_algorithm,
        dea_qvalue_threshold = input$dea_qval_thresh,
        dea_log2foldchange_threshold = input$dea_log2sig_thresh,
        diffdetect_min_peptides_observed = input$msdap_min_pept_dd, #
        diffdetect_min_samples_observed = input$msdap_min_samples_dd, #
        #diffdetect_min_fraction_observed = 0.5,
        pca_sample_labels = "auto",
        var_explained_sample_metadata = NULL,
        multiprocessing_maxcores = 10,
        output_abundance_tables = TRUE,
        output_qc_report = TRUE,
        output_dir = msdap_dir,
        output_within_timestamped_subdirectory = TRUE,
        dump_all_data = TRUE
    )

    #once this finishes, then create the quickomics output

    print("Now extracting information for Quickomics Input Files...")

    # Quickomics Processing --------------------------------------------------------

    #create the quickomics output directory
    quickomics_dir <- paste0(projectDir, "/quickomics_files/")
    dir.create(quickomics_dir)
    print(paste0("Created output directory for Quickomics Files at ", quickomics_dir))

    #process quickomics extraction for all contrasts
    quickomics_expression_sets(msdap_output_directory = msdap_dir, msdap_data = dataset)


    # #Step 1: Get the Normalized Protein Abundance Values - this is not in the dataset object, need to read in from the directory
    #
    # ##TO DO: Modify this to detect all expression sets. Scrape name from the expression set. Create separate output folders inside quickomics path for each expression set. Pull in quan for each instead of global.
    # #check that the protein quantitation file exists
    # #msdap_dir <- "./test_output/test_20240625163653/msdap_output/"
    # msdap_files <- list.files(msdap_dir, full.names=TRUE, recursive = TRUE)
    # prot_quan_file <- msdap_files[stringr::str_detect(msdap_files, pattern = "protein_abundance__global")]
    # if(file.exists(prot_quan_file)) {
    #     prot_quan <- data.table::fread(prot_quan_file)
    # } else (
    #     shinyalert::shinyalert(title = "Cannot Find MS-DAP Protein Quantitation File",
    #                            text = "Please check in the MS-DAP output folder for a file 'protein_abundance__global data filter.tsv",
    #                            type = "error")
    # )
    #
    # #Modify columns to match expected for quickomics upload
    # quickomics_data <- prot_quan |>
    #     dplyr::select(-fasta_headers, -gene_symbols_or_id) |>
    #     dplyr::rename("UniqueID" = "protein_id")
    # write.csv(quickomics_data, paste0(quickomics_dir, "quickomics_expression_data.csv"), row.names = FALSE)
    #
    # #Step 2: Get the DE test results and filter for test of interest
    # de_results <- dataset$de_proteins
    #
    # #To Do: select protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm. Filter for dea_algorithm of choice. Modify column names to match expected inputs for quickomics.
    # #Required column names are UniqueId (protein group), test (comparison), Adj.P.Value, P.Value, and logFC
    # quickomics_de <- de_results |>
    #     dplyr::select(protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm) |>
    #     dplyr::rename("UniqueID" = "protein_id",
    #                   "test" = "contrast",
    #                   "Adj.P.Value" = "qvalue",
    #                   "P.Value" = "pvalue",
    #                   "logFC" = "foldchange.log2") |>
    #     dplyr::mutate(test = stringr::str_remove_all(test, pattern = "contrast: "),
    #                   test = stringr::str_replace(test, pattern = " vs ", "-"),
    #                   logFC = -logFC) #this reverses the standard notation for FC that MS-DAP uses (FC = A/B instead of conventional FC = B/A)
    #
    # #Now, split out by test and save csv's for each test method
    # quickomics_de_list <- quickomics_de |> dplyr::group_split(dea_algorithm)
    # de_names <- lapply(quickomics_de_list, function(x) unique(x$dea_algorithm))
    # names(quickomics_de_list) <- as.character(unlist(de_names))
    #
    # #save out each csv file
    # paths <- lapply(de_names, function(x) paste0(quickomics_dir, x, "_dea_results.csv")) |> unlist()
    # purrr::map2(quickomics_de_list, paths, .f = function(x, y) write.csv(x, file = y, row.names = F))
    #
    # #Step 3: Get the Sample Metadata Table
    # sample_md <- dataset$samples
    #
    # #Make sure that the values match the column names in the protein expression data
    # sample_ids <- sample_md$sample_id
    # columns_present <- colnames(quickomics_data)
    # test <- length(intersect(sample_ids, columns_present)) == length(unique(sample_ids))
    #
    # #if not true, then throw an error and stop
    # if(!test) {
    #     stop("Columns in sample metadata not matching to names in the expression data. Check dataset$samples$sample_id values and compare with values in column names for protein abundance global filter dataset.")
    # }
    #
    # #Modify columns to match expected for quickomics upload
    # #Sample MetaData file. Must have columns "sampleid" and "group", "Order" and "ComparePairs" columns
    # quickomics_md <- jmv_mixedLengthDF(list(sampleid = sample_md$sample_id,
    #                                         group = sample_md$group,
    #                                         Order = unique(sample_md$group),
    #                                         ComparePairs = unique(quickomics_de$test)))
    #
    # write.csv(quickomics_md, paste0(quickomics_dir, "quickomics_sample_metadata.csv"), row.names=FALSE)
    #
    # #Step 4: Gene name table
    # gene_table <- dataset$proteins
    #
    # #Gene/Protein Name File - must have four columns: id (sequential numbers), UniqueID (must match the rownames in expression data and comparison data), Gene.Name, Protein.ID. Can have additional
    #
    # quickomics_gene_table <- gene_table |>
    #     dplyr::mutate(id = dplyr::row_number()) |>
    #     dplyr::rename("UniqueID" = "protein_id",
    #                   "Gene.Name" = "gene_symbols_or_id",
    #                   "Protein.ID" = "accessions")
    # write.csv(quickomics_gene_table, paste0(quickomics_dir, "quickomics_gene_protein_table.csv"), row.names = F)
    sink()
    waiter_hide()

    #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
    shinyalert::shinyalert(title = "Dataset Successfully Processed with MS-DAP",
                           text = paste0("MS-DAP ouput files are located at ", msdap_dir, ". See the next tab for overview QC plots and information."),
                           type = "success")


})

