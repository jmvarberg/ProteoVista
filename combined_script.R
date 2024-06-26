# Tab 1 Server Logic ------------------------------------------------------


# Specification of Sidebar UI elements --------------------------------------------

#UI for Selecting Primary Analysis Tool Used
output$primaryAnalysis <- renderUI({

    selectInput(
        inputId = "analysisTool",
        label = "Primary Analysis Software",
        choices = list("Spectronaut" = "spectronaut",
                       "FragPipe/IonQuant (DDA)" = "fp_ionquant",
                       "FragPipe/DIA-NN" = "fp_diann",
                       "ProteomeDiscoverer" = "pd"),
        selected = "spectronaut",
        width = "100%"
    )
})

#UI to dynamically change the presentation of the Upload File selection based on software selection
output$inputData <- renderUI({

    req(input$analysisTool)

    if(input$analysisTool == "fp_ionquant") {

        textInput(
            inputId = "data",
            label = "Path to FragPipe Output Directory",
            placeholder = "/n/proteomics/..."
        )

    } else if(input$analysisTool == "fp_diann") {

        fileInput(
            inputId = "data",
            label = "Select DIA-NN report TSV file (in diann-output folder in output directory)",
            multiple = F,
            accept = c(".tsv"),
            width = "100%"
        )
    } else if(input$analysisTool == "spectronaut") {

        fileInput(
            inputId = "data",
            label = "Select MS-DAP Formatted Spectronaut Report",
            multiple = F,
            accept = c(".tsv", ".csv"),
            width = "100%"
        )
    } else if(input$analysisTool == "pd") {

        fileInput(
            inputId = "data",
            label = "Select ProteomeDiscoverer PSMs.txt File",
            multiple = F,
            accept = c(".txt"),
            width = "100%"
        )
    }

})

#UI for Uploading FASTA database
output$inputDatabase <- renderUI({

    fileInput(inputId = "database",
              label = "Choose FASTA database file used for search",
              multiple = F,
              accept = ".fasta",
              width = "100%"
    )
})

#UI for adding project notes in textbox.
output$inputNotes <- renderUI({

    textAreaInput(inputId = "project_description",
                  label = "Project Description/Notes",
                  placeholder = "Sample information, processing notes, etc.",
                  rows = 10,
    )

})

#UI to handle primary analysis tool-specific MS-DAP upload parameters. Different ingestion functions have specific parameter options.
output$mdsap_input_params <- renderUI({

    req(input$analysisTool)
    print(input$analysisTool)

    if(input$analysisTool == "fp_ionquant") {
        tagList(
            h4("Selected Analysis Tool: FragPipe + IonQuant"),
            selectInput(inputId = "fp_mode",
                        label = "Select Acquisition Mode (DDA or DIA)",
                        choices = list("DDA" = "dda", "DIA" = "dia"),
                        multiple = FALSE),
            numericInput(inputId = "conf_thresh",
                         label = "Select PSM Identification Confidence (FDR) threshold value (default = 0.01)",
                         value = 0.01,
                         min = 0,
                         max = 1,
                         step = 0.01),
            selectInput(inputId = "pep_collapse",
                        label = "Select method to combine/collapse peptides (use modified or plain peptide sequences)?",
                        choices = list("Modified" = "mod", "Plain" = "plain"),
                        multiple = FALSE,
                        selected = "mod")

        )
    } else if(input$analysisTool == "spectronaut") {
        tagList(
            h4("Selected Analysis Tool: Spectronaut"),
            numericInput(inputId = "conf_thresh",
                         label = "Select Peptide Identification Confidence (FDR) threshold value (default = 0.01)",
                         value = 0.01,
                         min = 0,
                         max = 1,
                         step = 0.01),
            selectInput(inputId = "spec_irt",
                        label = "Use Standardized Retention Time (iRT, recommended), or Empirical Retention Time (RTEmpirical)",
                        multiple = F,
                        choices = list("iRT" = TRUE, "Emperical RT" = FALSE))
        )
    }  else if(input$analysisTool == "pd") {

        tagList(
            h4("Selected Analysis Tool: ProteomeDiscoverer"),
            numericInput(inputId = "conf_thresh",
                         label = "Select Peptide Identification Confidence (Percolator Q-value) threshold value (default = 0.01)",
                         value = 0.01,
                         min = 0, max = 1,
                         step = 0.01),
            selectInput(inputId = "pep_collapse",
                        label = "Select method to combine/collapse peptides (use modified or plain peptide sequences)?",
                        choices = list("Modified" = "mod", "Plain" = "plain"),
                        multiple = FALSE,
                        selected = "mod"),
            selectInput(inputId = "psm_per_precursor",
                        label = "Select how to to roll-up PSM quantitation to precursors.",
                        choices = list("Sum" = "sum", "Intensity" = "intensity", "Confidence" = "confidence"),
                        multiple = FALSE,
                        selected = "Sum"),
            #shinyBS::addTooltip(session = session, id = "psm_per_precursor", title = "This parameter allows you to control how abundance values from precursors matched by multiple PSM are handled, as this might depend on your ProteomeDiscoverer settings. If ProteomeDiscoverer performed peak integration and reports the same (redundant) peak intensity for each PSM of the same precursor, we suggest to use 'Intensity'. Set to 'Sum' to use the sum of all PSM intensity values per precursor*sample (default). Use 'one_psm_per_precursor = 'Intensity' to select the highest intensity value (within the subset of PSM where confidence < confidence_threshold). Use 'Confidence' to select the intensity value from the PSM with best/lowest confidence value",
            #placement = "bottom", options = list(container = "body")),
            shinyWidgets::switchInput(inputId = "remove_low_conf", label = "Remove low-confidence peptides?", value = TRUE, width = "auto", size = "small", labelWidth = "250px", handleWidth = "50px", inline = TRUE)

        )
    }  else if(input$analysisTool == "fp_diann") {
        tagList(
            h4("Selected Analysis Tool: FragPipe + DIA-NN"),
            numericInput(inputId = "conf_thresh",
                         label = "Select Peptide Identification Confidence (FDR) threshold value (default = 0.01)",
                         value = 0.01,
                         min = 0,
                         max = 1,
                         step = 0.01),
            selectInput(inputId = "spec_irt",
                        label = "Use Standardized Retention Time (iRT, recommended), or Empirical Retention Time (RTEmpirical)",
                        multiple = F,
                        choices = list("iRT" = TRUE, "Emperical RT" = FALSE))
        )
    }
})

# Specification of Main Panel UI Elements ---------------------------------

#HTML document for main panel
output$tab_1_md <- renderUI({
    includeHTML("markdown/proteovista_project_setup.html")
})

# Processing --------------------------------------------------------------

#make a reactive object for the project directory. This is used throughout for saving output and is based on the Project ID and date/time
project_dir <- reactive({
    #create output directory for the project using specified project name as directory in /n/proteomics/washburn/Joe/ProteoVista_output/
    date.time <- format(Sys.time(), "%Y%m%d%H%M%S")
    #outdir <- "/Volumes/proteomics/washburn/Joe/Utilities/ProteoVista_output/"
    outdir <- "./test_output/"
    projdir <- paste0(outdir, input$projectID, "_", date.time)


})


#Action: Reactive object for MS-DAP dataset object.
msdap_dataset <- NULL
makeReactiveBinding("msdap_dataset")

#Action: Once input is specified, run initial processing steps of MS-DAP upon user clicking "Submit Input Files" button
observeEvent(input$submit_input, {

    projdir <- project_dir()
    print(projdir)

    req(input$data, input$database, input$projectID)

    # #create output directory for the project using specified project name as directory in /n/proteomics/washburn/Joe/ProteoVista_output/
    # date.time <- format(Sys.time(), "%Y%m%d%H%M%S")
    # #outdir <- "/Volumes/proteomics/washburn/Joe/Utilities/ProteoVista_output/"
    # outdir <- "./test_output/"
    # projdir <- paste0(outdir, input$projectID, "_", date.time)
    dir.create(path = projdir)

    #create input_data subdirectory
    dir.create(path = paste0(projdir, "/input_data/"))

    #update the project_dir reactive object
    #project_dir <<- projdir

    #copy input files into the project directory
    file.copy(input$data$datapath, paste0(projdir, "/input_data/", input$data$name))
    file.copy(input$database$datapath, paste0(projdir, "/input_data/", input$database$name))

    if(input$analysisTool == "spectronaut") {

        waiter_show(html = data_ingest_waiting)

        # print(paste0("conf_thresh: ", input$conf_thresh))
        # print(paste0("iRT: ", input$spec_irt))
        # print(paste0(projdir, "/", input$data$name))

        #import the Spectronaut Report file that was uploaded by reading from the project directory
        dataset <- msdap::import_dataset_spectronaut(filename = paste0(projdir, "/input_data/", input$data$name), confidence_threshold = input$conf_thresh)

        #add the fasta database to the dataset
        dataset <- msdap::import_fasta(dataset, files = paste0(projdir, "/input_data/", input$database$name))

        #write out sample metadata template to projdir
        msdap::write_template_for_sample_metadata(dataset, filename = paste0(projdir, "/input_data/", input$projectID, "_sample_metadata_table.xlsx"), overwrite = F)

        waiter_hide()

        #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
        shinyalert::shinyalert(title = "Dataset Uploaded and Metadata Template Generated",
                               text = paste0("Please open and edit the metadata template file, located at ", paste0(projdir, "/", input$projectID, "_sample_metadata_table.xlsx"), " Once you have finished modifying the metadata table, save, and then proceed to Step 2 to continue your analysis."),
                               type = "success")

    } else if(input$analysisTool == "fp_diann") {

        waiter_show(html = data_ingest_waiting)

        #import the Spectronaut Report file that was uploaded by reading from the project directory
        dataset <- msdap::import_dataset_diann(filename = paste0(projdir, "/", input$data$name))

        #add the fasta database to the dataset
        dataset <- msdap::import_fasta(dataset, files = paste0(projdir, "/", input$database$name))

        #write out sample metadata template to projdir
        msdap::write_template_for_sample_metadata(dataset, filename = paste0(projdir, "/", input$projectID, "_sample_metadata_table.xlsx"), overwrite = F)

        waiter_hide()

        #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
        shinyalert::shinyalert(title = "Dataset Uploaded and Metadata Template Generated",
                               text = paste0("Please open and edit the metadata template file, located at ", paste0(projdir, "/", input$projectID, "_sample_metadata_table.xlsx"), " Once you have finished modifying the metadata table, save, and then proceed to Step 2 to continue your analysis."),
                               type = "success")

    } else if(input$analysisTool == "pd") {

        waiter_show(html = data_ingest_waiting)

        #import the Spectronaut Report file that was uploaded by reading from the project directory
        dataset <- msdap::import_dataset_diann(filename = paste0(projdir, "/", input$data$name))

        #add the fasta database to the dataset
        dataset <- msdap::import_fasta(dataset, files = paste0(projdir, "/", input$database$name))

        #write out sample metadata template to projdir
        msdap::write_template_for_sample_metadata(dataset, filename = paste0(projdir, "/", input$projectID, "_sample_metadata_table.xlsx"), overwrite = F)

        waiter_hide()

        #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
        shinyalert::shinyalert(title = "Dataset Uploaded and Metadata Template Generated",
                               text = paste0("Please open and edit the metadata template file, located at ", paste0(projdir, "/", input$projectID, "_sample_metadata_table.xlsx"), " Once you have finished modifying the metadata table, save, and then proceed to Step 2 to continue your analysis."),
                               type = "success")
    }

    msdap::print_dataset_summary(dataset)
    #assign("msdap_dataset", dataset, envir = .Global.env)

})


# Tab 2 Server Logic ------------------------------------------------------

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
            inputId = "msdap_min_detect",
            label = "Peptide Filter (Confidence Score): Minimum number of samples per group with high confidence ID",
            value = 0,
            width = "100%"
        ),
        numericInput(
            inputId = "msdap_min_quant",
            label = "Peptide Filter (Quantitation): Minimum number of samples per group with quantitative value",
            value = 0,
            width = "100%"
        ),
        numericInput(
            inputId = "msdap_fraction_detect",
            label = "Peptide Filter (Confidence Score): Minimum fraction of samples per group with high confidence ID",
            value = 0,
            width = "100%",
            step = 0.1,
            min = 0,
            max = 1
        ),
        numericInput(
            inputId = "msdap_fraction_quant",
            label = "Peptide Filter (Quantitation): Minimum fraction of samples per group with quantitative value",
            value = 0,
            min = 0,
            max = 1,
            step = 0.1,
            width = "100%"
        ),
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
        )

    ) #close tagList()

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

    dataset <- msdap_dataset()
    #import the final metadata and add to the dataset object
    dataset <- msdap::import_sample_metadata(dataset, filename = input$modified_metadata$datapath)

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
    #run quickstart analysis for MS-DAP
    dataset <- msdap::analysis_quickstart(
        dataset,
        filter_min_detect = input$msdap_min_detect, #inputId = msdap_min_detect
        filter_fraction_detect = input$msdap_fraction_detect, #inputId = msdap_fraction_detect
        filter_min_quant = input$msdap_min_quant, #inputId = msdap_min_quant
        filter_fraction_quant = input$msdap_fraction_quant, #inputId = msdap_fraction_quant
        filter_min_peptide_per_prot = input$msdap_min_pep_per_prot, #inputId = msdap_min_pep_per_prot
        filter_topn_peptides = 0, #inputId = ? not currently used, need to check if should be added. This is like topN for MaxLFQ?
        filter_by_contrast = input$msdap_filter_by_contrast, #inputId = msdap_filter_by_contrast
        norm_algorithm = norm_to_use,
        rollup_algorithm = "maxlfq",
        dea_algorithm = c("deqms", "msqrob", "msempire"),
        dea_qvalue_threshold = 0.01,
        dea_log2foldchange_threshold = 0.58,
        diffdetect_min_peptides_observed = 2,
        diffdetect_min_samples_observed = 2,
        diffdetect_min_fraction_observed = 0.5,
        pca_sample_labels = "auto",
        var_explained_sample_metadata = NULL,
        multiprocessing_maxcores = 8,
        output_abundance_tables = TRUE,
        output_qc_report = TRUE,
        output_dir = msdap_dir,
        output_within_timestamped_subdirectory = TRUE,
        dump_all_data = FALSE
    )

    #once this finishes, then create the quickomics output

    print("Now extracting information for Quickomics Input Files...")

    # Quickomics Processing --------------------------------------------------------

    #create the quickomics output directory
    quickomics_dir <- paste0(projectDir, "/quickomics_files/")
    dir.create(quickomics_dir)
    print(paste0("Created output directory for Quickomics Files at ", quickomics_dir))

    #Step 1: Get the Normalized Protein Abundance Values - this is not in the dataset object, need to read in from the directory
    #check that the protein quantitation file exists
    if(file.exists(paste0(msdap_dir, "protein_abundance__global data filter.tsv"))) {
        prot_quan <- data.table::fread(paste0(msdap_dir, "protein_abundance__global data filter.tsv"))
    } else (
        shinyalert::shinyalert(title = "Cannot Find MS-DAP Protein Quantitation File",
                               text = "Please check in the MS-DAP output folder for a file 'protein_abundance__global data filter.tsv",
                               type = "error")
    )

    #Modify columns to match expected for quickomics upload
    quickomics_data <- prot_quan |>
        dplyr::select(-fasta_headers, -gene_symbols_or_id) |>
        dplyr::rename("UniqueID" = "protein_id")
    write.csv(quickomics_data, paste0(quickomics_dir, "quickomics_expression_data.csv"), row.names = FALSE)

    #Step 2: Get the DE test results and filter for test of interest
    de_results <- dataset$de_proteins

    #To Do: select protein_id, pvalue, qvalue, contrast, foldchange.log2, dea_algorithm. Filter for dea_algorithm of choice. Modify column names to match expected inputs for quickomics.
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
                      logFC = -logFC) #this reverses the standard notation for FC that MS-DAP uses (FC = A/B instead of conventional FC = B/A)

    #Now, split out by test and save csv's for each test method
    quickomics_de_list <- quickomics_de |> dplyr::group_split(dea_algorithm)
    de_names <- lapply(quickomics_de_list, function(x) unique(x$dea_algorithm))
    names(quickomics_de_list) <- as.character(unlist(de_names))

    #save out each csv file
    paths <- lapply(de_names, function(x) paste0(quickomics_dir, x, "_dea_results.csv")) |> unlist()
    purrr::map2(quickomics_de_list, paths, .f = function(x, y) write.csv(x, file = y, row.names = F))

    #Step 3: Get the Sample Metadata Table
    sample_md <- dataset$samples

    #Make sure that the values match the column names in the protein expression data
    sample_ids <- sample_md$sample_id
    columns_present <- colnames(quickomics_data)
    test <- length(intersect(sample_ids, columns_present)) == length(unique(sample_ids))

    #if not true, then throw an error and stop
    if(!test) {
        stop("Columns in sample metadata not matching to names in the expression data. Check dataset$samples$sample_id values and compare with values in column names for protein abundance global filter dataset.")
    }

    #Modify columns to match expected for quickomics upload
    #Sample MetaData file. Must have columns "sampleid" and "group", "Order" and "ComparePairs" columns
    quickomics_md <- jmvtools::jmv_mixedLengthDF(list(sampleid = sample_md$sample_id,
                                                      group = sample_md$group,
                                                      Order = unique(sample_md$group),
                                                      ComparePairs = unique(quickomics_de$test)))

    write.csv(quickomics_md, paste0(quickomics_dir, "quickomics_sample_metadata.csv", row.names=FALSE))

    #Step 4: Gene name table
    gene_table <- dataset$proteins

    #Gene/Protein Name File - must have four columns: id (sequential numbers), UniqueID (must match the rownames in expression data and comparison data), Gene.Name, Protein.ID. Can have additional

    quickomics_gene_table <- gene_table |>
        dplyr::mutate(id = dplyr::row_number()) |>
        dplyr::rename("UniqueID" = "protein_id",
                      "Gene.Name" = "gene_symbols_or_id",
                      "Protein.ID" = "accessions")
    write.csv(quickomics_gene_table, paste0(quickomics_dir, "quickomics_gene_protein_table.csv"), row.names = F)

    waiter_hide()



    #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
    shinyalert::shinyalert(title = "Dataset Successfully Processed with MS-DAP",
                           text = paste0("MS-DAP ouput files are located at ", msdap_dir, ". See the next tab for overview QC plots and information."),
                           type = "success")

})


