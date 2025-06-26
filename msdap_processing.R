#msdap processing functions



#once user specifies, run MS-DAP processing using the uploaded dataset and metadata
#browser()

observeEvent(input$submit_msdap, {

    #require the modified metadata file, the data file, and the database to continue
    req(input$modified_metadata, input$data, input$database)

    dataset <- msdap_dataset$dataset

    #get the project_dir path
    projectDir <- project_dir()
    print(projectDir)

    #initialize HTML window for processing visualization
    waiter_show(html = data_processing_waiting)

    # #check that the metadata has been added to the msdap_dataset object. This isn't working - currently it always pops up at the end of processing. Something wrong with the timing.
    # check_samples <- length(dataset$samples)
    #
    # if(check_samples < 1) {
    #     shinyalert::shinyalert(title = "Metadata Missing?",
    #                            text = "Looks like you have tried running MS-DAP before you have submitted your sample metadata. Please upload and submit the final metadata table and then try again.",
    #                            type = "error")
    # }

    #prep input for normalization methods
    norm_to_use <- input$msdap_norm_method

    #if mbprot selected, then add to this
    if(input$add_mbprot) {
        norm_to_use <- c(norm_to_use, "modebetween_protein")
    }

    #if norm to use is "none", replace with ""
    if(input$msdap_norm_method == "none") {
        norm_to_use <- ""
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
        dataset <- msdap::remove_proteins_by_name(dataset=dataset, irt_peptides = F, regular_expression = input$regexFilter)

    } else {

        #filter dataset according to user provided regex and/or gene/symbol filters
        dataset <- msdap::remove_proteins_by_name(dataset=dataset, irt_peptides = F, regular_expression = input$regexFilter, gene_symbols = input$geneFilter)

    }

    #save out modified metadata file that user uploaded with grouping and other information
    write.csv(final_metadata(), file = paste0(projectDir, "/input_data/", stringr::str_remove(input$modified_metadata$name, ".xlsx"), "_final_uploaded_metadata_file.csv"), row.names=F)

    #import the final metadata and add to the dataset object
    dataset <- msdap::import_sample_metadata(dataset, filename = input$modified_metadata$datapath)

    #if technical replicates indicated by user in UI, filter the peptides on FDR value, merge them together
    if(input$msdap_tech_reps) {

        dataset$peptides <- dataset$peptides |> dplyr::filter(detect)

        dataset = msdap::merge_replicate_samples(dataset, colname = "merge", minsample = 1, rescale_intensities = T)
    }

    #define contrasts
    dataset <- msdap::setup_contrasts(dataset, contrast_list = selected_contrast_pairs(), random_variables = NULL)

    #perform peptide-level filtering and normalization
    dataset = msdap::filter_dataset(dataset,
                                    filter_min_detect = input$msdap_min_detect_dea,
                                    filter_min_peptide_per_prot = input$msdap_min_pep_per_prot,
                                    norm_algorithm = norm_to_use,
                                    rollup_algorithm = input$pep2protMethod, #"maxlfq",
                                    by_group = T,
                                    all_group = F,
                                    by_contrast = F)

    #now add peptide intensities from by_group filtering to a column that will be recognized/used by DEA function
    dataset$peptides$intensity_norm = dataset$peptides$intensity_by_group


    #Perform DEA analysis
    #first, have to set up multithreading for msqrob
    cl <<- initialize_multiprocessing(n_thread = 10)
    dataset = msdap::dea(dataset,
                         qval_signif = input$dea_qval_thresh,
                         fc_signif = input$dea_log2sig_thresh,
                         dea_algorithm = input$de_algorithm,
                         rollup_algorithm = "maxlfq",
                         output_dir_for_eset = msdap_dir)

    #Now, if selected, perform the Differential Detection analysis
    if(input$msdap_do_dd) {

        dataset = differential_detect(dataset,
                                      min_peptides_observed = input$msdap_min_pept_dd,
                                      min_samples_observed = input$msdap_min_samples_dd,
                                      count_mode = "auto",
                                      rescale_counts_per_sample = TRUE)
    }

    #Export all tables and output files
    export_peptide_abundance_matrix(dataset, output_dir = msdap_dir)
    export_protein_abundance_matrix(dataset, output_dir = msdap_dir, rollup_algorithm = "maxlfq")
    export_statistical_results(dataset, output_dir = msdap_dir)

    data.table::fwrite(dataset$peptides, path_append_and_check(msdap_dir, "peptides.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$proteins, path_append_and_check(msdap_dir, "proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    data.table::fwrite(dataset$samples, path_append_and_check(msdap_dir, "samples.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    if(is_tibble(dataset$de_proteins) && nrow(dataset$de_proteins) > 0) {
        data.table::fwrite(dataset$de_proteins, path_append_and_check(msdap_dir, "de_proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    }
    if(is_tibble(dataset$dd_proteins) && nrow(dataset$dd_proteins) > 0) {
        data.table::fwrite(dataset$dd_proteins, path_append_and_check(msdap_dir, "dd_proteins.tsv.gz"), sep="\t", col.names = T, row.names = F, quote = F, na = "")
    }

    #Save out dataset object with qs
    qs::qsave(dataset, file = paste0(msdap_dir, "/dataset.qs"))

    #Save out captured input parameters as YAML file
    #capture all user input params and save out as YAML file
    contrasts_used = selected_contrast_pairs()
    print(contrasts_used)
    #contrasts_used <- c(dataset$contrasts[[1]]$label_contrast, dataset$contrasts[[2]]$label_contrast)
    params_yaml <- list(
        #capture UI input params
        project_id = if (!is.null(input$projectID)) input$projectID else NULL,
        search_tool = if (!is.null(input$analysisTool)) input$analysisTool else NULL,
        regex_filter = if (!is.null(input$regexFilter)) input$regexFilter else NULL,
        gene_filter = if (!is.null(input$geneFilter)) input$geneFilter else NULL,
        msdap_report = if (!is.null(input$data)) input$data else NULL,
        fasta_db = if (!is.null(input$database)) input$database else NULL,
        project_notes = if (!is.null(input$project_description)) input$project_description else NULL,
        pep_conf_thresh = if (!is.null(input$conf_thresh)) input$conf_thresh else NULL,
        retention_time_units = if (!is.null(input$spec_irt)) input$spec_irt else NULL,
        pd_peptide_collapse = if (!is.null(input$pep_collapse)) input$pep_collapse else NULL,
        pd_psm_per_precursor = if (!is.null(input$psm_per_precursor)) input$psm_per_precursor else NULL,
        sn_min_pep_per_prot = if (!is.null(input$msdap_min_pep_per_prot)) input$msdap_min_pep_per_prot else NULL,
        sn_technical_replicates = if (!is.null(input$msdap_tech_reps)) input$msdap_tech_reps else NULL,
        peptide_normalization_method = if (!is.null(input$msdap_norm_method)) input$msdap_norm_method else NULL,
        peptide_to_protein_rollup_method = if (!is.null(input$pep2protMethod)) input$pep2protMethod else NULL,
        use_mode_between_protein_norm = if (!is.null(input$add_mbprot)) input$add_mbprot else NULL,
        msdap_dea_qval_sig_thresh = if (!is.null(input$dea_qval_thresh)) input$dea_qval_thresh else NULL,
        msdap_dea_log2FC_sig_thesh = if (!is.null(input$dea_log2sig_thresh)) input$dea_log2sig_thresh else NULL,
        msdap_run_diff_detect = if (!is.null(input$msdap_do_dd)) input$msdap_do_dd else NULL,
        msdap_dea_min_peptide_detect_n_samples = if (!is.null(input$msdap_min_detect_dea)) input$msdap_min_detect_dea else NULL,
        msdap_dd_zscore_threshold = if (!is.null(input$dd_abs_zscore_thresh)) input$dd_abs_zscore_thresh else NULL,
        msdap_remove_dd_in_de = if (!is.null(input$filter_dd_in_de)) input$filter_dd_in_de else NULL,
        msdap_dd_min_peptide_per_prot = if (!is.null(input$msdap_min_pept_dd)) input$msdap_min_pept_dd else NULL,
        msdap_dd_min_samples_detected = if (!is.null(input$msdap_min_samples_dd)) input$msdap_min_samples_dd else NULL,
        msdap_dea_algorithms = if (!is.null(input$de_algorithm)) input$de_algorithm else NULL,
        msdap_contrasts = if (!is.null(contrasts_used)) contrasts_used else NULL
    )

    # Convert to YAML and save
    yaml_params <- as.yaml(params_yaml)
    writeLines(yaml_params, paste0(msdap_dir, "/msdap_params.yaml"))

    #once this finishes, then create the quickomics output

    print("Now extracting information for Quickomics Input Files...")

    # Quickomics Processing --------------------------------------------------------

    #create the quickomics output directory
    quickomics_dir <- paste0(projectDir, "/quickomics_files/")
    dir.create(quickomics_dir)
    print(paste0("Created output directory for Quickomics Files at ", quickomics_dir))

    # #process quickomics extraction for all contrasts
    # quickomics_expression_sets(msdap_output_directory = msdap_dir, msdap_data = dataset)

    #export quickomics files #NEED TO WRAP THESE EACH IN TRYCATCH WITH DISRUPTION OF HTML IF FAIL
    quickomics_expression(msdap_output_directory = msdap_dir, quickomics_output_directory = quickomics_dir)
    quickomics_de_dd_results(msdap_output_directory = msdap_dir, quickomics_output_directory = quickomics_dir, remove_dd_in_de = input$filter_dd_in_de, dd_zscore_thresh = input$dd_abs_zscore_thresh) #update this eventually once UI options are set up for input parameter passing
    quickomics_metadata(msdap_output_directory = msdap_dir, quickomics_output_directory = quickomics_dir)
    quickomics_gene_table(msdap_output_directory = msdap_dir, quickomics_output_directory = quickomics_dir)


    #now run the SIMR excel report function

    excel_report <- simr_excel_report(msdap_output_directory = msdap_dir, projName = input$projectID, msdap_data = dataset, sig.thresh = input$dea_qval_thresh, dd_zscore_thresh = input$dd_abs_zscore_thresh)

    #now, generate the gzipped output for download
    output$downloadResults<- renderUI({
        downloadButton(outputId = "downloadProteoVistaResults", label = "Download Excel Report of Results")
    })

    print("Creating output Excel report for download via download button....")
    # temp_zip <- tempfile(fileext = ".gz")
    # zip::zip(zipfile = temp_zip, files = project)

    #now assign the tempXL template to the download button
    output$downloadProteoVistaResults <- downloadHandler(
        filename = function() {
            paste0(input$projectID, "_ProteoVista_results.xlsx")
        },
        content = function(file) {
            openxlsx::write.xlsx(excel_report, file)

        }
    )

    # Capture session information
    session_info <- sessioninfo::session_info(to_file = paste0(msdap_dir, '/R_session_info.txt'))

    # # Extract relevant parts
    # pkg_info <- session_info$otherPkgs  # Load only non-base packages
    #
    # # Convert package info into a data frame
    # pkg_df <- data.frame(
    #     Package = names(pkg_info),
    #     Version = sapply(pkg_info, function(pkg) pkg$Version),
    #     stringsAsFactors = FALSE
    # )
    # write.csv(pkg_df, paste0(msdap_dir, '/R_session_info.csv'), row.names=F)

    #Generate PDF Report
    print("Now saving out the MS-DAP and ProteoVista PDF reports. This may take some time...")

    tryCatch({
        msdap::generate_pdf_report(dataset, output_dir = msdap_dir)
    }, error = function(e) {
        waiter_hide()

        shinyalert::shinyalert(
            title = "Error Generating MS-DAP Report PDF File",
            text = paste(e[["call"]][[2]]),
            type = "error"
        )
        return(NULL)
    })

    #generate_par_report(dataset, output_dir = msdap_dir)

    sink()
    waiter_hide()

    #Now, want to show pop-up window telling user that the processing is completed.
    shinyalert::shinyalert(title = "Dataset Successfully Processed with MS-DAP",
                           text = paste0("MS-DAP ouput files are located at ", msdap_dir, ". You can download as a zipped file using the 'Download Results' button."),
                           type = "success")


})

observeEvent(input$save_button, {
    # Capture all input values
    input_list <- reactiveValuesToList(input)

    # Convert to a data frame
    df <- data.frame(
        id = names(input_list),
        value = sapply(input_list, toString),  # handle vectors/lists
        stringsAsFactors = FALSE
    )

    # Save to CSV
    write.csv(df, file = "inputs_snapshot.csv", row.names = FALSE)

    message("✅ Inputs saved to inputs_snapshot.csv")
})
