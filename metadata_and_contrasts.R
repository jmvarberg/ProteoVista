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
        bsPopover("msdap_min_pep_per_prot", title = "DE/Protein Quant Minimum Peptide Filter", content = "This sets the minimum number of unique peptides a protein must have to be retained to be included in DEA analysis and rolled-up to protein level quantitation with MaxLFQ.", placement = "bottom", trigger = "hover"),
        # checkboxInput(
        #     inputId = "msdap_filter_by_contrast",
        #     label = "Filter By Contrast (If unchecked, filtering and normalizations applied globally; if checked, filtering and normalization applied on subset of relevant samples within each specified contrast pair)",
        #     value = TRUE,
        #     width = "100%"
        # ),
        checkboxInput(
            inputId = "msdap_tech_reps",
            label = "Combine technical replicates? If checked, must include information in 'merge' column of uploaded metadata file",
            value = FALSE,
            width = "100%"
        ),
        bsPopover("msdap_tech_reps", title = "Merge Technical Replicates", content = "If selected, precursors from each replicate are filtered for FDR, then averaged across replicates to derive a single precursor intensity value.", placement = "bottom", trigger = "hover"),
        selectInput(
            inputId = "msdap_norm_method",
            label = "Normalization Method",
            width = "100%",
            choices = list(
                "None" = "none",
                "Variance Stabilizing Normalization (VSN)" = "vsn",
                "Variation Within, Mode Between" = "vwmb",
                "Mode Within, Mode Between" = "mwmb",
                "Median" = "median",
                "Loess" = "loess",
                "Robust Linear Regression (RLR)" = "rlr",
                "MS-EmpiRe" = "msempire"
                # "Mode Between (Peptide)" = "mb_pep",
                # "Mode Between (Protein)" = "mb_prot"
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
        bsPopover("dea_qval_thresh", title = "DE FDR Threshold", content = "This selection only influences the hits that are classified as significant in MS-DAP PDF report and SIMR Excell ProteoVista Summary Report. All values are returned in DE results exported as CSV for Quickomics.", placement = "bottom", trigger = "hover"),
        numericInput(
            inputId = "dea_log2sig_thresh",
            label = "Absolute log2-Fold Change Significance Threshold",
            value = 0,
            width = "100%",
            min = 0
        ),
        bsPopover("dea_log2sig_thresh", title = "DE Significance Threshold", content = "This selection only influences the hits that are classified as significant in MS-DAP PDF report and SIMR Excell ProteoVista Summary Report. All values are returned in DE results exported as CSV for Quickomics. Defalut value returns all hits based solely on FDR threshold (no logFC requirements).", placement = "bottom", trigger = "hover"),
        checkboxInput(
            inputId = "msdap_do_dd",
            label = "Do you want to perform Differential Detection analysis?",
            value = T,
            width = "100%"
        ),
        bsPopover("msdap_do_dd", title = "Differential Detection (Optional)", content = "MS-DAP uses Differential Detection to identify proteins that are unique or significantly enriched in one group vs. the other that would otherwise be missed by DEA stats testing.", placement = "bottom", trigger = "hover")
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
        bsPopover("msdap_min_detect_dea", title = "DEA Peptide Filter", content = "Peptides must be detected in at least N samples in at least one group to be retained for protein-level quantitation using MaxLFQ, and for DE testing. Minimum recommended value is 3, can be increased if there are more replicates within each group.", placement = "bottom", trigger = "hover")
    ) #close tagList

}) #close renderUI()

#UI element to specify peptide filtering criteria for Differential Detection analysis (dd)
output$mdsap_dd_filt_params <- renderUI({

    #diffdetect_min_peptides_observed = 2, #minimum number of peptides for a protein to be included in diff. detection analyses (protein level filter)
    #diffdetect_min_samples_observed = 3, #minimum number of samples within a condition that a peptide must be detected in to be kept for diff. detection (peptide level filter)

    tagList(
        numericInput(
            inputId = "dd_abs_zscore_thresh",
            label = "Absolute z-score filter for Differential Detection Hits",
            value = 4,
            width = "100%",
            min = 0
        ),
        bsPopover("dd_abs_zscore_thresh", title = "Differential Detection Z-score Filtering", content = "This value sets the minimum |z-score| value threshold. Any DD hits with values less than this are not included in DEA/DD results output in ProteoVista Summary Excel document or for Quickomics Files. Full DD lists are available in msdap_output directory.", placement = "bottom", trigger = "hover"),
        checkboxInput(
            inputId = "filter_dd_in_de",
            label = "Remove any hits in DD that have results in DE?",
            value = T,
            width = "100%"
        ),
        bsPopover("filter_dd_in_de", title = "Differential Detection Filtering (Optional)", content = "Occasionally, proteins that are used for DEA also show up as significant hits in Differential Detection. This optionally removes any DD hits that have results in DE to avoid two logFC values for a single protein/contrast.", placement = "bottom", trigger = "hover"),
        numericInput(
            inputId = "msdap_min_pept_dd",
            label = "DD Protein Filter:",
            value = 2,
            width = "100%"
        ),
        bsPopover("msdap_min_pept_dd", title = "Differntial Detection Filtering", content = "Minimum number of peptides for a protein to pass filtering rules (i.e. otherwise, no z-score is computed)", placement = "bottom", trigger = "hover"),
        numericInput(
            inputId = "msdap_min_samples_dd",
            label = "DD Samples Filter:",
            value = 3,
            width = "100%"
        ),
        bsPopover("msdap_min_samples_dd", title = "Differential Detection Filtering", content = "Minimum number of samples where a protein should be observed with at least the specified minimum number of peptides (in either group) when comparing a contrast of group A vs B", placement = "bottom", trigger = "hover")

    ) #close tagList

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
        choices = list("limma" = "ebayes",
                       "DEqMS" = "deqms",
                       "msEmpiRe" = "msempire",
                       "MSqRob" = "msqrob"),
        multiple=TRUE,
        selected = c("ebayes", "deqms", "msempire", "msqrob"),
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



