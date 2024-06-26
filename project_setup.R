#Step 1: Project Setup Server Side

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
msdap_dataset <- reactiveValues(dataset = NULL)

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
        msdap_dataset_2 <<- dataset
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

    msdap_dataset$dataset <- dataset

})

