#Step 1: Project Setup Server Side

# Specification of Sidebar UI elements --------------------------------------------

#UI for Selecting Primary Analysis Tool Used
output$primaryAnalysis <- renderUI({

    selectInput(
        inputId = "analysisTool",
        label = "Primary Analysis Software",
        choices = list("Spectronaut" = "spectronaut"),
                       #"FragPipe/IonQuant (DDA)" = "fp_ionquant",
                       #"FragPipe/DIA-NN" = "fp_diann",
                       #"ProteomeDiscoverer" = "pd"),
        selected = "spectronaut",
        width = "100%"
    )
})

#UI to dynamically change the presentation of the Upload File selection based on software selection
output$inputData <- renderUI({

    req(input$analysisTool)

    if(input$analysisTool == "fp_ionquant") {

        #multifile upload - need combined_protein.tsv, MSstats.csv, and all psm.tsv files in samples subdirectories.
        fileInput(
            inputId = "data",
            label = "Upload FragPipe Files: 1) combined_protein.tsv, 2) MSstats.csv, and 3) all psm.tsv files - these are located in sample subdirectories, need one per sample.",
            multiple = T,
            width = "100%",
            accept = c(".tsv", ".csv")
        )
        # textInput(
        #     inputId = "data",
        #     label = "Path to FragPipe Output Directory",
        #     placeholder = "/n/proteomics/..."
        # )

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
            bsTooltip("conf_thresh", "This is based on probability scores calculated by PeptideProphet and stored as PeptideProphet.Probability value in each psm.tsv report", "bottom", options = list(container = "body")),
            #bsTooltip("conf_thresh", "Confidence score threshold at which a peptide is considered 'identified', should be a numeric value between 0 and 1 (peptides with a '1 - PeptideProphet.Probability' value in psm.tsv that is <= this threshold are classified as 'detected')", "right", options = list(container = "body")),
            selectInput(inputId = "pep_collapse",
                        label = "Select method to combine/collapse peptides (use modified or plain peptide sequences)?",
                        choices = list("Modified" = "sequence_modified", "Plain" = "sequence_plain"),
                        multiple = FALSE,
                        selected = "mod"),
            bsTooltip("pep_collapse", "If multiple data points are available for a peptide in a sample, at what level should these be combined?", "bottom", options = list(container = "body"))
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
            bsTooltip("conf_thresh", "This is based on precursor-level Q-value from Spectronaut", "right", options = list(container = "body")),
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
            bsTooltip("pep_collapse", "If multiple data points are available for a peptide in a sample, at what level should these be combined?", "bottom", options = list(container = "body")),
            selectInput(inputId = "psm_per_precursor",
                        label = "Select how to to roll-up PSM quantitation to precursors.",
                        choices = list("Sum" = "", "Intensity" = "intensity", "Confidence" = "confidence"),
                        multiple = FALSE,
                        selected = "Sum"),
            bsTooltip("psm_per_precursor", "If ProteomeDiscoverer performed peak integration and reports the same (redundant) peak intensity for each PSM of the same precursor, we suggest to use Intensity. Use Sum to use the sum of all PSM intensity values per precursor*sample (default). Intensity to select the highest intensity value of PSMs above confidence threshold. Confidence to select the intensity value from the PSM with best/lowest confidence value", placement = "bottom", options = list(container = "body")),
            shinyWidgets::switchInput(inputId = "remove_low_conf", label = "Remove low-confidence peptides (as defined by PD)?", value = TRUE, width = "auto", size = "small", labelWidth = "250px", handleWidth = "50px", inline = TRUE)
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

#UI to print dataset summary text
output$msdap_summary <- renderUI({

    #if user has not submitted yet, print "awaiting dataset upload"
    if(input$submit_input == 0) {
        renderText("Waiting for dataset to be uploaded.")
    } else (
        renderText(capture.output(msdap::print_dataset_summary(msdap_dataset$dataset)))
    )
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
    #outdir <- "/n/proteomics/washburn/Joe/Utilities/ProteoVista_output/"
    outdir <- "./ProteoVista_output/"
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

    #save out project notes as text file
    proj_notes <- input$project_description
    sink(paste0(projdir, "/input_data/", input$projectID, "_project_description_notes.txt"))
    cat(proj_notes)
    sink()

    #copy input files into the project directory - moved data copy into individual if statements to work with fragpipe_ionquant requirement for folder path instead of individual files
    #file.copy(input$data$datapath, paste0(projdir, "/input_data/", input$data$name))
    file.copy(input$database$datapath, paste0(projdir, "/input_data/", input$database$name))

    if(input$analysisTool == "spectronaut") {
        waiter_show(html = data_ingest_waiting)

        file.copy(input$data$datapath, paste0(projdir, "/input_data/", input$data$name))

        # print(paste0("conf_thresh: ", input$conf_thresh))
        # print(paste0("iRT: ", input$spec_irt))
        # print(paste0(projdir, "/", input$data$name))

        #Step 1: Import the Spectronaut Report file that was uploaded by reading from the project directory
        dataset <- tryCatch({
                msdap::import_dataset_spectronaut(filename = paste0(projdir, "/input_data/", input$data$name),
                                                  confidence_threshold = input$conf_thresh
                )
        }, error = function(e) {
            waiter_hide()

            shinyalert::shinyalert(
                title = "Error Importing Spectronaut File",
                text = paste(e[["call"]][[2]]),
                type = "error"
            )
            return(NULL)
        })

        if (is.null(dataset)) return()

        #Step 2: Import FASTA database file and add to the dataset
        dataset <- tryCatch({
            msdap::import_fasta(dataset,
                                files = paste0(projdir, "/input_data/", input$database$name)
                                )
        }, error = function(e) {
            waiter_hide()

            shinyalert::shinyalert(
                title = "Error Importing FASTA File",
                text = paste(e[["call"]][[2]]),
                type = "error"
            )
            return(NULL)
        }
        )

        #if there is an error return the dataset as a null object. User can re-upload and try again by re-clicking submit button
        if (is.null(dataset)) return()

        #write out sample metadata template to projdir
        msdap::write_template_for_sample_metadata(dataset, filename = paste0(projdir, "/input_data/", input$projectID, "_sample_metadata_table.xlsx"), overwrite = T)
        #msdap_dataset_2 <<- dataset
        waiter_hide()

        #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
        shinyalert::shinyalert(title = "Dataset Uploaded and Metadata Template Generated",
                               text = paste0("Please go to Step 2, download the metadata template file, and upload after completing group assignments."),
                               type = "success")

    } else if(input$analysisTool == "fp_diann") {

        waiter_show(html = data_ingest_waiting)

        file.copy(input$data$datapath, paste0(projdir, "/input_data/", input$data$name))

        #import the Spectronaut Report file that was uploaded by reading from the project directory
        dataset <- msdap::import_dataset_diann(filename = paste0(projdir, "/", input$data$name))

        #add the fasta database to the dataset
        dataset <- msdap::import_fasta(dataset, files = paste0(projdir, "/", input$database$name))

        #write out sample metadata template to projdir
        msdap::write_template_for_sample_metadata(dataset, filename = paste0(projdir, "/", input$projectID, "_sample_metadata_table.xlsx"), overwrite = T)

        waiter_hide()

        #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
        shinyalert::shinyalert(title = "Dataset Uploaded and Metadata Template Generated",
                               text = paste0("Please go to Step 2, download the metadata template file, and upload after completing group assignments."),
                               type = "success")

    } else if(input$analysisTool == "pd") {

        waiter_show(html = data_ingest_waiting)

        file.copy(input$data$datapath, paste0(projdir, "/input_data/", input$data$name))

        #Step 1: import the ProteomeDiscoverer Report file that was uploaded by reading from the project directory
        dataset <- tryCatch({
            msdap::import_dataset_proteomediscoverer_txt(filename = paste0(projdir, "/", input$data$name),
                                                         collapse_peptide_by = input$pep_collapse,
                                                         confidence_threshold = input$conf_thresh,
                                                         remove_lowconf = input$remove_low_conf,
                                                         one_psm_per_precursor = input$psm_per_precursor
            )
        }, error = function(e) {
            waiter_hide()

            shinyalert::shinyalert(
                title = "Error Importing ProteomeDiscoverer PSM.txt File",
                text = paste(e[["call"]][[2]]),
                type = "error"
            )
            return(NULL)
        }
        )

        if (is.null(dataset)) return()

        #Step 2: Import FASTA database file and add to the dataset
        dataset <- tryCatch({
            msdap::import_fasta(dataset,
                                files = paste0(projdir, "/input_data/", input$database$name)
            )
        }, error = function(e) {
            waiter_hide()

            shinyalert::shinyalert(
                title = "Error Importing FASTA File",
                text = paste(e[["call"]][[2]]),
                type = "error"
            )
            return(NULL)
        }
        )

        #if there is an error return the dataset as a null object. User can re-upload and try again by re-clicking submit button
        if (is.null(dataset)) return()

        #write out sample metadata template to projdir
        msdap::write_template_for_sample_metadata(dataset, filename = paste0(projdir, "/", input$projectID, "_sample_metadata_table.xlsx"), overwrite = T)

        waiter_hide()

        #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
        shinyalert::shinyalert(title = "Dataset Uploaded and Metadata Template Generated",
                               text = paste0("Please go to Step 2, download the metadata template file, and upload after completing group assignments."),
                               type = "success")
    } else if(input$analysisTool == "fp_ionquant") {

        waiter_show(html = data_ingest_waiting)

        # #Create temp directory to store uploaded files
        # temp_fp_dir <- file.path(tempdir(), "uploaded_files")
        # if (!dir.exists(temp_fp_dir)) {
        #     dir.create(temp_fp_dir)
        # }

        #store all of the uploaded files in an object
        uploaded_files <- input$data
        print(uploaded_files)

        input_dir <- paste0(projdir, "/input_data")
        print(input_dir)

        # Copy each uploaded file to the projdir/input_data directory
        #For psm.tsvs - msdap import helper function cannot detect if they have prefix. it's also built to look for them in all subdirectories of the path you provide.
        #So, to get this to work, specifically for psm.tsv files, have it extract the prefix, create a subdir in the input_data directory with that name, then copy that respective psm.tsv file inside it with the prefix removed.
        #all other non psm.tsv files just go straight into the input_data directory.
        #Note, this requires users to manually modify the psm.tsv files in the FragPipe output folders to have the correct name for the prefix - no great work around as temp dir from Shiny upload not aware of parent folder on client side.
        for (i in seq_len(nrow(uploaded_files))) {
            filename <- uploaded_files$name[i]
            from_path <- uploaded_files$datapath[i]

            if (grepl("psm\\.tsv$", filename)) {
                # Extract prefix
                prefix <- sub("[-_]?psm\\.tsv$", "", filename)

                # Create subfolder
                subfolder <- file.path(input_dir, prefix)
                if (!dir.exists(subfolder)) dir.create(subfolder, recursive = TRUE)

                # Define target path
                to_path <- file.path(subfolder, "psm.tsv")

                # Read, modify, and write the data
                df <- read.delim(from_path, sep = "\t", stringsAsFactors = FALSE)

                # Modify a column here — for example, add 1 to intensity
                if ("Spectrum.File" %in% colnames(df)) {
                    df$Spectrum.File <- stringr::str_match(df$Spectrum.File, "interact-(.*?)\\.pep\\.xml")[, 2]

                }

                # Write modified file

                write.table(df, file = to_path, sep = "\t", quote = FALSE, row.names = FALSE)

            } else {
                # Copy other files unchanged
                to_path <- file.path(input_dir, filename)
                file.copy(from = from_path, to = to_path, overwrite = TRUE)
            }
        }

        #now, pass the input_data dir with all of the uploaded files to the msdap import function
        print("Importing FragPipe Files")
        dataset <- tryCatch({
            msdap::import_dataset_fragpipe_ionquant(path = paste0(projdir, "/input_data"),
                                                    acquisition_mode = input$fp_mode,
                                                    confidence_threshold = input$conf_thresh,
                                                    collapse_peptide_by = input$pep_collapse
            )
        }, error = function(e) {
            waiter_hide()

            shinyalert::shinyalert(
                title = "Error Importing ProteomeDiscoverer PSM.txt File",
                text = paste(e[["call"]][[2]]),
                type = "error"
            )
            return(NULL)
        }
        )

        if (is.null(dataset)) return()

        #Step 2: Import FASTA database file and add to the dataset
        dataset <- tryCatch({
            msdap::import_fasta(dataset,
                                files = paste0(projdir, "/input_data/", input$database$name)
            )
        }, error = function(e) {
            waiter_hide()

            shinyalert::shinyalert(
                title = "Error Importing FASTA File",
                text = paste(e[["call"]][[2]]),
                type = "error"
            )
            return(NULL)
        }
        )

        if (is.null(dataset)) return()

        #write out sample metadata template to projdir
        print("Writing metadata template file")
        msdap::write_template_for_sample_metadata(dataset, filename = paste0(input_dir, "/", input$projectID, "_sample_metadata_table.xlsx"), overwrite = T)

        waiter_hide()

        #Now, want to show pop-up window telling user to find the sample metadata file, edit, and save, then move to Step 2 tab.
        shinyalert::shinyalert(title = "Dataset Uploaded and Metadata Template Generated",
                               text = paste0("Please go to Step 2, download the metadata template file, and upload after completing group assignments."),
                               type = "success")

    }

    #make the download button appear for the metadata template
    output$metaTemplate <- renderUI({
        downloadButton(outputId = "downloadMetaTemplate", label = "Download Metadata Template File")
    })

    #now assign the tempXL template to the download button
    output$downloadMetaTemplate <- downloadHandler(
        filename = function() {
            paste0(input$projectID, "_sample_metadata_table.xlsx")
        },
        content <- function(file) {
            #now save the metadata template to a temp file
            tempXL <- tempfile(fileext = ".xlsx")
            msdap::write_template_for_sample_metadata(dataset, filename = tempXL, overwrite=T)
            md <- readxl::read_xlsx(path = tempXL)
            openxlsx::write.xlsx(md, file)
        }
    )

    msdap_dataset$dataset <- dataset

})

