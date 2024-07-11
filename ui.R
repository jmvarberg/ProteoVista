navbarPage(
    # Tab 1: Project Setup and Data Input -------------------------------------
    "ProteoVista",
    tabPanel("Step 1: Input Files", value = "setup",
             autoWaiter(),
             sidebarLayout(
                 sidebarPanel(
                     h3("Project Setup"),
                     textInput(inputId = "projectID",
                               label = "Project Name (PROT)",
                               placeholder = "PROT-XXX..."
                     ),
                     uiOutput(outputId = "primaryAnalysis"),
                     uiOutput(outputId = "inputData"),
                     uiOutput(outputId = "inputDatabase"),
                     uiOutput(outputId = "inputNotes"),
                     br(),
                     br(),
                     h3("Select MS-DAP Import Parameters"),
                     uiOutput(outputId = "mdsap_input_params"),
                     br(),
                     actionButton(inputId = "submit_input", label = "Submit Input Files"),
                     br(),
                     h3("MS-DAP Summary:"),
                     uiOutput(outputId = "msdap_summary")
                 ),
                 mainPanel(
                     uiOutput(outputId = "tab_1_md")
                 )
             )),
    # Tab 2: Parameters and Contrast Definitions ------------------------------
    tabPanel("Step 2: Set Up Analysis", value = "params",
             sidebarLayout(
                 sidebarPanel(
                     #Upload the modified sample metadata from projdir
                     h3("Select Metadata"),
                     uiOutput(outputId = "inputMetadata"),
                     actionButton(inputId = "submit_metadata", "Submit Metadata File"),
                     br(),
                     h3("Select MS-DAP Processing Parameters"),
                     uiOutput(outputId = "mdsap_run_params"),
                     br(),
                     h3("Specify Contrasts for Differential Testing"),
                     uiOutput(outputId = "contrastMethod"),
                     #uiOutput(outputId = "referenceSelect"),
                     conditionalPanel("input.contrast_method=='contrast_ref'", uiOutput(outputId = "referenceSelect")),
                     actionButton(inputId = "submit_msdap", "Run MS-DAP")
                 ),
                 mainPanel(
                     h1("Sample Metadata Table"),
                     br(),
                     DT::dataTableOutput(outputId = "finalMeta"),
                     br(),
                     fluidRow(
                         column(6,
                                h3("Groups Detected: "),
                                htmlOutput(outputId = "groupsCheck")),
                         column(6,
                                h3("Selected Comparisons To Be Used for Differential Testing"),
                                htmlOutput(outputId = "finalContrasts"))
                     ),
                     br(),
                     conditionalPanel("input.contrast_method=='constrast_custom'", uiOutput(outputId = "contrastSelection")),
                     br(),
                 )
             )),
    # Tab 3: Summary QC Plots -------------------------------------------------
    tabPanel("Step 3: Summary Plots", value = "summary_plots"), #to build here: number of peptides and proteins per group; percent data complete; CV distributions; data points per peak(?); PCA plot; Heatmap top X variable features; Table with # differentially expressed per contrast? DE results table?


    # Tab 4: Comparison of Two ProteoVista Results ----------------------------

    tabPanel("Compare Two Completed Results", value = "compare_results"),


    # Help Section ------------------------------------------------------------

    navbarMenu("About/Help",
               tabPanel("MS-DAP", value = "help_msdap"),
               tabPanel("Normalization Methods", value = "help_norm"),
               tabPanel("Differential Expression Methods", value = "help_dea",
                        fluidPage(
                            uiOutput(outputId = "dea_md")
                        )
               ),
               tabPanel("Quickomics", value = "help_quick")
    )
)

