function(input, output, session) {
    source('./global.R', local = TRUE)
    source('./project_setup.R', local = T)
    source('./metadata_and_contrasts.R', local = T)
    source('./msdap_processing.R', local=T)
    source('./help_sections.R', local=T)
    source('./visualizations.R', local=T)
    source('./quickomics_processing.R', local=T)
    #source('./combined_script.R', local = TRUE)
}
