## Setup - source after installation to setup necessary output files and directories

#create output folder directory if it doesn't exist
if(!dir.exists("ProteoVista_output")) {
    message("Initiallizing ProteoVista Output Subfolder...")
    dir.create("./ProteoVista_output")
}

#install any missing packages based on renv dependencies in renv.lock file
if (requireNamespace("renv", quietly = TRUE)) {
    renv::restore(prompt = FALSE)
} else {
    message("Package 'renv' is not installed. Please run install.packages('renv') and then re-run source(setup.R)"
}
