# ProteoVista

ProteoVista is a Shiny dashboard for proteomics analysis. It provides a GUI for setting up and processing mass spectrometry data with functions provided by the R package [MS-DAP](https://github.com/ftwkoopmans/msdap). 
After processing by MS-DAP, output files are formatted for uploading to Quickomics Shiny dashboard for interactive data visualization and exploration.

ProteoVista is currently in development, with a beta version available for processing Spectronaut DIA datasets. Support for additional data types is currently being built and tested.

## Installation

Stable beta releases are made available at the following link:

https://github.com/jmvarberg/ProteoVista/releases/

Please use the most recent release version and download the source code files. Once downloaded, unzip and save the folder in your desired location.

After installation, and before the first use, please open R from within the ProteoVista directory and use the following command to setup the directory for saving outputs:

```{r}
source("setup.R")
```

This will create an output directory where results files will be stored. It will also use the `renv` package function `renv::restore()` to install all required package dependencies.

Once installation of package dependencies is complete, you can launch the app with the command: 

```{r}
shiny::runApp()
```

This should launch ProteoVista in a browser window.
