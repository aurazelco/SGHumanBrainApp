# HumanBrainSexSingleCell

This repository contains the code for the web app [HumanBrainSexSingleCell](link) and the files used to build the app. 

The repo contains the following content:
- [app.R](app.R)
- [data](data/) folder
- [custom_scripts](custom_scripts/) folder
- [www](www/) folder

## app.R

Main R script to build the web application using Shiny and shinydashboard. It organizes the information in multiple tabs, and gives theuser the possibility to change the adjusted p-value and fold change (FC) threhsold, to see how the different plots vary accordingly. 

## Data folder

In this folder, all the files needed to build the plots are present. The [Unfiltered_DEGs](data/Unfiltered_DEGs/) contains all the DEGs as calculated from Seurat, organized by dataset and cell type (original annotation). The [Filtered_DEGs](data/Filtered_DEGs/) instead has the same DEGs, but filtered by several combinations of adjusted p-value and fold change thresholds, again organized by datasets and cell types. 

The [extra_files](data/extra_files/) contains the files needed to compare or analyze the DEGs. 


## Custom_scripts

The scripts are used to generate all the different files and plots, and to build the Rshiny app. 

Brief description of the custom scripts:
- [Generate_output_files.R](custom_scripts/Generate_output_files.R): script to generate the files and plots used and displayed in the web application
- [Generate_output_files_func.R Source.R](custom_scripts/Generate_output_files_func.R): custom functions, sourced in Generate_output_files.R
- [About.R](custom_scripts/About.R): "Summary" tab user interface and server code
- [DEGs.R](custom_scripts/DEGs.R): "Differentially expressed genes (DEGs)" tab user interface and server code
- [DsInfo.R](custom_scripts/DsInfo.R): "Dataset information" tab user interface and server code
- [Source.R](custom_scripts/Source.R): "Source" tab user interface and server code


## www

The [www](www/) folder contains all plots and images, mainly as PNG, which are displayed for the user to explore the data. The images were either generated using [Generate_output_files.R](custom_scripts/Generate_output_files.R), while the number of samples etc were generated in the main analysis, which can be found [here](https://github.com/aurazelco/MSc_Thesis_scripts/tree/main/Suppl_files). 