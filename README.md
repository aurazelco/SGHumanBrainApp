# SGHumanBrainApp

This repository contains the code for the web app [SGHumanBrainApp](link) and the files used to build the app. The link for the SGHumanBrainApp web app is https://joshiapps.cbu.uib.no/GenderRank_app. 


## Table of Contents

The repo contains the following content:
- [app.R](app.R)
- [data](data/) folder
- [custom_scripts](custom_scripts/) folder
- [www](www/) folder

----------------------------------------------------------------------------------------------------------

## app.R

Main R script to build the web application using Shiny and shinydashboard. It organizes the information in multiple tabs, and gives the user the possibility to change the adjusted p-value and fold change (FC) threhsold, to see how the different plots vary accordingly. 

## Data folder

### Single-cell RNA-seq
In this folder, all the files needed to build the plots are present. The [Unfiltered_DEGs](data/Unfiltered_DEGs/) contains all the DEGs as calculated from Seurat, organized by dataset and cell type (original annotation). The [Filtered_DEGs](data/Filtered_DEGs/) instead has the same DEGs, but filtered by several combinations of adjusted p-value and fold change thresholds, again organized by datasets and cell types. 

The [extra_files](data/extra_files/) contains the files needed to generate the output plots. 

### Bulk RNA-seq
Some of the results from the bulk RNA-seq analysis are already calculated, while other plots are generated according to the user's input. The two tabs in the app have two corresponding data folder, namely [Healthy](data/Healthy/) and [pre_result_allreg](data/pre_result_allreg/). [Healthy](data/Healthy/) contains the DEGs lists for each study used in each brain region, and they can be downloaded. [pre_result_allreg](data/pre_result_allreg/) instead contains the plots and data from the comparison across brain regions. 

## Custom_scripts

The scripts are used to generate all the different files and plots, and to build the Rshiny app. 

Brief description of the custom scripts:
- [Generate_output_files.R](custom_scripts/Generate_output_files.R): script to generate the files and plots used and displayed in the web application
- [Generate_output_files_func.R Source.R](custom_scripts/Generate_output_files_func.R): custom functions, sourced in Generate_output_files.R
- [About.R](custom_scripts/About.R): "Summary" tab user interface and server code
- [DEGs.R](custom_scripts/DEGs.R): "ScRNA-seq - DEGs Analysis" tab user interface and server code
- [DsInfo.R](custom_scripts/DsInfo.R): "ScRNA-seq - Dataset information" tab user interface and server code
- [Healthy_reg.R](custom_scripts/Healthy_reg.R): "Bulk RNA-seq - Specific brain region" tab user interface and server code
- [Healthy_allreg2.R](custom_scripts/Healthy_allreg2.R): "Bulk RNA-seq - Across brain regions" tab user interface and server code
- [Source.R](custom_scripts/Source.R): "Source" tab user interface and server code
- [enrichr.R](custom_scripts/enrichr.R), [Fun_Rank_disease.R](custom_scripts/Fun_Rank_disease.R) and [global.R](custom_scripts/global.R): scripts containing functions used in the bulk RNA-seq tabs

## www

### Single-cell RNA-seq
The [www](www/) folder contains all plots and images, mainly as PNG, which are displayed for the user to explore the data. The images in [Plots](wwww/Plots/) were either generated using [Generate_output_files.R](custom_scripts/Generate_output_files.R), while the number of samples etc were generated in the main analysis, which can be found [here](https://github.com/aurazelco/MSc_Thesis_scripts/tree/main/Suppl_files). 

### Bulk RNA-seq
Images with information about the datasets per region can be found [here](wwww/bulk/). 

