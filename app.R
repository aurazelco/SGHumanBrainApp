# Author: Aura Zelco
# Brief description:
  # This script is used to build a web app to explore the sex differences in cell types from single-cell/nucleus human brain samples at different ages.
  # The web app allows the user to select the threshold for fold change (FC) and p-values, to explore how the results differ if different thresholds are applied
# Brief procedure:
  # 1. Imports the necessary libraries and custom functions
  # 2. Imports the unfiltered DEGs
  # 3. Sets up the user interface (UI), with multiple tabs and widgets
  # 4. Sets up the server, where the calculations and plots are going ot be generated
  # 5. Using default or user-defined parameters for FC and p-values, plots the results of the analyses in the UI

#---------------------------------------------------------------------------------------------------

############################################################################################

# 1. Imports the necessary libraries and custom functions

library(shiny)
library(shinydashboard)

# General summary pages
source("custom_scripts/About.R")
source("custom_scripts/Source.R")

# Single cell datasets and DEGs
source("custom_scripts/DsInfo.R")
source("custom_scripts/DEGs.R")

# Bulk RNA-seq datasets and DEGs
source('custom_scripts/Healthy_reg.R')
source('custom_scripts/Healthy_allreg2.R')


############################################################################################


############################################################################################

# 3. UI

ui <- dashboardPage(
  skin = "purple",
  # title of web app
  dashboardHeader(title = "SGHumanBrain", titleWidth = 300),
  
  # tabs
  dashboardSidebar(
    width = 300,
    sidebarMenu(
    menuItem("Summary", tabName = "summary", icon = icon("dashboard")),
    menuItem("ScRNA-seq - Datasets information", tabName = "ds_info", icon = icon("folder-open", lib="glyphicon")),
    menuItem("ScRNA-seq - DEGs Analysis", tabName = "degs", icon = icon("stats", lib="glyphicon")),
    menuItem("Bulk RNA-seq - Specific brain region", tabName = "Healthy1", icon = icon("stats", lib="glyphicon")),
    menuItem("Bulk RNA-seq - Across brain regions", tabName = "Healthy2", icon = icon("stats", lib="glyphicon")),
    menuItem("Open source code", tabName = "source_code", icon = icon("cog", lib="glyphicon"))
  )),
  
  dashboardBody(
      tabItems(
        tabItem(tabName = "summary",
                h2("SGHumanBrain"),
                AboutUI("Summary")
        ),
        
        tabItem(tabName = "ds_info",
                h2("ScRNA-seq - Datasets information"),
                DsUI("Datasets information")
        ),
        
        tabItem(tabName = "degs",
                h2("ScRNA-seq - DEGs Analysis"),
                DEGsUI("DEGs")
        ),
        
        tabItem(tabName = "Healthy1",
                h2("Sex-biased genes rank for healthy brain samples"),
                Healthy_UI('Healthy1')
        ),
        
        tabItem(tabName = "Healthy2",
                h2("Sex-biased genes rank across brain regions for healthy brain samples"),
                HealthyAllUI('Healthy2')
        ),
        
        tabItem(tabName = "source_code",
                h2("Open source code"),
                SourceUI("Open source code")
        )
        
      )
    )
)



############################################################################################

# 4. Server

server <- function(input, output, session) {
  DEGsServer("DEGs")
  Healthy_Server('Healthy1')
  HealthyAllServer('Healthy2')
}


############################################################################################

# 5. Runs app

shinyApp(ui, server)
