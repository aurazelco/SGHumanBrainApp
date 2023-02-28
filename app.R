# Author: Aura Zelco
# Brief description:
  # This script is used to build a web app to explore the sex differences in cell types from single-cell/nucleus human brain samples at different ages.
  # The web app allows the user to select the threshold for fold change (FC) and p-values, to explore how the results differ if different thresholds are applied
# Brief procedure:
  # 1. Imports the necessary libraries and custom functions
  # 2. Imports the unfitlered DEGs
  # 3. Sets up the user interface (UI), with multiple tabs and widgets
  # 4. Sets up the server, where the calculations and plots are going ot be generated
  # 5. Using default or user-defined parameters for FC and p-values, plots the results of the analyses in the UI

#---------------------------------------------------------------------------------------------------

############################################################################################

# 1. Imports the necessary libraries and custom functions

library(shiny)
library(shinydashboard)
source("custom_scripts/ImportDEGs.R")
source("custom_scripts/About.R")
source("custom_scripts/Source.R")
source("custom_scripts/DsInfo.R")
source("custom_scripts/DEGs.R")

############################################################################################

# 2. Imports the unfiltered DEGs and generates another list, plus defines common annotation and group order to be followed in plots

# manually decided how to combine the sub-celltypes
unified_annotation <- c("CXCL14 IN" = "Interneurons",
                        "EC" = "Endothelial cells",
                        "fibrous astrocyte"  = "Astrocytes",
                        "L2_3 EN" = "Excitatory neurons", 
                        "L4 EN" = "Excitatory neurons",
                        "L5 EN" = "Excitatory neurons",
                        "L5_6 EN" = "Excitatory neurons",
                        "L5b EN" = "Excitatory neurons",
                        "L6 EN" = "Excitatory neurons",                
                        "microglia" = "Microglia", 
                        "Oligodendrocyte" =  "Oligodendrocytes",      
                        "OPC" = "OPCs",                  
                        "PLCH1 L4_5 EN" = "Excitatory neurons", 
                        "protoplasmic astrocyte" = "Astrocytes",
                        "PVALB IN"  = "Interneurons",            
                        "pyramidal neuron"  = "Excitatory neurons",
                        "SST IN" = "Interneurons",   
                        "SV2C IN"  = "Interneurons",   
                        "TSHZ2 L4_5 EN" = "Excitatory neurons",  
                        "VIP IN" = "Interneurons",
                        "Mesenchymal" = "Mesenchymal",      
                        "Neuroepithelial" =     "Neuroepithelial",
                        "Neuronal" = "Neurons",            
                        "Other"    = "Other",                
                        "Radial Glial"     = "Radial Glia",       
                        "Astrocytes" = "Astrocytes",        
                        "Excitatory neurons"  = "Excitatory neurons",
                        "Interneurons"   = "Interneurons",     
                        "Microglia"  = "Microglia",         
                        "Oligodendrocytes" = "Oligodendrocytes",
                        "OPCs" = "OPCs",            
                        "Unknown" = "Unknown",           
                        "Vascular cells" = "Vascular cells",     
                        "Dorsal progenitors"  = "Dorsal progenitors" ,   
                        "Ventral progenitors" = "Ventral progenitors")
names(unified_annotation) <- tolower(names(unified_annotation))

# defines the order in which to organize the presence heatmaps, so the groups are in developmental order, with the last groups as diseases
groups_order <- c("Velmeshev_2nd_trimester",        
                  "Velmeshev_3rd_trimester",        
                  "Velmeshev_0_1_years",            
                  "Velmeshev_1_2_years",            
                  "Velmeshev_2_4_years",           
                  "Velmeshev_10_20_years",        
                  "Velmeshev_Adults",   
                  "GSE157827_Healthy",              
                  "GSE174367_Healthy",              
                  "PRJNA544731_Healthy",    
                  "GSE157827_Alzheimer's disease",  
                  "GSE174367_Alzheimer's disease",  
                  "PRJNA544731_Multiple Sclerosis" 
)

all_degs <- ImportDatasets("data/DEGs/")

all_degs_common <- UnifyAnnotation(all_degs, unified_annotation)

############################################################################################

# 3. UI

ui <- dashboardPage(
  skin = "purple",
  # title of web app
  dashboardHeader(title = "HumanBrainSexSingleCell", titleWidth = 300),
  
  # tabs
  dashboardSidebar(
    width = 300,
    sidebarMenu(
    menuItem("Summary", tabName = "summary", icon = icon("dashboard")),
    menuItem("Datasets information", tabName = "ds_info", icon = icon("folder-open", lib="glyphicon")),
    menuItem("Differential Expressed Genes Analysis", tabName = "degs", icon = icon("stats", lib="glyphicon")),
    menuItem("Open source code", tabName = "source_code", icon = icon("cog", lib="glyphicon"))
  )),
  
  dashboardBody(
      tabItems(
        tabItem(tabName = "summary",
                h2("HumanBrainSexSingleCell"),
                AboutUI("Summary")
        ),
        
        tabItem(tabName = "ds_info",
                h2("Datasets information"),
                DsUI("Datasets information")
        ),
        
        tabItem(tabName = "degs",
                h2("Differential Expressed Genes Analysis"),
                DEGsUI("DEGs"),
                column(6, textOutput("selected_pval")),
                column(6, textOutput("selected_FC"))
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

server <- function(input, output) { 
  
  output$selected_pval <- renderText({ 
    paste0("You have selected this p-value threshold: ", input$pval)
  })
  
  output$selected_FC <- renderText({ 
    paste0("You have selected this FC threshold: ", input$FC)
  })
  
  }


############################################################################################

# 5. Runs app

shinyApp(ui, server)
