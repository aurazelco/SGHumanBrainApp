library(shiny)
library(shinydashboard)

source("custom_scripts/DEGs_func.R")

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

brewer_palette <- c(colorRampPalette(c("white", "#228B22"))(8), "#FF0000")
custom_palette <- c(
  "Velmeshev_2nd_trimester"=brewer_palette[1],           
  "Velmeshev_3rd_trimester"=brewer_palette[2], 
  "Velmeshev_0_1_years"=brewer_palette[3],                
  "Velmeshev_1_2_years"=brewer_palette[4],            
  "Velmeshev_2_4_years"=brewer_palette[5],  
  "Velmeshev_10_20_years"=brewer_palette[6],      
  "Velmeshev_Adults"=brewer_palette[7],
  "GSE157827_Healthy"=brewer_palette[8],              
  "GSE174367_Healthy"=brewer_palette[8],               
  "PRJNA544731_Healthy"=brewer_palette[8], 
  "GSE157827_Alzheimer's disease"=brewer_palette[9],
  "GSE174367_Alzheimer's disease"=brewer_palette[9],
  "PRJNA544731_Multiple Sclerosis"=brewer_palette[9] 
)


DEGsUI <- function(id, label = "degs"){
  fluidPage(
    # First Section
    fluidRow(
      br(),
      p("Please select the threshold for the p-value and fold change (FC). The FC threshold will be then transformed to log10 values."),
      p("Only the DEGs which have both the p-value and the FC less than the respective thresholds will be considered in the analyses below."),
      br(),
      p("Please wait few moments for all graphs to be generated."),
      br(),
      column(6, 
             numericInput(NS(id,"pval"), 
                          p("P-value threshold"), 
                          value = 0.05)),
      column(6, 
             numericInput(NS(id,"FC"), 
                          p("FC threshold"), 
                          value = 1.2))
    ),
    fluidRow(column(6, textOutput(NS(id,"selected_pval"))),
             column(6, textOutput(NS(id,"selected_FC")))
    ),
    br(),
    br(),
    box( height = 1300, width = 900,
        title='The number of sex-biased genes',
        plotOutput(NS(id,"num_degs_plot"))
    ),
    br(),
    fluidRow(column(3, downloadButton(NS(id,'save_num_degs_plot'), 'Download plot as PNG'))
    )
  )
}

DEGsServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
  
    output$selected_pval <- renderText({ 
      paste0("You have selected this p-value threshold: ", input$pval)
    })
    
    output$selected_FC <- renderText({ 
      paste0("You have selected this FC threshold: ", input$FC)
    })
    
    filt_DEGs <- reactive({
      FilterDs(all_degs, input$pval, input$FC)
    })
    
    presence_df <- reactive({
      CreateSexDf(filt_DEGs(), unified_annotation)
    })
    
    num_degs <- reactive({
      NumDEGsAcrossGroups(presence_df(), groups_order)
    })
    
    num_degs_plot <- reactive({
      PlotNumDEGsFaceted(num_degs_final(), custom_palette)
    })
    
    num_degs_plot <- reactive({
      if (input$pval!= 0.05 | input$pval!= 1.2 )  {
        filt_DEGs_filt <- FilterDs(all_degs, input$pval, input$FC) 
        presence_df_filt <- CreateSexDf(filt_DEGs_filt, unified_annotation)
        num_degs_filt <- NumDEGsAcrossGroups(presence_df_filt, groups_order)
        num_degs_plot_filt <- PlotNumDEGsFaceted(num_degs_filt, custom_palette)
      } 
      return(num_degs_plot_filt)
    })
    
    output$num_degs_plot <- renderPlot({
      num_degs_plot()
    }, height = 1200, width = 1000)
    
    output$save_num_degs_plot <- downloadHandler(
      filename = "Number_of_DEGs.png",
      content = function(file) {
        png(file, res = 300, units = 'in', height = 22, width = 15)
        plot(num_degs_plot())
        dev.off()
      }
    )
    
  })
  
}