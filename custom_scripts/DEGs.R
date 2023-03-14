library(shiny)
library(shinydashboard)

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

#all_degs <- ImportDatasets("data/Unfiltered_DEGs/")

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
      p("Please select the threshold for the p-value and fold change (FC). The FC threshold will be then transformed to log2 values."),
      p("Only the DEGs which have both the p-value and the FC less than the respective thresholds will be considered in the analyses below."),
      br(),
      column(4, 
             selectInput(NS(id,"pval"), 
                          p("Adjusted P-value threshold"), 
                          choices = list("NS"=1, "0.05"=0.05, "0.01"=0.01, "0.001"=0.001), 
                          selected = 0.05)),
      column(4, 
             selectInput(NS(id,"FC"), 
                          p("FC threshold"), 
                         choices = list("NS"=1, "1.2"=1.2, "1.5"=1.5, "2"=2), 
                         selected = 1.2)),
      column(4, 
             p("Thresholds Help"),
             helpText("Note: in the options on the left, NS stands for non-significant;",
                      "therefore, the DEGs will not be filtered for either adjusted p-value",
                      "or fold change"))
    ),
    fluidRow(column(4, textOutput(NS(id,"selected_pval"))),
             column(4, textOutput(NS(id,"selected_FC")))
    ),
    br(),
    p("download CSV files of filtered DEGs?", style = "color:red"),
    br(),
    p("Please wait few moments for all graphs to be generated."),
    br(),
    br(),
    box( height = 1300, width = 900,
        title='The number of sex-biased DEGs',
        imageOutput(NS(id,"num_degs_plot"))
    ),
    br(),
    fluidRow(column(3, downloadButton(NS(id,'save_num_degs_plot'), 'Download plot as PNG'))
    ),
    br(),
    box( height = 1300, width = 900,
         title='Shared DEGs among datasets, colored by chromosome annotation',
         imageOutput(NS(id,"chr_fraction"))
    ),
    br(),
    fluidRow(column(3, downloadButton(NS(id,'save_chr_fraction_plot'), 'Download plot as PNG'))
    )
  )
}

DEGsServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
  
    output$selected_pval <- renderText({ 
      paste0("You have selected this adjusted p-value threshold: ", input$pval)
    })
    
    output$selected_FC <- renderText({ 
      paste0("You have selected this FC threshold: ", input$FC)
    })
    
    
    filt_DEGs <- reactive({
      filt_files_path <-paste0("data/Filtered_DEGs/pval_", 
                                      str_replace(input$pval, "\\.",  ","), 
                                      "_FC_", 
                                      str_replace(input$FC, "\\.",  ","), "/")
      filt_DEGs <- list.files(filt_files_path, "*.csv", recursive = T)
      return(filt_DEGs)
    })
    
    
    output$num_degs_plot <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                "/Number_of_DEGs.png")))
      list(src = filename, height = 1200, width = 750)
    }, deleteFile = FALSE)
    
    
    output$save_num_degs_plot <- downloadHandler(
      filename = paste0("Number_of_DEGs_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Number_of_DEGs.png"))), file)
      }
    )
    
    output$chr_fraction <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Chr_fractions.png")))
      list(src = filename, height = 1200, width = 750)
    }, deleteFile = FALSE)
    
    output$save_chr_fraction_plot <- downloadHandler(
      filename = paste0("Chr_fraction_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Chr_fractions.png"))), file)
      }
    )  
    
  })
  
}