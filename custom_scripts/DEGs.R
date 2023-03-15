library(shiny)
library(shinydashboard)
#library(zip)

DEGsUI <- function(id, label = "degs"){
  fluidPage(
    
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
    downloadButton(NS(id,'save_all_plots'), 'Download all plots as a ZIP folder'),
    br(),
    br(),
    fluidRow(
         column(6,
                box( height = 1300, width = 450,
                strong('The number of sex-biased DEGs'),
                imageOutput(NS(id,"num_degs_plot"))),
                downloadButton(NS(id,'save_num_degs_plot'), 'Download plot as PNG')),
         column(6,
                box( height = 1300, width = 450,
                strong('Shared DEGs among datasets, colored by chromosome annotation'),
                imageOutput(NS(id,"chr_fraction"))),
                downloadButton(NS(id,'save_chr_fraction_plot'), 'Download plot as PNG'))
    ),
    br(),
    br(),
    box( height = 900, width = 1200,
        strong('Top 20 most differentially present genes'),
        imageOutput(NS(id,"mostdiffgenes"))),
    column(3,
           downloadButton(NS(id,'save_mostdiffgenes_plot'), 'Download plot as PNG')
           ),
    br(),
    br(),
    box( height = 900, width = 1200,
        strong('Mitochondrial genes'),
        imageOutput(NS(id,"MTgenes"))
        ),
    column(3,
           downloadButton(NS(id,'save_MTgenes_plot'), 'Download plot as PNG')
    ),
    br(),
    br(),
    box( height = 900, width = 1200,
        strong('X-escaping genes'),
        imageOutput(NS(id,"Xescapinggenes"))
    ),
    column(3,
           downloadButton(NS(id,'save_Xescapinggenes_plot'), 'Download plot as PNG')
    ),
    br(),
    br(),
    box( height = 900, width = 1200,
         strong('Percentage of cell type markers from McKenzie et al. 2018'),
         imageOutput(NS(id,"McKenzie"))
    ),
    column(3,
           downloadButton(NS(id,'save_McKenzie_plot'), 'Download plot as PNG')
    ),
    br(),
    br(),
    box( height = 900, width = 1200,
         strong("Presence heatmaps of disease markers from Chlamydas et al. 2022"),
         imageOutput(NS(id,"Chlamydas"))
    ),
    column(3,
           downloadButton(NS(id,'save_Chlamydas_plot'), 'Download plot as PNG')
    ),
    br(),
    br(),
    box( height = 1300, width = 900,
         strong("Hormone targets enrichment from Jadhav et al. 2022"),
         imageOutput(NS(id,"Hormones"))
    ),
    column(3,
           downloadButton(NS(id,'save_Hormones_plot'), 'Download plot as PNG')
    ),
    br(),
    br(),
    fluidRow(
      column(6,
             box( height = 900, width = 500,
                  strong('Percentages of ARE sites among the sex-biased DEGs'),
                  imageOutput(NS(id,"ARE"))),
             downloadButton(NS(id,'save_ARE_plot'), 'Download plot as PNG')),
      column(6,
             box( height = 900, width = 500,
                  strong('Percentages of ERE sites among the sex-biased DEGs'),
                  imageOutput(NS(id,"ERE"))),
             downloadButton(NS(id,'save_ERE_plot'), 'Download plot as PNG'))
    ),
    br(),
    br(),
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
    
    output$mostdiffgenes <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/top_20_most_diff_genes.png")))
      list(src = filename, height = 800, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_mostdiffgenes_plot <- downloadHandler(
      filename = paste0("Top_20_most_diff_genes_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/top_20_most_diff_genes.png"))), file)
      }
    )  
    
    output$MTgenes <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/MT_genes.png")))
      list(src = filename, height = 800, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_mostdiffgenes_plot <- downloadHandler(
      filename = paste0("MT_genes_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/MT_genes.png"))), file)
      }
    )  
    
    output$Xescapinggenes <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/X_escaping_genes.png")))
      list(src = filename, height = 800, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_Xescapinggenes_plot <- downloadHandler(
      filename = paste0("X_escaping_genes_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/X_escaping_genes.png"))), file)
      }
    )  
    
    
    output$McKenzie <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/McKenzie_perc.png")))
      list(src = filename, height = 800, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_McKenzie_plot <- downloadHandler(
      filename = paste0("McKenzie_percentages_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/McKenzie_perc.png"))), file)
      }
    )  
    
    output$Chlamydas <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Chlamydas_hmp.png")))
      list(src = filename, height = 800, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_Chlamydas_plot <- downloadHandler(
      filename = paste0("Chlamydas_hmp_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Chlamydas_hmp.png"))), file)
      }
    )  
    
    output$Hormones <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Hormone_target_enrichment.png")))
      list(src = filename, height = 1200, width = 800)
    }, deleteFile = FALSE)
    
    output$save_Hormones_plot <- downloadHandler(
      filename = paste0("Hormone_target_enrichment_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Hormone_target_enrichment.png"))), file)
      }
    )  
    
    output$ARE <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/ARE.png")))
      list(src = filename, height = 800, width = 700)
    }, deleteFile = FALSE)
    
    output$save_ARE_plot <- downloadHandler(
      filename = paste0("ARE_sites_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/ARE.png"))), file)
      }
    )  
    
    output$ERE <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/ERE.png")))
      list(src = filename, height = 800, width = 700)
    }, deleteFile = FALSE)
    
    output$save_ERE_plot <- downloadHandler(
      filename = paste0("ERE_sites_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/ERE.png"))), file)
      }
    )  

    
    output$save_all_plots <- downloadHandler(
      filename = paste0("All_plots_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".zip"),
        contentType = "application/zip",
      content = function(file){
        plot_ls <- c(
          "Number_of_DEGs.png",
          "Chr_fractions.png",
          "top_20_most_diff_genes.png",
          "MT_genes.png",
          "X_escaping_genes.png",
          "McKenzie_perc.png",
          "Chlamydas_hmp.png",
          "Hormone_target_enrichment.png",
          "ARE.png",
          "ERE.png"
        )
        files <- NULL
        for (i in plot_ls){
          files <- c(paste0("www/Plots/pval_", 
                    str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                    "/", i), files)
        }
        #create the zip file
        zip(file, files)
      }
    )
    
    
  
  
  })
  
}