library(shiny)
library(shinydashboard)
library(stringr)
library(shinyjs)

linebreaks <- function(n){HTML(strrep(br(), n))}

DEGsUI <- function(id, label = "degs"){
  fluidPage(
    shinyjs::useShinyjs(),
    
    fluidRow(
      br(),
      p("Please select the threshold for the p-value and fold change (FC). The FC threshold will be then transformed to log2 values."),
      p("Only the DEGs which have both the p-value and the FC less than the respective thresholds will be considered in the analyses below."),
      br(),
      h4("Parameter thresholds:"),
      column(4, 
             selectInput(NS(id,"pval"), 
                          p("Adjusted P-value threshold"), 
                          choices = list("NS"=1, "0.05"=0.05, "0.01"=0.01, "0.001"=0.001), 
                          selected = 0.05),
             br(),
             textOutput(NS(id,"selected_pval"))),
      
      column(4, 
             selectInput(NS(id,"FC"), 
                          p("FC threshold"), 
                         choices = list("NS"=1, "1.2"=1.2, "1.5"=1.5, "2"=2), 
                         selected = 1.2),
             br(),
             textOutput(NS(id,"selected_FC"))),
      
      column(4, 
             p("Thresholds Help"),
             helpText("Note: in the options on the left, NS stands for non-significant;",
                      "therefore, the DEGs will not be filtered for either adjusted p-value",
                      "or fold change"))
    ),
    linebreaks(4),
    h4("Download files as ZIP folder?"),
    fluidRow(
      column(3,
             checkboxGroupInput(NS(id, "zips"), "Files to be downloaded as ZIP:",
                                choices = c("Plots"="plots",
                                            "Unfiltered DEGs" = "unfiltered DEGs",
                                            "Filtered DEGs" = "filtered DEGs"))),
      
      column(5, downloadButton(NS(id,'save_zip'), 'Download selected files as a ZIP folder')),
      column(6, 
             helpText("Please note that the few moments may be required to download the files, ",
                      "especially if all three options are chosen"))
    ),
    linebreaks(5),
    h4("Output plots"),
    linebreaks(2),
    fluidRow(
      tabBox(width = 18,
             height = 1400,
             title = "Overview pof sex-biased DEGs",
             id = "Overview_DEGs",
             tabPanel(strong("The number of sex-biased DEGs"),
                      p("Below the number of sex-biased DEGs in each group for each sex and cell type can be found. "),
                      linebreaks(2),
                      imageOutput(NS(id,"num_degs_plot")),
                      linebreaks(45),
                      downloadButton(NS(id,'save_num_degs_plot'), 'Download plot as PNG'),
             ),
             
             tabPanel(strong("Shared DEGs among datasets, colored by chromosome annotation"),
                      p("Below the fractions of which chromosome (Autosome, X or Y) the sex-biased DEGs belong to, ordered by how many groups share the DEGs. The fractions are in logaritmic scale to enhance the lower fractions. "),
                      linebreaks(2),
                      imageOutput(NS(id,"chr_fraction")),
                      linebreaks(45),
                      downloadButton(NS(id,'save_chr_fraction_plot'), 'Download plot as PNG')
             ),
             
             tabPanel(strong('Sex-biased gene locations'),
                      p("In the plot, the frequency of cellular location is indicated by the size, while the color represents the cellular compartment. "),
                      linebreaks(2),
                      imageOutput(NS(id,"GeneLocation")),
                      linebreaks(45),
                      downloadButton(NS(id,'save_GeneLocation_plot'), 'Download plot as PNG')
             )
      )
    ),
    linebreaks(5),
    fluidRow(
      tabBox(width = 18,
             height = 1000,
             title = "Presence heatmaps of specific genes of interest",
             id = "Ind_Hmps",
             tabPanel(strong("Top 20 most differentially present genes"),
                      imageOutput(NS(id,"mostdiffgenes")),
                      linebreaks(20),
                      downloadButton(NS(id,'save_mostdiffgenes_plot'), 'Download plot as PNG')
              ),
             
             tabPanel(strong("Mitochondrial genes"),
                      imageOutput(NS(id,"MTgenes")),
                      linebreaks(20),
                      downloadButton(NS(id,'save_MTgenes_plot'), 'Download plot as PNG')
              ),
             
             tabPanel(strong('X-escaping genes'),
                      imageOutput(NS(id,"Xescapinggenes")),
                      linebreaks(20),
                      downloadButton(NS(id,'save_Xescapinggenes_plot'), 'Download plot as PNG')
              )
      )
    ),
    linebreaks(2),
    box( height = 900, width = 1200,
         strong('Percentage of cell type markers from McKenzie et al. 2018'),
         imageOutput(NS(id,"McKenzie"))
    ),
    column(3,
           downloadButton(NS(id,'save_McKenzie_plot'), 'Download plot as PNG')
    ),
    linebreaks(2),
    box( height = 900, width = 1200,
         strong("Presence heatmaps of disease markers from Chlamydas et al. 2022"),
         imageOutput(NS(id,"Chlamydas"))
    ),
    column(3,
           downloadButton(NS(id,'save_Chlamydas_plot'), 'Download plot as PNG')
    ),
    linebreaks(2),
    box( height = 1300, width = 900,
         strong("Hormone targets enrichment from Jadhav et al. 2022"),
         imageOutput(NS(id,"Hormones"))
    ),
    column(3,
           downloadButton(NS(id,'save_Hormones_plot'), 'Download plot as PNG')
    ),
    linebreaks(2),
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
    linebreaks(5),
    fluidRow(
      tabBox(width = 18,
             height = 2000,
             title = "Functional Enrichment",
             id = "Functional_enrichment",
             tabPanel(strong("Gene Ontology (GO) - Biological Processes"),
                      imageOutput(NS(id,"GOBP")),
                      linebreaks(70),
                      downloadButton(NS(id,'save_GOBP_plot'), 'Download plot as PNG')
             ),
             
             tabPanel(strong("Gene Ontology (GO) - Cellular Components"),
                      imageOutput(NS(id,"GOCC")),
                      linebreaks(70),
                      downloadButton(NS(id,'save_GOCC_plot'), 'Download plot as PNG')
             ),
             
             tabPanel(strong("Gene Ontology (GO) - Molecular Functions"),
                      imageOutput(NS(id,"GOMF")),
                      linebreaks(70),
                      downloadButton(NS(id,'save_GOMF_plot'), 'Download plot as PNG')
             ),
             
             tabPanel(strong("Kyoto Encyclopedia of Genes and Genomes (KEGG) - pathway enrichment"),
                      imageOutput(NS(id,"KEGG")),
                      linebreaks(70),
                      downloadButton(NS(id,'save_KEGG_plot'), 'Download plot as PNG')
             )
      )
    ),
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
    
    observe({
      shinyjs::toggleState("save_zip", condition = !is.null(input$zips))
    })
 
    output$save_zip <- downloadHandler(
      filename = paste0("HumanBrainSexSingleCell_download.zip"),
      contentType = "application/zip",
      content = function(file){
        
        plt_dir <- paste0("www/Plots/pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","))
        plt_files <- list.dirs(plt_dir, full.names = T, recursive = T)
        
        unfilt_DEGs_dir <- "data/Unfiltered_DEGs"
        unfilt_files <- list.dirs(unfilt_DEGs_dir, full.names = T, recursive = T)
        
        filt_DEGs_dir <- paste0("data/Filtered_DEGs/pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","))
        filt_files <- list.dirs(filt_DEGs_dir, full.names = T, recursive = T)
        
        files_ls <- list("plots"=plt_files,
                         "unfiltered DEGs"=unfilt_files,
                         "filtered DEGs"=filt_files)
        
        files <- unlist(files_ls[names(files_ls) %in% input$zips], use.names = F)

        #create the zip file
        zip(file, files)
        
      }
    )
    
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
    
    output$GeneLocation <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Location.png")))
      list(src = filename, height = 1200, width = 750)
    }, deleteFile = FALSE)
    
    output$save_GeneLocation_plot <- downloadHandler(
      filename = paste0("Location_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/Location.png"))), file)
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
    
    output$save_MTgenes_plot <- downloadHandler(
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
    
    output$GOBP <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/GO_BP.png")))
      list(src = filename, height = 1600, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_GOBP_plot <- downloadHandler(
      filename = paste0("GO_BP_sites_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/GO_BP.png"))), file)
      }
    )  
    
    output$GOCC <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/GO_CC.png")))
      list(src = filename, height = 1600, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_GOCC_plot <- downloadHandler(
      filename = paste0("GO_CC_sites_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/GO_CC.png"))), file)
      }
    )  
    
    output$GOMF <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/GO_MF.png")))
      list(src = filename, height = 1600, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_GOMF_plot <- downloadHandler(
      filename = paste0("GO_MF_sites_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/GO_MF.png"))), file)
      }
    )  
    
    output$KEGG <- renderImage({
      filename <- normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/KEGG.png")))
      list(src = filename, height = 1600, width = 1400)
    }, deleteFile = FALSE)
    
    output$save_KEGG_plot <- downloadHandler(
      filename = paste0("KEGG_sites_pval_", str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), ".png"),
      contentType = "image/png",
      content = function(file) {
        file.copy(normalizePath(file.path(paste0("www/Plots/pval_", 
                                                 str_replace(input$pval, "\\.",  ","), "_FC_", str_replace(input$FC, "\\.",  ","), 
                                                 "/KEGG.png"))), file)
      }
    )  
    
  
  })
  
}