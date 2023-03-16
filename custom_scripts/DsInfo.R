library(shiny)
library(shinydashboard)

linebreaks <- function(n){HTML(strrep(br(), n))}

DsUI <- function(id, label = "ds_info"){
  fluidPage(
    # First Section
    fluidRow(
      br(),
      p("This analysis includes two major datasets, ",
        tags$a(href = "https://www.biorxiv.org/content/10.1101/2022.10.24.513555v1", "Velmeshev et al. 2022"),
        " and ",
        tags$a(href = "https://www.immunesinglecell.org", "Li et al. 2022 (DISCO v1.0)"),
        "."),
      p("For each sub-dataset, as described below, we compared the genes expressed in females and males for each cell type, therefore obtaining sex-biased differentially expressed gene lists (DEGs) in each cell type and sub-dataset."),
      p("Furthermore, we also perfomed a Gene Regulatory Network analysis (GRN), to investigate if any transcription factor was predicted to be sex-biased. "),
      br(),
      h4("Velmeshev et al. 2022"),
      p("This single-nucleus RNA-seq dataset collects healthy human brain samples, from the second trimester of gestation until adulthood (defined as more than 20 years of age)."),
      p("In our analysis, we analyzed the ages separately. "),
      br(),
      h4("Li et al. 2022 - DISCO v1.0"),
      p("This dataset collects human brain samples from previously published data. Here only three projects are included (GSE157827, GSE174367 and PRJNA544731), since they had complete sex metadata and samples from both females and males. "),
      p("The included samples were from adult individuals (35-94 years of age), and the individuals were either healthy or patients suffering from Alzheimer's disease or multiple sclerosis."),
      linebreaks(2), 
      h4("Samples information"),
      p("Below the number of samples and cells in each dataset is plotted, as well as the number of cells per cell type. "),

      ),
    br(),
    fluidRow(
      tabBox(width = 18,
             height = 1000,
             title = "Dataset sample information",
             id = "ds_meta",
             tabPanel(strong("Number of samples in the sub-datasets, divided by sex"),
                      linebreaks(2),
                      p("All datasets but one (Velmeshev 4-10 years) have more than 3 samples in both female and male groups. "),
                      p("The dashed line represent n=3, which we used a sminimum number of samples per each group and sex. "),
                      br(),
                      img(src = "ds_info/num_samples.png", height="55%", width="55%", style="display: block; margin-left: auto; margin-right: auto;")
             ),
             
             tabPanel(strong("Number of cells in the sub-datasets, divided by sex"),
                      linebreaks(2),
                      p("Here the number of cells in each dataset is howed. The datasets are divided by helathy/disease. "),
                      br(),
                      img(src = "ds_info/num_cells.png", height="40%", width="40%", style="display: block; margin-left: auto; margin-right: auto;")    
             ),
             
             tabPanel(strong("Number of cells in the sub-datasets, divided by sex and cell type"),
                      linebreaks(2),
                      p("The dashed line corresponds to 100 cells, which was used as a threshold for determining which cell types to analyze. "),
                      br(),
                      img(src = "ds_info/num_cells_per_ct.png", height="50%", width="50%", style="display: block; margin-left: auto; margin-right: auto;")
             )
      )
    )
  )
}