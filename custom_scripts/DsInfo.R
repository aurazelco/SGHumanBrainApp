library(shiny)
library(shinydashboard)

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
      p("For each sub-dataset, as described below, we compared the genes expressed in females and males for each cell type, therefore obtaining sex-biased gene lists in each cell type and sub-dataset."),
      br(),
      h4("Velmeshev et al. 2022"),
      p("This single-nucleus RNA-seq dataset collects healthy human brain samples, from the second trimester of gestation until adulthood (defined as more than 20 years of age)."),
      p("In our analysis, we analyzed the ages separately. "),
      br(),
      h4("Li et al. 2022 - DISCO v1.0"),
      p("This dataset collects human brain samples from previously published data. Here only three projects are included (GSE157827, GSE174367 and PRJNA544731), since they had complete sex metadata and samples from both females and males. "),
      p("The included samples were from adult individuals (35-94 years of age), and the individuals were either healthy or patients suffering from Alzheimer's disease or multiple sclerosis."),
      br(),
      br(),
      h4("Samples information"),
      p("Below the number of samples and cells in each dataset is plotted, as well as the numebr of cells per cell type. "),
      column(6,
             h5("Number of samples in the sub-datasets, divided by sex"),
             img(src = "num_samples.png", height="400px", width="100%")),
      column(6,
             h5("Number of cells in the sub-datasets, divided by sex"),
             img(src = "num_cells.png", height="400px", width="70%"))
    ),
    fluidRow(
      br(),
      br(),
      h5("Number of cells in the sub-datasets, divided by sex and cell type"),
      img(src = "num_cells_per_ct.png", height="50%", width="50%", style="display: block; margin-left: auto; margin-right: auto;"),
      p("The dashed line corresponds to 100 cells, which was used as a threshold for determinng which cell types to analyze. ")
    )
  )
}