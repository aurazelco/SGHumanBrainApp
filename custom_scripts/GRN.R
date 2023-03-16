library(shiny)
library(shinydashboard)

linebreaks <- function(n){HTML(strrep(br(), n))}

GRNUI <- function(id, label = "GRN"){
  fluidPage(
    
    fluidRow(
      br(),
      p("The GRN analysis was performed using SCENIC (ref). A minimum of 100 cells per sex and cell type was used to filter out cell types with too few cells. "),
      p("The cell types with 100 cels or more were randomly sampled 3 times, 100 cells selected every time for each cell type. "),
      p("An expression matrix of the 2000 most variable genes was extracted for each set of randomly sampled cells. "),
      p("Thus, we ran through SCENIC 3 randomly sampled expression matrices, 100 celss per cell type after filtering, for each sex separately. ")
      
    ),

  )
}

GRNServer <- function(id) {
  
  moduleServer(id, function(input, output, session) {
  
  
  })
  
}