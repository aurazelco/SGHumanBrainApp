library(shiny)
library(shinydashboard)




DEGsUI <- function(id, label = "degs"){
  fluidPage(
    # First Section
    fluidRow(
      br(),
      p("Please select the threshold for the p-value and fold change (FC). The FC threshold will be then transformed to log10 values."),
      p("Only the DEGs which have both the p-value and the FC less than the respective thresholds will be considered in the analyses below."),
      br(),
      column(6, 
             sliderInput("pval", p("P-value threshold"),
                         min = 0, max = 1, value = 0.05)),
      column(6, 
           sliderInput("FC", p("FC threshold"),
                       min = 1, max = 2, value = 1.2))
             
    )
  )
}