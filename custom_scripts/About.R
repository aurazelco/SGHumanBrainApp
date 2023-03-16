
library(shiny)
library(shinydashboard)

AboutUI <- function(id, label = "summary"){
  fluidPage(
    # First Section
    fluidRow(
      br(),
      p(strong("HumanBrainSexSingleCell"), " is a web application to explore the results from sex-biased differentially expressed genes (DEGs) in the human brain."),
      p("Please refer to the abstract below and the paper for more information."),
      br(),
      h4("Abstract"),
      p("add text", style = "color:red")
      
    )
    
  )
}