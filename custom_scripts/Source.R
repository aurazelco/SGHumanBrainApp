library(shiny)
library(shinydashboard)

SourceUI <- function(id, label = "source_code"){
  fluidPage(
    # First Section
    fluidRow(
      br(),
      p("The open source code for this application is available on Github :"),
      tags$a(href = "https://github.com/aurazelco/HumanBrainSexSingleCell", "https://github.com/aurazelco/HumanBrainSexSingleCell")
      
    )
    
  )
}

