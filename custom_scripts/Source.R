library(shiny)
library(shinydashboard)

SourceUI <- function(id, label = "source_code"){
  fluidPage(
    # First Section
    fluidRow(
      br(),
      p("The open source code for this application is available on Github :"),
      tags$a(href = "https://github.com/aurazelco/SGHumanBrainApp.git", "https://github.com/aurazelco/SGHumanBrainApp.git")
      
    ),
    
    fluidRow(
      br(),
      p("The bulk RNA-seq application was also developed indipendently, please check the original source at: ", ),
      tags$a(href = "https://github.com/PattaWapee/SexRankBrain", "https://github.com/PattaWapee/SexRankBrain")
      
    )
    
  )
}

