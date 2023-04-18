
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
      p(strong("Background:"),  " Sex and gender have only lately been systematically regarded as a biological variable in both pre-clinical and clinical investigations. The impact of sex and gender on a wide range of biological and psychological variables remains largely unknown during brain development and ageing disorders."),
      p(strong("Methods:"),  " To systematically evaluate sex and gender (SG) differences at different development stages and ageing disorders at a single cell level, we gathered publicly available single-nucleus RNA-sequencing studies through human life span from second trimester of gestation until geriatric age in healthy individuals and from Alzheimer's disease (AD) and Multiple Sclerosis (MS) patients. In summary, we collected single cell data for a total of 419885 single nuclei from 161 human brain samples (72 females and 89 males) to identify and characterise SG-biased genes."),
      p(strong("Results:"),  " We identified SG-biased genes in both females and males in 11 major brain cell types across 7 developmental stages and two brain disorders. SG-biased genes were located mostly on the autosomes, with an enrichment for Y-linked genes in males in some cell types. SG-biased genes in many cell types were enriched for cell type markers. Accordingly, SG-biased genes showed little overlap across cell types and developmental stages. Interestingly, there was extensive functional overlap across SG-biased genes in developmental stages. Female-biased genes were enriched for brain-related functions and processes, and male-biased genes were enriched for metabolic pathways. Common female-biased genes across cell types and developmental stages contained many mitochondrial genes. Investigation of hormonal influence identified thymosin targets enriched in male-biased genes, and SG-biased genes for both males and females were highly enriched for androgen (not oestrogen) response elements across cell types"),
      p(strong("Conclusion:"),  " We systematically characterised SG differences in brain development and brain-related disorders at a single cell level, leading to the identification of hormonal influences likely establishing these differences as well as enriched pathways and functional categories likely contributing to the SG differences in brain-related disorders. We have further made the entire analysis available as a web resource for the scientific community by developing a web application which can be found here."),
      p(strong("Keywords:"), " sex; gender; sex and gender differences; human brain; transcriptome; single-nucleus RNA-seq")
    )
  )
}