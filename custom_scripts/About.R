
library(shiny)
library(shinydashboard)

AboutUI <- function(id, label = "summary"){
  fluidPage(
    # First Section
    fluidRow(
      br(),
      p(strong("SGHumanBrainApp"), " is a web application to explore the results from sex-biased differentially expressed genes (DEGs) in the human brain."),
      p("The web application contains DEGs analysis from both single-cell and bulk RNA-seq, performed in two separate papers. The tabs indicate from which analysis the data come from."),
      p("Please refer to the abstracts below and the papers for more information."),
      br(),
      h3("Single-cell RNA-seq paper",  style = "color:purple"),
      h3("Exploring the sex and gender differences in the human brain at single cell level"), 
      h4("Abstract"),
      p(strong("Background:"),  " Sex and gender have only lately been systematically regarded as a biological variable in both pre-clinical and clinical investigations. The impact of sex and gender on a wide range of biological and psychological variables remains largely unknown during brain development and ageing disorders."),
      p(strong("Methods:"),  " To systematically evaluate sex and gender (SG) differences at different development stages and ageing disorders at a single cell level, we gathered publicly available single-nucleus RNA-sequencing studies through human life span from second trimester of gestation until geriatric age in healthy individuals and from Alzheimer's disease (AD) and Multiple Sclerosis (MS) patients. In summary, we collected single cell data for a total of 419885 single nuclei from 161 human brain samples (72 females and 89 males) to identify and characterise SG-biased genes."),
      p(strong("Results:"),  " We identified SG-biased genes in both females and males in 11 major brain cell types across 7 developmental stages and two brain disorders. SG-biased genes were located mostly on the autosomes, with an enrichment for Y-linked genes in males in some cell types. SG-biased genes in many cell types were enriched for cell type markers. Accordingly, SG-biased genes showed little overlap across cell types and developmental stages. Interestingly, there was extensive functional overlap across SG-biased genes in developmental stages. Female-biased genes were enriched for brain-related functions and processes, and male-biased genes were enriched for metabolic pathways. Common female-biased genes across cell types and developmental stages contained many mitochondrial genes. Investigation of hormonal influence identified thymosin targets enriched in male-biased genes, and SG-biased genes for both males and females were highly enriched for androgen (not oestrogen) response elements across cell types"),
      p(strong("Conclusion:"),  " We systematically characterised SG differences in brain development and brain-related disorders at a single cell level, leading to the identification of hormonal influences likely establishing these differences as well as enriched pathways and functional categories likely contributing to the SG differences in brain-related disorders. We have further made the entire analysis available as a web resource for the scientific community by developing a web application which can be found here."),
      p(strong("Keywords:"), " sex; gender; sex and gender differences; human brain; transcriptome; single-nucleus RNA-seq"),
      br(),
      h3("Bulk RNA-seq paper",  style = "color:purple"),
      h3("Integrated analysis of robust sex-biased gene signatures in human brain"), 
      h4("Abstract"),
      p(strong("Background:"),  " Sex dimorphism is highly prominent in mammals with
                         many systematic differences between male and female 
                         form of the species. Sex differences arise from genetic
                         and environmental factors. Accordingly, the fundamental 
                         social and cultural stratification factors for humans is 
                         sex. It distinguishes individuals most prominently on the 
                         reproductive traits, but also affects many of the other 
                         related traits. This can lead to different disease 
                         susceptibilities and treatment responses across sexes. 
                         Sex differences in brain have raised a lot of controversy 
                         due to small and sometimes contradictory sex-specific 
                         effects. Many studies have been published to identify 
                         sex-biased genes in one or several brain regions but 
                         the assessment of the robustness of these studies is 
                         missing. We therefore collected huge amount of publicly
                         available transcriptomic data to first estimate whether
                         consistent sex differences exists and further explore 
                         their likely origin and functional significance."),
      p(strong("Results:"),  " In order to systematically characterise sex specific 
                         differences across human brain regions, we collected 
                         transcription profiles for more than 16000 samples from 
                         over 46 datasets across eleven brain regions. By systematic 
                         integration of the data from multiple studies, we identified 
                         robust transcription level differences in human brain across eleven 
                         brain regions and classified male-biased and female-biased genes. 
                         Firstly, both male and female-biased genes were highly conserved across
                         primates and showed a high overlap with sex-biased genes in other species. 
                         Female-biased genes were enriched for neuron-associated processes while 
                         male-biased genes were enriched for membranes and nuclear structures. 
                         Male-biased genes were enriched on the Y chromosome while female-biased 
                         genes were enriched on the X chromosome, which included X chromosome inactivation 
                         escapees explaining the origins of some sex differences. We noted that age is 
                         a major co-variate of sex-differences. Finally, many more female-biased genes 
                         were affected by adverse drug reactions than male-biased genes."),
      p(strong("Conclusion:"), " In summary, by building a comprehensive resource of sex differences across 
                         human brain regions at gene expression level, we explored their likely origin and 
                         functional significance. We have also developed a web resource to make the entire
                         analysis available for the scientific community for further exploration."),
      p(strong("Keywords:"), "  sex difference; human brain; gene regulation; hormones; data integration; conservation; brain disorders; drug response")
    )
  )
}