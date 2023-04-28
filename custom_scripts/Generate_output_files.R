
######################### Initial set-up and definition of core variables
wd <- "HumanBrainSexSingleCell/"
source(paste0(wd, "custom_scripts/Generate_output_files_func.R"))

# manually decided how to combine the sub-celltypes
unified_annotation <- c("CXCL14 IN" = "Interneurons",
                        "EC" = "Endothelial cells",
                        "fibrous astrocyte"  = "Astrocytes",
                        "L2_3 EN" = "Excitatory neurons", 
                        "L4 EN" = "Excitatory neurons",
                        "L5 EN" = "Excitatory neurons",
                        "L5_6 EN" = "Excitatory neurons",
                        "L5b EN" = "Excitatory neurons",
                        "L6 EN" = "Excitatory neurons",                
                        "microglia" = "Microglia", 
                        "Oligodendrocyte" =  "Oligodendrocytes",      
                        "OPC" = "OPCs",                  
                        "PLCH1 L4_5 EN" = "Excitatory neurons", 
                        "protoplasmic astrocyte" = "Astrocytes",
                        "PVALB IN"  = "Interneurons",            
                        "pyramidal neuron"  = "Excitatory neurons",
                        "SST IN" = "Interneurons",   
                        "SV2C IN"  = "Interneurons",   
                        "TSHZ2 L4_5 EN" = "Excitatory neurons",  
                        "VIP IN" = "Interneurons",
                        "Mesenchymal" = "Mesenchymal",      
                        "Neuroepithelial" =     "Neuroepithelial",
                        "Neuronal" = "Neurons",            
                        "Other"    = "Other",                
                        "Radial Glial"     = "Radial Glia",       
                        "Astrocytes" = "Astrocytes",        
                        "Excitatory neurons"  = "Excitatory neurons",
                        "Interneurons"   = "Interneurons",     
                        "Microglia"  = "Microglia",         
                        "Oligodendrocytes" = "Oligodendrocytes",
                        "OPCs" = "OPCs",            
                        "Unknown" = "Unknown",           
                        "Vascular cells" = "Vascular cells",     
                        "Dorsal progenitors"  = "Dorsal progenitors" ,   
                        "Ventral progenitors" = "Ventral progenitors")
names(unified_annotation) <- tolower(names(unified_annotation))

# defines the order in which to organize the presence heatmaps, so the groups are in developmental order, with the last groups as diseases
groups_order <- c("Velmeshev_2nd_trimester",        
                  "Velmeshev_3rd_trimester",        
                  "Velmeshev_0_1_years",            
                  "Velmeshev_1_2_years",            
                  "Velmeshev_2_4_years",           
                  "Velmeshev_10_20_years",        
                  "Velmeshev_Adults",   
                  "GSE157827_Healthy",              
                  "GSE174367_Healthy",              
                  "PRJNA544731_Healthy",    
                  "GSE157827_Alzheimer's disease",  
                  "GSE174367_Alzheimer's disease",  
                  "PRJNA544731_Multiple Sclerosis" 
)

cts_order <- c(
  "Excitatory neurons",  
  "Interneurons",        
  "Astrocytes",          
  "Microglia" ,          
  "Oligodendrocytes",    
  "OPCs",   
  "Endothelial cells", 
  "Vascular cells",      
  "Dorsal progenitors",  
  "Ventral progenitors",
  "Unknown"   
)

cts_sex_order <- paste(rep(cts_order, each=length(c("F", "M"))), c("F", "M"), sep = "/")

# Imports original dex-biased DEGs
all_degs <- ImportDatasets(paste0(wd, "data/Unfiltered_DEGs/"))

######################### Generate Filtered DEGs

# Adjusted p-value and FC thresholds
pval_ls <- c(1, 0.05, 0.01, 0.001)
FC_ls <- c(1, 1.2, 1.5, 2)

# Generates filtered files and saves them in specific sub-folders
for (pval_x in pval_ls) {
  for (FC_x in FC_ls) {
    print(paste0("p-value threshold: ", pval_x, "; FC threshold: ", FC_x))
    
    filt_degs <- FilterDs(all_degs, pval_x, FC_x)
    
    out_path <- paste0(wd, "data/Filtered_DEGs/pval_", str_replace(pval_x, "\\.",  ","), "_FC_", str_replace(FC_x, "\\.",  ","), "/")
    dir.create(out_path, recursive = T, showWarnings = F)
    
    for (group_id in names(filt_degs)) {
      for (ct in names(filt_degs[[group_id]])) {
        out_folder <- paste0(out_path, group_id, "/", ct, "/")
        dir.create(out_folder, recursive = T, showWarnings = F)
        lapply(1:length(names(filt_degs[[group_id]][[ct]])), function(x) write.csv(filt_degs[[group_id]][[ct]][[x]],
                                                                                   paste0(out_folder, names(filt_degs[[group_id]][[ct]])[x], "_filt.csv")
        ))
      }
    }
    
  }
}

######################### Extra files needed for plots

# X-escaping genes
x_escapees <- read.table(paste0(wd, "data/extra_files/escape_Xchr.txt"), sep="\t", skip = 2)

# McKenzie et al 2018 - cell type markers
McKenzie <- as.data.frame(read_xlsx(paste0(wd, "data/extra_files/McKenzie_2018_suppl.xlsx"),
                                    sheet = 'top_human_enrich',
                                    skip = 2))

McKenzie_ct_names <- c(
  "ast" = "Astrocytes", 
  "end" = "Endothelial Cells",
  "mic" = "Microglia",
  "neu" = "Neurons",
  "oli" = "Oligodendrocytes"
)

# Chlamydas et al 2022 - table for disease markers
Chlamydas <- as.data.frame(read_xlsx(paste0(wd, "data/extra_files/Chlamydas_2022.xlsx"), skip = 1))
colnames(Chlamydas) <- str_replace_all(colnames(Chlamydas), " ", "_")
Chlamydas <- Chlamydas[, c(1,2,4)]
Chlamydas <- drop_na(Chlamydas)

# Jadhav et al 2022 - Hormone targets 
hormones <- fromJSON(file=paste0(wd, "data/extra_files/hgv1_hormone_genes.json"))
hormones_filt <- hormones[names(which(lapply(hormones, length)>=10))]

# AREs and EREs - Wilson et al. 2016, Bourdeau et al 2004, 
ARE <- read_excel(paste0(wd, "data/extra_files/AREsitesHuman.xlsx"),skip=1)
colnames(ARE) <- c("fullsites", "halfsites")
ARE_sites <- list("full" = unique(setdiff(ARE$fullsites, ARE$halfsites)), 
               "half" = unique(setdiff(ARE$halfsites, ARE$fullsites)),
               "hf" = unique(intersect(ARE$fullsites, ARE$halfsites)))

ERE <- read_excel(paste0(wd, "data/extra_files/Hs_allEREs.xls"))
EREgene <- unique(ERE$`Hs Gene Name`)

# Sets the palette for a plot generated later on
brewer_palette <- c(colorRampPalette(c("white", "#228B22"))(8), "#FF0000")
custom_palette <- c(
  "Velmeshev_2nd_trimester"=brewer_palette[1],           
  "Velmeshev_3rd_trimester"=brewer_palette[2], 
  "Velmeshev_0_1_years"=brewer_palette[3],                
  "Velmeshev_1_2_years"=brewer_palette[4],            
  "Velmeshev_2_4_years"=brewer_palette[5],  
  "Velmeshev_10_20_years"=brewer_palette[6],      
  "Velmeshev_Adults"=brewer_palette[7],
  "GSE157827_Healthy"=brewer_palette[8],              
  "GSE174367_Healthy"=brewer_palette[8],               
  "PRJNA544731_Healthy"=brewer_palette[8], 
  "GSE157827_Alzheimer's disease"=brewer_palette[9],
  "GSE174367_Alzheimer's disease"=brewer_palette[9],
  "PRJNA544731_Multiple Sclerosis"=brewer_palette[9] 
)

location_df <- read_xlsx(paste0(wd, "data/extra_files/Thul_2017_sm_table_s6.xlsx" ))
loc_counts <- c(colSums(location_df[which(location_df$Reliability!="Uncertain"), c(4:32)]),
                "Uncertain"=sum(colSums(location_df[which(location_df$Reliability=="Uncertain"), c(4:32)])))

tot_genes <- 20000


######################### Plots

# Adjusted p-value and FC thresholds
for (pval_x in pval_ls) {
  for (FC_x in FC_ls) {
    print(paste0("p-value threshold: ", pval_x, "; FC threshold: ", FC_x))
    
    plot_path <- paste0(wd, "www/Plots/pval_", str_replace(pval_x, "\\.",  ","), "_FC_", str_replace(FC_x, "\\.",  ","), "/")
    dir.create(plot_path, recursive = T, showWarnings = F)
    
    filt_degs <- ImportFiltDatasets(paste0(wd, "data/Filtered_DEGs/"), pval_x, FC_x)
    presence_df_filt <- CreateSexDf(filt_degs, unified_annotation)
    all_genes_filt <- do.call(rbind, presence_df_filt)
    all_genes_filt$ct <- gsub("\\..*", "", rownames(all_genes_filt))
    
    # Number of DEGs
    num_degs_filt <- NumDEGsAcrossGroups(presence_df_filt, groups_order)
    num_degs_plot_filt <- PlotNumDEGsFaceted(num_degs_filt, custom_palette)
    png(paste0(plot_path, "Number_of_DEGs.png"), res = 300, units = 'in', height = 22, width = 15)
    print(num_degs_plot_filt)
    dev.off()
    
    # Chr fractions
    gene_counts_filt <- CreateCountDfs(presence_df_filt)
    faceted_chr_fraction <- PlotChrFraction(gene_counts_filt)
    png(paste0(plot_path, "Chr_fractions.png"), res = 300, units = 'in', height = 22, width = 15)
    print(faceted_chr_fraction)
    dev.off()
    
    # Sex-biased gene location counts
    genes_loc <- ExtractLocation(presence_df_filt, location_df, groups_order)
    enriched_loc <- HyperGeomLocation(genes_loc, tot_genes, loc_counts)
    location_results <- merge(genes_loc, enriched_loc, by=c("groups", "ct", "sex", "locations"))
    location_plt <- PlotEnrichedPvalues(location_results, groups_order, cts_order)
    png(paste0(plot_path, "Location.png"), res = 300, units = 'in', height = 20, width = 15)
    print(location_plt)
    dev.off()
    
    # Top 20 most different genes - presence heatmap
    mostdiff <- PlotTop20DiffGenes(all_genes_filt, groups_order)
    png(paste0(plot_path, "top_20_most_diff_genes.png"), res = 300, units = 'in', height = 8, width = 14)
    print(mostdiff)
    dev.off()
    
    # Presence of mitochondrial genes (MT)
    MTgenes <- PlotMTgenes(all_genes_filt, groups_order)
    png(paste0(plot_path, "MT_genes.png"), res = 300, units = 'in', height = 8, width = 14)
    print(MTgenes)
    dev.off()
    
    # Presence of X-escaping genes
    Xescapees <- PlotXescapees(all_genes_filt, groups_order, x_escapees)
    png(paste0(plot_path, "X_escaping_genes.png"), res = 300, units = 'in', height = 8, width = 14)
    print(Xescapees)
    dev.off()
    
    # Percentage of McKenzie cell type markers
    mckenzie_perc <- PlotPercRef(all_genes_filt, McKenzie, groups_order[7:13], McKenzie_ct_names)
    png(paste0(plot_path, "McKenzie_perc.png"), res = 300, units = 'in', height = 8, width = 15)
    print(mckenzie_perc)
    dev.off()
    
    # Chlamydas disease markers
    chl_deg <- CreateDiseaseDf(Chlamydas, all_genes_filt)
    chl_hmp <- PlotDisDegGroup(chl_deg, "Neuropsychiatric diseases", groups_order)
    png(paste0(plot_path, "Chlamydas_hmp.png"), res = 300, units = 'in', height = 8, width = 15)
    print(chl_hmp)
    dev.off()
    
    # Hormone targets enrichment
    hormones_df <- CreateHormonesDf(all_genes_filt, hormones_filt, groups_order)
    hormones_pval <- HormoneEnrichment(hormones_df)
    hmp_hormones <- HmpHormoneEnrichment(hormones_pval, groups_order)
    png(paste0(plot_path, "Hormone_target_enrichment.png"), res = 300, units = 'in', height = 15, width = 10)
    print(hmp_hormones)
    dev.off()
    
    # ARE sites plots
    ARE_plt <- ARE_ERE_plots(all_genes_filt, "ARE", ARE_sites, EREgene, groups_order)
    png(paste0(plot_path, "ARE.png"), res = 300, units = 'in', height = 12, width = 10)
    print(ARE_plt)
    dev.off()
    
    # ERE sites plots
    ERE_plt <- ARE_ERE_plots(all_genes_filt, "ERE", ARE_sites, EREgene, groups_order)
    png(paste0(plot_path, "ERE.png"), res = 300, units = 'in', height = 12, width = 10)
    print(ERE_plt)
    dev.off()
    
    # GO Biological Processes enrichment
    GOBP <- DBClusterProfiler(all_genes_filt[which(all_genes_filt$presence=="Yes"), ], "GO", "BP", gene_thresh = 100, groups_ordered = groups_order, rotate_x_axis = T, adj_pval_thresh =  0.01, cts_sex_order)
    png(paste0(plot_path, "GO_BP.png"), res = 300, units = 'in', height = 30, width = 30)
    print(ggarrange(plotlist = GOBP, ncol = 4, nrow = 5))
    dev.off()
    
    # GO Cellular Component enrichment
    GOCC <- DBClusterProfiler(all_genes_filt[which(all_genes_filt$presence=="Yes"), ], "GO", "CC", gene_thresh = 100, groups_ordered = groups_order, rotate_x_axis = T, adj_pval_thresh =  0.01, cts_sex_order)
    png(paste0(plot_path, "GO_CC.png"), res = 300, units = 'in', height = 30, width = 30)
    print(ggarrange(plotlist = GOCC, ncol = 4, nrow = 5))
    dev.off()
    
    # GO Molecular Function enrichment
    GOMF <- DBClusterProfiler(all_genes_filt[which(all_genes_filt$presence=="Yes"), ], "GO", "MF", gene_thresh = 100, groups_ordered = groups_order, rotate_x_axis = T, adj_pval_thresh =  0.01, cts_sex_order)
    png(paste0(plot_path, "GO_MF.png"), res = 300, units = 'in', height = 30, width = 30)
    print(ggarrange(plotlist = GOMF, ncol = 4, nrow = 5))
    dev.off()
    
    # KEGG enrichment
    KEGG_res <- DBEnrichR(all_genes_filt[which(all_genes_filt$presence=="Yes"), ], "EnrichR",  "KEGG_2021_Human", groups_order, cts_sex_order)
    png(paste0(plot_path, "KEGG.png"), res = 300, units = 'in', height = 30, width = 30)
    print(ggarrange(plotlist = KEGG_res, ncol = 4, nrow = 5))
    dev.off()
    
  }
}

