###################################


######################### Generate Filtered DEGs
wd <- "/Users/aurazelco/Desktop/Lund_MSc/Thesis/HumanBrainSexSingleCell/"
source(paste0(wd, "custom_scripts/Generate_output_files_func.R"))


# Imports the unfiltered DEGs and generates another list, plus defines common annotation and group order to be followed in plots

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

all_degs <- ImportDatasets(paste0(wd, "data/Unfiltered_DEGs/"))
x_escapees <- read.table(paste0(wd, "data/extra_files/escape_Xchr.txt"), sep="\t", skip = 2)


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

pval_ls <- c(1, 0.05, 0.01, 0.001)
FC_ls <- c(1, 1.2, 1.5, 2)

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

######################### Plots

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
    #num_degs_filt <- NumDEGsAcrossGroups(presence_df_filt, groups_order)
    #num_degs_plot_filt <- PlotNumDEGsFaceted(num_degs_filt, custom_palette)
    #png(paste0(plot_path, "Number_of_DEGs.png"), res = 300, units = 'in', height = 22, width = 15)
    #print(num_degs_plot_filt)
    #dev.off()
    
    # Chr fractions
    #gene_counts_filt <- CreateCountDfs(presence_df_filt)
    #faceted_chr_fraction <- PlotChrFraction(gene_counts_filt)
    #png(paste0(plot_path, "Chr_fractions.png"), res = 300, units = 'in', height = 22, width = 15)
    #print(faceted_chr_fraction)
    #dev.off()
    
    # Top 20 most different genes - presence heatmap
    #mostdiff <- PlotTop20DiffGenes(all_genes_filt, groups_order)
    #png(paste0(plot_path, "top_20_most_diff_genes.png"), res = 300, units = 'in', height = 8, width = 14)
    #print(mostdiff)
    #dev.off()
    
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
    
  }
}


