# 0. Import libraries
library(biomaRt) # to query to which chromosome the shared genes belong to
library(scales) # to set the palette to be used in the PlotDEGsOverlap function
library(stringr)  # to modify and harmonize names
library(RColorBrewer)  # to set a palette for the number of DEGs palette
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs

# 1. Import data for each ct
  # Input: the path where to find the DEGs for each cell type within a dataset, file extension, where to find row names
  # Return: list of dfs

ImportDEGs <- function(path, ext, row_col) {
  deg_files <- list.files(path = path, pattern = paste0("\\.",ext,"$"),full.names = TRUE)
  deg <- lapply(deg_files, read.csv, row.names=row_col)
  names_deg <- list.files(path = path, pattern = "\\.csv$",full.names = FALSE)
  names(deg) <- substr(names_deg, 1, nchar(names_deg)-4)
  return(deg)
}

# 2. Import All DEGs from F and M for all ct; slight different folder structure requires different inputs
  # Input: directory where to find ct sub-folders, file extension, where to find row names
  # Return: list of 2 lists, one for F and one for M dfs

ImportCts <- function(main_dir, ext, row_col) {
  path <- paste0(main_dir, "/Unfiltered_DEGs/")
  sub_ct <- list.dirs(path, recursive=F, full.names = F)
  ds_degs <- list()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDEGs(paste(path, sub_ct[ct], sep="/"), ext, row_col)
    ds_degs <- append(ds_degs, list(deg))}
  names(ds_degs) <- sub_ct
  return(ds_degs)
}


# 3. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, file extension, where to find row names
  # Return: unfiltered DEGs as nested list

ImportDatasets <- function(main_dir, ext="csv", row_col=1) {
  ds_list <- list()
  group_names <- vector()
  folder_list <- list.dirs(main_dir, recursive = F, full.names = F)
  for (folder in folder_list) {
    ds_list <- append(ds_list, list(ImportCts(paste0(main_dir, folder), ext, row_col)))
    group_names <- c(group_names, folder)
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}

# 4. Import All DEGs from F and M for all ct; slight different folder structure requires different inputs
  # Input: directory where to find ct sub-folders, file extension, where to find row names
  # Return: list of 2 lists, one for F and one for M dfs

ImportFiltDEGs <- function(main_dir, ext, row_col) {
  sub_ct <- list.dirs(main_dir, recursive=F, full.names = F)
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDEGs(paste(main_dir, sub_ct[ct], sep="/"), ext, row_col)
    for (i in names(deg)) {
      if (grepl("F", i, fixed=TRUE)){
        df_F <- append(df_F, list(deg[[i]]))
        names_F <- c(names_F, sub_ct[ct])
      } else {
        df_M <- append(df_M, list(deg[[i]]))
        names_M <- c(names_M, sub_ct[ct])
      }
    }
  }
  names(df_F) <- names_F
  names(df_M) <- names_M
  return(list("F"=df_F, "M"=df_M))
}

# 5. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, file extension, where to find row names
  # Return: unfiltered DEGs as nested list

ImportFiltDatasets <- function(main_dir, pval_thresh, FC_thresh, ext="csv", row_col=1) {
  ds_list <- list()
  group_names <- vector()
  filt_folder <- paste0(main_dir, "pval_", str_replace(pval_thresh, "\\.",  ","), "_FC_", str_replace(FC_thresh, "\\.",  ","), "/")
  folder_list <- list.dirs(filt_folder, recursive = F, full.names = F)
  for (folder in folder_list) {
    ds_list <- append(ds_list, list(ImportFiltDEGs(paste0(filt_folder, folder), ext, row_col)))
    group_names <- c(group_names, folder)
  }
  names(ds_list) <- group_names
  ds_list[lengths(ds_list) != 0]
  return(ds_list)
}


# 6. Function to filter out ns genes and too low FC, and order based on FC 
  # Input: dataframe of DEGs
  # Return: gene list of significant genes as data.frame

Filter_gene <- function( order.gene.df, pval, FC) {
  logFC <- log2(FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val_adj"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  if(nrow(gene.sig) > 10) {
    gene.sig$index <- seq.int(nrow(gene.sig))
    gene.sig <- data.frame("Genes"=rownames(gene.sig))
  } else {
    gene.sig <- data.frame()
  }
  return(gene.sig)
}

# 7. Function to filter the non-signifcant DEGs in all datasets
  # Input: list of DEGs, pvalue and FC thresholds
  # Return: list of significant DEGs 

FilterDs <- function(list_ds, pval, FC) {
  filt_ds <- list()
  for (group_id in names(list_ds)) {
    group_ls <- list()
    for (ct_id in names(list_ds[[group_id]])) {
      ct_ls <- list("F" = Filter_gene(list_ds[[group_id]][[ct_id]][["F"]], pval, FC),
                    "M" = Filter_gene(list_ds[[group_id]][[ct_id]][["M"]], pval, FC))
      ct_ls <- ct_ls[lapply(ct_ls,length)>0]
      group_ls <- append(group_ls, list(ct_ls))
    }
    names(group_ls) <- names(list_ds[[group_id]])
    group_ls <- group_ls[lapply(group_ls,length)>0]
    filt_ds <- append(filt_ds, list(group_ls))
  }
  names(filt_ds) <- names(list_ds)
  filt_ds <- filt_ds[lapply(filt_ds,length)>0]
  return(filt_ds)
}

# 8. Creates the df for the input ct so that we know if a DEG is found in a certain group or not -> used to generate hmps
  # Input: list of ct dfs, which sex and ct to analyze
  # Return: df with info whether each gene is present in all group groups in which the ct is found

CreatePresenceCtDf <- function(sex_dfs, sex, ct) {
  sub_ct <- sex_dfs[[sex]][which(sex_dfs[[sex]]$common_annot==ct), ]
  ct_sex <-(rep(unique(sub_ct$gene_id), length(unique(sub_ct$groups))))
  ct_sex <- cbind(as.data.frame(ct_sex), rep(unique(sub_ct$groups), each=length(unique(sub_ct$gene_id))))
  colnames(ct_sex) <- c("gene_id", "groups")
  ct_sex$sex <- rep(sex, nrow(ct_sex))
  ct_sex$presence <- rep("No", nrow(ct_sex))
  for (group_id in unique(ct_sex$groups)) {
    ct_sex[which(ct_sex$groups==group_id), "presence"] <- ifelse(unique(ct_sex$gene_id) %in% sub_ct[which(sub_ct$groups==group_id), "gene_id"], "Yes", "No")
  }
  return(ct_sex)
}

# 9. Creates all PresenceDfs for all cts
  # Input: list of lists, each list corresponding to a specific group-sex-ct combination
  # Return: list of ct dfs, with information on presence of each gene across all groups

CreatePresenceDf <- function(sex_dfs) {
  ct_df_list <- list()
  for (ct in unique(sex_dfs[["F"]]$common_annot)) {
    f_ct <- CreatePresenceCtDf(sex_dfs, "F", ct)
    m_ct <- CreatePresenceCtDf(sex_dfs, "M", ct)
    df_ct <- rbind(f_ct, m_ct)
    ct_df_list <- append(ct_df_list, list(df_ct))    
  }
  names(ct_df_list) <- unique(sex_dfs[["F"]]$common_annot)
  return(ct_df_list)
}

# 10. Groups cts according to common annotation, then creates the presence dfs
  # Input: list of lists generated from ImportDEGs, and the named vector used to harmonize the annotation
  # Return: list of presence dfs, one per each ct

CreateSexDf <- function(list_ds, common_annot) {
  all <- unlist(list_ds, recursive = F)
  sex_dfs <- list()
  for (sex in c("F", "M")) {
    sex_list <- unlist(all[names(all)[which(grepl(paste0("\\.", sex), names(all)))]])
    sex_ct <- data.frame()
    for (ct in unique(common_annot)) {
      sex_filt <- list()
      for (ct_class in names(common_annot[which(common_annot==ct)])) {
        sex_filt <- append(sex_filt, sex_list[names(sex_list)[which(grepl(ct_class, tolower(names(sex_list))))]])
      }
      if (length(sex_filt)>0) {
        sex_df <- data.frame("groups" = rep(names(sex_filt), sapply(sex_filt, length)),
                             "gene_id" = unlist(sex_filt))
        rownames(sex_df) <- NULL
        sex_df <- separate(sex_df, groups, into=c("groups", "sex", "ct", "gene_num"), sep="\\.", remove=T)
        sex_df$gene_num <- NULL
        sex_df$common_annot <- rep(ct, nrow(sex_df))
        sex_ct <- rbind(sex_ct, sex_df)
      }
    } 
    sex_dfs <- append(sex_dfs, list(sex_ct))
  }
  names(sex_dfs) <- c("F", "M")
  ct_df_list <- CreatePresenceDf(sex_dfs)
  return(ct_df_list)
}

# 11. Count DEGs for each ct in each age
  # Input: list of presence dfs, one per each ct, order of the group
  # Return: dataframe with num of DEGs for each ct and group

NumDEGsAcrossGroups <- function(ct_df_list, groups_ordered) {
  ct <- vector()
  groups <- vector()
  count_degs <- vector()
  for (ct_id in names(ct_df_list)) {
    sub_ct <- ct_df_list[[ct_id]]
    for (group_id in groups_ordered) {
      ct <- c(ct, rep(ct_id, 2))
      groups <- c(groups, rep(group_id, 2))
      if (group_id %in% unique(sub_ct$groups)) {
        count_degs <- c(count_degs, nrow(sub_ct[which(sub_ct$groups==group_id & sub_ct$sex=="F" & sub_ct$presence=="Yes"), ]))
        count_degs <- c(count_degs, nrow(sub_ct[which(sub_ct$groups==group_id & sub_ct$sex=="M" & sub_ct$presence=="Yes"), ]))
      } else {
        count_degs <- c(count_degs, rep(0, 2))
      }
    }
  }
  sex <- rep(c("F", "M"), length(groups)/2)
  num_deg_df <- data.frame(ct, groups, sex, count_degs)
  num_deg_df$groups <- factor(num_deg_df$groups, groups_ordered[which(groups_ordered %in% unique(num_deg_df$groups))])
  num_deg_df <- num_deg_df[order(num_deg_df$groups), ]
  return(num_deg_df)
}

# 12. Plots the total number of DEGs across groups, faceted by ct and sex
  # Input: the dataframe with num of DEGs for each ct and group, the palette to use
  # Return: nothing, saves the plot instead

PlotNumDEGsFaceted <- function(num_deg, col_palette) {
  num_degs_plot <- ggplot(num_deg, aes(groups, count_degs, fill=groups)) +
    geom_bar(stat = "identity", show.legend = T, color="black") +
    labs(x="Datasets", y="Number of DEGs", fill="Datasets") +
    scale_fill_manual(values = col_palette) + 
    facet_grid(ct ~ sex, scales = "free", switch = "y", drop = T) +
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = NA, color = "black"), 
      panel.spacing.x = unit(0.5, "lines"),
      plot.title = element_text(size=16, face="bold", colour = "black"),
      axis.line = element_line(colour = "black"),
      axis.title.y = element_text(size=16, face="bold", colour = "black"),
      axis.text.y = element_text(size=14, colour = "black", vjust = 0.7, hjust=0.5),
      axis.title.x = element_text(size=16, face="bold", colour = "black"),
      axis.text.x = element_blank(), 
      axis.ticks.x = element_blank(), 
      strip.text = element_text(size = 10, face="bold", colour = "black"),
      legend.position = "bottom", 
      legend.text = element_text(size=12, colour = "black"),
      legend.title = element_text(size=10, face="bold", colour = "black"))
  return(num_degs_plot)
}

# 13. Creates dfs which counts in how many groups we find each gene, per sex and ct combo
  # Input: ct df obtained previously, and the sex to analyze
  # Return: df containing for each gene the number of groups which had that gene in their DEGs

GroupsSharingGenes <- function(ct_df, sex_id) {
  ct_df <- ct_df[which(ct_df$sex==sex_id),]
  gene_id <- vector()
  groups_count <- vector()
  for (id_gene in unique(ct_df$gene_id)) {
    gene_id <- c(gene_id, id_gene)
    groups_count <- c(groups_count, length(unique(ct_df[which(ct_df$gene_id==id_gene & ct_df$presence=="Yes"), "groups"])))
  }
  sex <- rep(sex_id, length(gene_id))
  return(data.frame(gene_id, sex, groups_count))
}

# 14. Create Count Dfs for all cts
  # Input: presence df list
  # Return: list of df containing the number of groups for each gene, for each ct

CreateCountDfs <- function(ct_df_list) {
  gene_count_dfs <- list()
  for (ct in names(ct_df_list)) {
    f_df <- GroupsSharingGenes(ct_df_list[[ct]], "F")
    m_df <- GroupsSharingGenes(ct_df_list[[ct]], "M")
    gene_count_dfs <- append(gene_count_dfs, list(rbind(f_df, m_df)))
  }
  names(gene_count_dfs) <- names(ct_df_list)
  return(gene_count_dfs)
}

# 15. Function to get chromosome number from gene symbol
  # Input: the genes as vector
  # Return: the annotated genes

Annot.chr.name <- function(gene.list){
  # define biomart object
  mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "www")
  Annot_idf <- getBM(attributes = c("hgnc_symbol",
                                    "chromosome_name"),
                     filters = c("hgnc_symbol") ,
                     values=list(gene.list),
                     mart = mart)
  #delete chromosome name with CHR label
  Annot_df <- Annot_idf[!str_detect(Annot_idf$chromosome_name,  "CHR"),]
  return(Annot_df)
}

# 16. Map genes from intersected genes against chromosome
  # Input: the list of filtered genes, the annotated genes
  # Return: merged dataframe

map_chr <- function(gene_count_filt, Annot_df){
  map_chr_df <- merge(gene_count_filt, Annot_df, by.x= "gene_id", by.y= "hgnc_symbol")
  return(map_chr_df)
}


# 17. Plot fractions of shared chromosomes
  # Input: shared gene df, the color palette to use
  # Return: bar plot with fractions of chr genes shared in how many groups

PlotNumSharedGenesChr <- function(shared_genes_chr, col_palette) {
  chr_plot <- ggplot(shared_genes_chr, aes(groups_count, fill=chr_simplified)) +
    geom_bar(position = "stack", color="black") +
    scale_y_log10() +
    facet_grid(ct ~ sex, scales = "free") +
    scale_fill_manual(values = c("X" = col_palette[1],
                                 "Y"= col_palette[3],
                                 "Autosome"= col_palette[2])) +
    scale_x_continuous(breaks=seq(min(shared_genes_chr$groups_count), max(shared_genes_chr$groups_count), by=1)) +
    labs(x = "Number of groups sharing genes", y="Log10 counts of shared genes", fill="Chromosomes") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text = element_text(size = 12, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(chr_plot)
}

# 18. Plot all count dfs for all cts
  # Input: main directory where to save the plots, the list of count dfs
  # Return: bar plot with fractions of chr genes shared in how many groups

PlotChrFraction <- function(gene_count_dfs) {
  chr_results <- lapply(gene_count_dfs, function(x) Annot.chr.name(x$gene_id))
  shared_genes_chr <- list()
  for(i in names(gene_count_dfs)){
    shared_genes_ct <- map_chr(gene_count_dfs[[i]], chr_results[[i]])
    shared_genes_chr <- append(shared_genes_chr, list(shared_genes_ct))
  }
  names(shared_genes_chr) <- names(gene_count_dfs)
  shared_genes_chr <- do.call(rbind, shared_genes_chr)
  shared_genes_chr$chr_simplified <- shared_genes_chr$chromosome_name
  shared_genes_chr$chr_simplified[which(shared_genes_chr$chr_simplified!= "X" & shared_genes_chr$chr_simplified!= "Y")] <- "Autosome"
  shared_genes_chr <- cbind("ct" = gsub('\\..*', '', rownames(shared_genes_chr)), shared_genes_chr)
  rownames(shared_genes_chr) <- NULL
  col_palette <- hue_pal()(3)
  chr_fractions <- PlotNumSharedGenesChr(shared_genes_chr, col_palette)
  return(chr_fractions)
}

