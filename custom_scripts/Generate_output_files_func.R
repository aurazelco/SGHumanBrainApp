# 0. Import libraries
library(readxl) # to import Excel files
library(biomaRt) # to query to which chromosome the shared genes belong to
library(scales) # to set the palette to be used in the PlotDEGsOverlap function
library(stringr)  # to modify and harmonize names
library(RColorBrewer)  # to set a palette for the number of DEGs palette
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs
library(rjson) # to import the json files containing the genes associated with the corresponding hormone


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
    group_names <- vector()
    for (ct_id in names(list_ds[[group_id]])) {
      ct_ls <- list("F" = Filter_gene(list_ds[[group_id]][[ct_id]][["F"]], pval, FC),
                    "M" = Filter_gene(list_ds[[group_id]][[ct_id]][["M"]], pval, FC))
      ct_ls <- ct_ls[lapply(ct_ls,length)>0]
      if (length(ct_ls)==2) {
        group_ls <- append(group_ls, list(ct_ls)) 
        group_names <- c(group_names, ct_id)
      }
    }
    names(group_ls) <- group_names
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


# 19. Calculates the genes whihc are most expressed in one sex compared to the other
  # Input: all genes with presence info, whihc sex to analyze, the pther sex to compare with
  # Return: df with the genes ordered by descending number of difference between sex and the other sex

CalculateMostDiffGenes <- function(genes_all_presence, sex, other_sex) {
  sex_count <- as.data.frame(table(genes_all_presence[which(genes_all_presence$presence=="Yes" & genes_all_presence$sex==sex), "gene_id"]))
  sex_count$sex <- rep(sex, nrow(sex_count))
  other_sex_count <- as.data.frame(table(genes_all_presence[which((genes_all_presence$gene_id %in% unique(sex_count$Var1)) & genes_all_presence$sex==other_sex & genes_all_presence$presence=="Yes"), "gene_id"]))
  other_sex_count$sex <- rep(other_sex, nrow(other_sex_count))
  sex_count <- rbind(sex_count, other_sex_count)
  colnames(sex_count) <- c("gene_id", "count", "sex")
  sex_count <- complete(sex_count, gene_id, sex)
  sex_count[which(is.na(sex_count$count)), "count"] <- 0
  
  abs_diff_sex <- data.frame("gene_id"=unique(sex_count$gene_id))
  abs_diff_sex$sex_diff <- rep(NA, nrow(abs_diff_sex))
  for (i in abs_diff_sex$gene_id) {
    abs_diff_sex[which(abs_diff_sex$gene_id==i), "sex_diff"] <- sex_count[which(sex_count$sex==sex & sex_count$gene_id==i), "count"] - sex_count[which(sex_count$sex==other_sex & sex_count$gene_id==i), "count"]
  }
  abs_diff_sex <- abs_diff_sex[which(abs_diff_sex$sex_diff > 0 ),]
  abs_diff_sex <- abs_diff_sex[order(abs_diff_sex$sex_diff, decreasing = T), ]
  abs_diff_sex$gene_id <- factor(abs_diff_sex$gene_id, unique(abs_diff_sex$gene_id))
  return(abs_diff_sex)
}

# 20. Extracts the top 20 genes (10 per sex) that are most different in presence
  # Input: the abs diff df for F, for M, all genes with presence info
  # Return: df with the top 20 most different genes

ExtractTop20DiffGenes <- function(abs_diff_F, abs_diff_M, genes_all_presence) {
  top10F <- as.character(abs_diff_F[1:10, "gene_id"])
  top10M <- as.character(abs_diff_M[1:10, "gene_id"])
  common_genes_sex <- intersect(top10F, top10M)
  `%!in%` <- Negate(`%in%`)
  
  if (length(common_genes_sex) > 0 ) {
    top10F_u <- top10F[which(top10F %!in% common_genes_sex)]
    top10M_u <- top10M[which(top10M %!in% common_genes_sex)]
    while (length(top10F_u) < 10 & length(top10M_u) < 10) {
      for (i  in common_genes) {
        if (abs_diff_F[which(abs_diff_F$gene_id==i), "sex_diff"] > abs_diff_M[which(abs_diff_M$gene_id==i), "sex_diff"]) {
          top10F_u <- c(top10F_u, i)
        } else {
          top10M_u <- c(top10M_u, i)
        }
      }
    }
  } else {
    top10F_u <- top10F
    top10M_u <- top10M
  }
  
  if (length(top10F_u) < 10) {
    missing_genes <- as.character(abs_diff_F$gene_id[11:(11+10 - length(top10F_u) - 1)])
    top10F_u <- c(top10F_u, missing_genes)
  }
  if (length(top10M_u) < 10) {
    missing_genes <- as.character(abs_diff_M$gene_id[11:(11+10 - length(top10M_u) - 1)])
    top10M_u <- c(top10M_u, missing_genes)
  }
  
  most_diff_genes <- c(top10F_u, top10M_u)
  most_diff_genes <- complete(genes_all_presence[which(genes_all_presence$gene_id %in% most_diff_genes), ], gene_id, groups, sex, ct)
  most_diff_genes$gene_id <- factor(most_diff_genes$gene_id, rev(c(top10F_u, top10M_u)))
  return(most_diff_genes)
}

# 21. Calculates the most different genes and plots the presence heatmap
  # Input: the presence unlisted df, the order of the groups
  # Return: presence heatmap for top 20 most different genes

PlotTop20DiffGenes <- function(genes_all_presence, groups_ordered) {
  abs_diff_F <- CalculateMostDiffGenes(genes_all_presence, "F", "M")
  abs_diff_M <- CalculateMostDiffGenes(genes_all_presence, "M", "F")
  most_diff_genes_df <- ExtractTop20DiffGenes(abs_diff_F, abs_diff_M, genes_all_presence)
  mostdiff_plot <- ggplot(most_diff_genes_df,
                          aes(factor(groups, groups_ordered[which(groups_ordered %in% unique(groups))]), gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(sex ~  ct, scales = "free") +
    labs(x="Groups", y="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text = element_text(size = 8, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(mostdiff_plot)
}

# 22. Presence heatmap for the mitochondrial genes
  # Input: the presence unlisted df, the order of the groups
  # Return: mitochondrial genes heatmap

PlotMTgenes <- function(genes_all_presence, groups_ordered) {
  mit_genes_ids <- c(
                     unique(genes_all_presence$gene_id[which(grepl("^TIMM", genes_all_presence$gene_id))]),
                     unique(genes_all_presence$gene_id[which(grepl("^TOMM", genes_all_presence$gene_id))]),
                     unique(genes_all_presence$gene_id[which(grepl("^MT-", genes_all_presence$gene_id))])
                     )
  mit_genes_ids <- c("XIST", mit_genes_ids)
  
  mit_genes <- genes_all_presence[which(genes_all_presence$gene_id %in% mit_genes_ids), ]
  mit_gene_count <- as.data.frame(table(mit_genes[which(mit_genes$presence=="Yes"), "gene_id"]))
  mit_gene_count <- mit_gene_count[order(mit_gene_count$Freq, decreasing = T),]
  mit_genes <- complete(mit_genes, gene_id, groups,sex,ct)
  MTplot <- 
    ggplot(mit_genes, 
           aes(factor(groups, groups_ordered[which(groups_ordered %in% unique(groups))]), factor(gene_id, rev(unique(mit_genes_ids))), fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(sex ~ ct , scales = "free") +
    labs(x="Groups", y="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text = element_text(size = 8, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(MTplot)
}

# 23. Presence heatmap of X-escaping genes
  # Input: the presence unlisted df, the order of the groups
  # Return: presence heatmap fo X-escaping genes

PlotXescapees <- function(genes_all_presence, groups_ordered, x_escapees_df) {
  Xescaping_plot <- 
  ggplot(complete(genes_all_presence[which(genes_all_presence$gene_id %in% x_escapees_df$V2), ], gene_id, groups, sex, ct), 
         aes(factor(groups, groups_ordered[which(groups_ordered %in% unique(groups))]), gene_id, fill=presence)) +
    geom_tile() +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "#00BFC4",
                      guide = guide_legend(reverse = TRUE)) +
    facet_grid(sex ~ ct , scales = "free") +
    labs(x="Groups", y="Genes", fill="Genes found") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text = element_text(size = 8, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(Xescaping_plot)
}

# 24. Calculates the % of known markers in the DEGs
  # Input: the presence df, the ct to plot, the gene lists
  # Return: the percent df

RefPerc <- function(ref_presence_df, ref_ct_id, sex_df, plot_titles) {
  pos_markers <- ref_presence_df[which(ref_presence_df$ref_ct==ref_ct_id), ]
  tot_genes <- vector()
  tot_names <- vector()
  num_pos <- vector()
  for (id in unique(pos_markers$groups)) {
    for (ct in unique(pos_markers[which(pos_markers$groups==id), "group_id_ct"])) {
      tot_names <- c(tot_names, paste(id, ct, "F", sep = "/"), paste(id, ct, "M", sep = "/"))
      num_pos <- c(num_pos, length(pos_markers[which(pos_markers$groups==id & pos_markers$group_id_ct==ct & pos_markers$sex=="F" & pos_markers$presence=="Yes"), "gene_ids"]))
      num_pos <- c(num_pos, length(pos_markers[which(pos_markers$groups==id & pos_markers$group_id_ct==ct & pos_markers$sex=="M" & pos_markers$presence=="Yes"), "gene_ids"]))
      tot_genes <- c(tot_genes, length(sex_df[which(sex_df$groups==id & sex_df$ct==ct & sex_df$sex=="F"), "gene_id"]))
      tot_genes <- c(tot_genes, length(sex_df[which(sex_df$groups==id & sex_df$ct==ct & sex_df$sex=="M"), "gene_id"]))
    }
  }
  ref_perc <- data.frame("ref_ct"=rep(plot_titles[ref_ct_id], length(tot_names)), tot_names, num_pos, tot_genes)
  ref_perc <- separate(ref_perc, tot_names, into = c("groups", "ct", "sex"), sep = "/")
  ref_perc$perc <- ref_perc$num_pos * 100 / ref_perc$tot_genes
  return(ref_perc)
}

# 25. Plot the presence number of genes as percentage
# Input: the presence df, the groups order
# Return: the plot

PlotBarPlotRefPerc <- function(ref_perc, groups_ordered) {
  ref_plot <- ggplot(ref_perc, aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), perc, fill = ref_ct)) +
    geom_bar(stat="identity", color="black", position = "dodge") +
    facet_grid(sex ~ ct, scales = "free") +
    labs(x="Groups", y="Markers %", fill="Reference cell types") +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          strip.text = element_text(size = 8, face="bold", colour = "black"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black"),
          axis.ticks.y = element_blank(),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(ref_plot)
}

# 26. Plots if gebes from a reference df are found or not in the DEGs
  # Input: the dataframe containing all DEGs, the reference df, the order in which plot the groups 
    # (and which groups to plot), the vector to use for plot titles
  # Return: plot

PlotPercRef <- function(sex_df, ref_df, groups_ordered, plot_titles){
  sex_df_filt <- sex_df[which(sex_df$groups %in% groups_ordered), ]
  presence <- vector()
  ids <- vector()
  gene_ids <- vector()
  for (sex_id in c("F", "M")) {
    for (ct in unique(ref_df$Celltype)) {
      ref_genes <- ref_df[which(ref_df$Celltype==ct), "gene"]
      for (group_id in unique(sex_df_filt[which(sex_df_filt$sex==sex_id), "groups"])) {
        for (ct_id in unique(sex_df_filt[which(sex_df_filt$sex==sex_id & sex_df_filt$groups==group_id), "ct"])) {
          presence <- c(presence, 
                        ifelse(ref_df[which(ref_df$Celltype==ct), "gene"] %in% sex_df_filt[which(sex_df_filt$sex==sex_id & sex_df_filt$groups==group_id & sex_df_filt$ct==ct_id), "gene_id"],
                               "Yes", "No"))
          ids <- c(ids, 
                   rep(paste(sex_id, ct, group_id, ct_id, sep = "/"), length(ref_genes)))
          gene_ids <- c(gene_ids, ref_genes)
          
        }
      }
    }
  }
  ref_presence_df <- data.frame(ids, gene_ids, presence)
  ref_presence_df <- separate(ref_presence_df, ids, into=c("sex", "ref_ct", "groups", "group_id_ct"), sep = "/")
  ref_presence_df$groups <- factor(ref_presence_df$groups, groups_ordered[which(groups_ordered %in% unique(ref_presence_df$groups))])
  ref_presence_df <- ref_presence_df[order(ref_presence_df$groups), ]
  ref_perc <- lapply(1:length(unique(ref_presence_df$ref_ct)), function(x) RefPerc(ref_presence_df, unique(ref_presence_df$ref_ct)[x], sex_df, plot_titles))
  ref_perc <- do.call(rbind, ref_perc)
  ref_plot <- PlotBarPlotRefPerc(ref_perc, groups_ordered)
  return(ref_plot)
}

# 27. Creates a dataframe with only the gens related to know diseases
  # Input:the reference dataframe with the diseases and genes, 
    # one dataframe containing all DEGs, the reference name to be used for the output folder
  # Return: the DEG dataframe with the presence of genes-associated genes and saves the results to CSV file

CreateDiseaseDf <- function(ref, sex_dfs) {
  dis_genes <- vector()
  group_id <- vector()
  for (dis_family in unique(ref$Disease_group)) {
    for (dis in unique(ref[which(ref$Disease_group==dis_family), "Disease"])) {
      dis_genes <- c(dis_genes, unlist(str_split(ref[which(ref$Disease_group==dis_family & ref$Disease==dis), "Affected_gene"], pattern = ", ")))
      group_id <- c(group_id, rep(paste(dis_family, dis, sep = "/"), length(unlist(str_split(ref[which(ref$Disease_group==dis_family & ref$Disease==dis), "Affected_gene"], pattern = ", ")))))
    }
  }
  ref_df <- data.frame(group_id, dis_genes)
  sex_dfs$id <- paste(sex_dfs$groups, sex_dfs$sex, sex_dfs$ct, sep = "/")
  dis_presence <- vector()
  dis_names <- vector()
  deg_ids <- vector()
  genes_ids <- vector()
  for (id in unique(ref_df$group_id)) {
    for (deg in unique(sex_dfs$id)) {
      genes_ids <- c(genes_ids, ref_df[which(ref_df$group_id==id), "dis_genes"])
      deg_presence <- ifelse(ref_df[which(ref_df$group_id==id), "dis_genes"] %in% sex_dfs[which(sex_dfs$id==deg), "gene_id"], "Yes", "No")
      dis_presence <- c(dis_presence, deg_presence)
      dis_names <-c(dis_names, rep(id, length(deg_presence)))
      deg_ids <- c(deg_ids, rep(deg, length(deg_presence)))
    }
  }
  ref_deg <- data.frame(dis_names, deg_ids, genes_ids, dis_presence)
  ref_deg <- separate(ref_deg, dis_names, into = c("disease_group", "disease"), sep = "/")
  ref_deg <- separate(ref_deg, deg_ids, into = c("groups", "sex", "ct"), sep = "/")
  ref_deg$dis_gene_id <- paste(ref_deg$disease, ref_deg$genes_ids, sep =  " - ")
  return(ref_deg)
}

# 28. Plot heatmap with the results of which disease-associated genes are found in the degs
  # Iput: the DEG data frame with the presence of genes-associated genes, the disease group to plot
  # Return: the faceted plot

PlotDisDegGroup <- function(ref_deg, dis_id, groups_ordered) {
  dis_plot <- ggplot(complete(ref_deg[which(ref_deg$disease_group==dis_id),]), aes(factor(groups, groups_ordered[which(groups_ordered %in% groups)]), dis_gene_id, fill=dis_presence)) +
    geom_tile(color="white") +
    facet_grid(sex ~ ct, scales = "free") +
    labs(x="Groups", y="Disease-associated genes", fill="Genes found", title =dis_id) +
    scale_fill_manual(values = c("Yes"="#F8766D",
                                 "No"="#00BFC4"),
                      na.value = "grey",
                      guide = guide_legend(reverse = TRUE)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          panel.spacing.x=unit(0, "lines"),
          plot.title = element_text(size=12, face="bold", colour = "black"),
          strip.text = element_text(size = 8, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_text(size=12, face="bold", colour = "black"),
          axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
          axis.ticks.x=element_blank(),
          axis.title.y = element_text(size=12, face="bold", colour = "black"),
          axis.text.y = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5),
          axis.ticks.y = element_blank(),
          legend.position = "bottom", 
          legend.title = element_text(size=12, face="bold", colour = "black"))
  return(dis_plot)
}

# 29. Generates a df with the nunmber of found hormone targets as absoluet numbers, percetages fo the hormone gene lists and percentages of  the sex-biased DEGs
  # Input: one dataframe with all sex-biased DEGs, the hormone gene lists in a list, and the order of the groups
  # Return: dataframe with all the counts and percentages

CreateHormonesDf <- function(sex_dfs, ref_hormones, groups_ordered) {
  hormone_ls <- list()
  for (horm in names(ref_hormones)) {
    group_ids <- vector()
    hormone_tgs <- vector()
    bg_genes <- vector()
    for (group_id in unique(sex_dfs$groups)) {
      for (ct in unique(sex_dfs[which(sex_dfs$groups==group_id), "ct"])) {
        for (sex in c("F", "M")) {
          group_ids <- c(group_ids, paste(group_id, ct, sex, sep = "/"))
          hormone_tgs <- c(hormone_tgs, 
                           length(intersect(
                             ref_hormones[[horm]], 
                             tolower(unique(sex_dfs[which(sex_dfs$groups==group_id & sex_dfs$ct==ct &sex_dfs$sex==sex & sex_dfs$presence=="Yes"), "gene_id"])))))
          bg_genes <- c(bg_genes, length(unique(sex_dfs[which(sex_dfs$groups==group_id & sex_dfs$ct==ct &sex_dfs$sex==sex & sex_dfs$presence=="Yes"), "gene_id"])))
        }
      }
    }
    tot_horm <- rep(length(ref_hormones[[horm]]), length(group_ids))
    horm_df <- data.frame(group_ids, hormone_tgs, tot_horm, bg_genes)
    horm_df$no_tgs <- horm_df$bg_genes - horm_df$hormone_tgs
    hormone_ls <- append(hormone_ls, list(horm_df))
  }
  names(hormone_ls) <- str_to_title(str_replace_all(names(ref_hormones), "/", "_"))
  hormones_df <- do.call(rbind, hormone_ls)
  hormones_df <- cbind("hormones"=gsub("\\..*", "", rownames(hormones_df)), hormones_df)
  rownames(hormones_df) <- NULL
  hormones_df <- separate(hormones_df, group_ids, into = c("groups", "ct", "sex"), sep = "/", remove = T)
  hormones_df$perc_hormones <- hormones_df$hormone_tgs * 100 / hormones_df$tot_horm
  hormones_df$perc_degs <- hormones_df$hormone_tgs * 100 / hormones_df$bg_genes
  return(hormones_df)
}

# 30. Calculates the enrichment of hormones in each sub group with a hypergeometric test
  # Input: the dataframe with all the counts and percentages, the significant pvalue threshold, the minimum number of groups to keep
  # Return: dataframe with all the significant pvalues

HormoneEnrichment <- function(hormones_df, pval_thresh=0.05, min_num_cond=1) {
  pval_names <- vector()
  pvalues <- vector()
  for (horm in unique(hormones_df$hormones)) {
    for (ct in unique(hormones_df[which(hormones_df$hormones==horm), "ct"])) {
      for (cond in unique(hormones_df[which(hormones_df$hormones==horm & hormones_df$ct==ct), "groups"])) {
        for (sex in c("F", "M")) {
          pval <- phyper(
            as.numeric(hormones_df[which(hormones_df$hormones==horm & hormones_df$ct==ct & hormones_df$groups==cond & hormones_df$sex==sex), "hormone_tgs"]) - 1,
            as.numeric(hormones_df[which(hormones_df$hormones==horm & hormones_df$ct==ct & hormones_df$groups==cond & hormones_df$sex==sex), "bg_genes"]),
            10000 - as.numeric(unique(hormones_df[which(hormones_df$hormones==horm), "tot_horm"])),
            as.numeric(unique(hormones_df[which(hormones_df$hormones==horm), "tot_horm"])),
            lower.tail= FALSE
          )
          pvalues <- c(pvalues, pval)
          pval_names <- c(pval_names, paste(horm, ct, cond, sex, sep="/"))
        }
      }
    }
  }
  pval_df <- data.frame(pval_names, pvalues)
  pval_df <- separate(pval_df, pval_names, into = c("hormone_id", "ct", "groups", "sex"), sep = "/", remove = T)
  pval_df<- pval_df[which(pval_df$pvalues < pval_thresh), ]
  pval_df <- pval_df %>% group_by(hormone_id) %>% filter(n() > min_num_cond)
  return(pval_df)
}

# 31. Plots the results of the hormone analysis as one faceted plot
  # Input: dataframe with all the significant pvalues, the order of the groups
  # Return: p-value heatmap

HmpHormoneEnrichment <- function(pval_df, groups_ordered) {
  hormone_hmp <- 
    ggplot(pval_df, aes(factor(groups, groups_ordered[which(groups_ordered %in% unique(groups))]), ct, fill=pvalues)) +
      geom_tile() +
      facet_grid(hormone_id ~ sex, scales = "free") +
      scale_fill_gradient(low="red", high="blue", na.value = "grey") +
      labs(x="Groups", y="Cell types", fill="P-values") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            strip.text.x = element_text(size = 8, face="bold", colour = "black"),
            strip.text.y = element_text(size=8, face="bold", colour = "black",angle = 0),
            axis.title.x = element_text(size=12, face="bold", colour = "black"),
            axis.text.x = element_text(size=8, colour = "black", vjust = 0.7, hjust=0.5, angle = 90),
            axis.ticks.x=element_blank(),
            axis.title.y = element_text(size=12, face="bold", colour = "black"),
            legend.position = "bottom", 
            legend.title = element_text(size=12, face="bold", colour = "black"))
  return(hormone_hmp)
}