
# 0. Import Libraries

#library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs
#library(ggpubr) # to assemble plots together before saving
#library(biomaRt) # to query to which chromosome the shared genes belong to
#library(scales) # to set the palette to be used in the PlotDEGsOverlap function
library(RColorBrewer) # to set a palette for the number of DEGs palette

# 1. Function to filter out ns genes and too low FC, and order based on FC 
  # Input: dataframe of DEGs
  # Return: gene list of significant genes as data.frame

Filter_gene <- function( order.gene.df, pval, FC) {
  logFC <- log2(FC)
  gene.sig <- order.gene.df[  order.gene.df[["p_val_adj"]] <= pval
                              & order.gene.df[["avg_log2FC"]] >= logFC, ]
  return(data.frame("Genes"=rownames(gene.sig)))
}

# 2. Function to filter the non-signifcant DEGs in all datasets
  # Input: list of DEGs, pvalue and FC thresholds
  # Return: list of significant DEGs 

FilterDs <- function(list_ds, pval, FC) {
  filt_ds <- list()
  for (group_id in names(list_ds)) {
    sex_ls <- list("F" = lapply(1:length(names(list_ds[[group_id]][["F"]])), function(x) Filter_gene(list_ds[[group_id]][["F"]][[x]], pval, FC)),
                   "M" = lapply(1:length(names(list_ds[[group_id]][["M"]])), function(x) Filter_gene(list_ds[[group_id]][["M"]][[x]], pval, FC)))
    names(sex_ls[["F"]]) <- names(list_ds[[group_id]][["F"]])
    names(sex_ls[["M"]]) <- names(list_ds[[group_id]][["M"]])
    filt_ds <- append(filt_ds, list(sex_ls))
  }
  names(filt_ds) <- names(list_ds)
  return(filt_ds)
}

# 3. Creates the df for the input ct so that we know if a DEG is found in a certain group or not -> used to generate hmps
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

# 4. Creates all PresenceDfs for all cts
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

# 5. Groups cts according to common annotation, then creates the presence dfs
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

# 6. Count DEGs for each ct in each age
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

# 7. Plots the total number of DEGs across groups, faceted by ct and sex
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