
# 0. Import Libraries

library(stringr) # to modify and harmonize names
library(ggplot2) # to plot
library(tidyr) # to clean and re-organize dfs
#library(ggpubr) # to assemble plots together before saving
#library(biomaRt) # to query to which chromosome the shared genes belong to
#library(scales) # to set the palette to be used in the PlotDEGsOverlap function
library(RColorBrewer) # to set a palette for the number of DEGs palette

# 1. Import All DEGs from F and M for all ct; slight different folder structure requires different inputs
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
    #names(deg) <- str_remove_all(names(deg), "_filt")
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

# 2. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
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