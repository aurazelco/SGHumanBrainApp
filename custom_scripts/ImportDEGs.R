
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
  df_F <- list()
  df_M <- list()
  names_F <- vector()
  names_M <- vector()
  for (ct in 1:length(sub_ct)) {
    deg <- ImportDEGs(paste(path, sub_ct[ct], sep="/"), ext, row_col)
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

# 3. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: main directory, file extension, where to find row names
  # Return: unfiltered DEGs as nexted list

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

# 4. Imports DISCO and UCSC datasets; slight different folder structure requires different inputs
  # Input: unfiltered DEGs as nexted list, the common annotation to compare cell types
  # Return: unfiltered DEGs as nexted list with new annotation

UnifyAnnotation <- function(ds_list, common_annotation) {
  ds_common_annot <- list()
  for (group in names(ds_list)) {
    group_ls <- list()
    for (sex in names(ds_list[[group]])) {
      sex_ls <- list()
      new_names <- vector()
      for (ct in unique(common_annotation)) {
        common_ct <- names(common_annotation[which(common_annotation==ct)])
        if(any(common_ct %in% tolower(names(ds_list[[group]][[sex]])))) {
          ct_ls <- ds_list[[group]][[sex]][names(ds_list[[group]][[sex]])[which(tolower(names(ds_list[[group]][[sex]])) %in% common_ct)]]
          ct_ls <- do.call(rbind, ct_ls)
          sex_ls <- append(sex_ls, list(ct_ls))
          new_names <- c(new_names, ct)
        }
      }
      names(sex_ls) <- new_names
      group_ls <- append(group_ls, list(sex_ls))
    }
    names(group_ls) <- names(ds_list[[group]])
    ds_common_annot <- append(ds_common_annot, list(group_ls))
  }
  names(ds_common_annot) <- names(ds_list)
  return(ds_common_annot)
}