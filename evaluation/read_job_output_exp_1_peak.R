library(tidyverse)
library(ggplot2)
library(openxlsx)
library(kableExtra)
library(stringr)


#exp_dir <- "C:/Users/msavery/OneDrive - UGent/Documents/ghent_phd_spatial_doe/data/globus_hpc_collection"
exp_dir <- "/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration"

exp_name <- "experimental_runs/exp_1_peak"
result_dir <- file.path(exp_dir, exp_name)
print(result_dir)
# Look at all files corresponding to this set of experiments
file_list <- list.files(result_dir)
print("available runs")
print(file_list)
n_f <- length(file_list)
print("number of files")
print(n_f)
#exp_col <- c("model-1_n-5", "model-1_n-10", "model-3_n-5", "model-3_n-10")
#exp_col <- c("model-1_n-1", "model-1_n-5", "model-3_n-1", "model-3_n-5", "model-4_n-1", "model-4_n-5", "model-5_n-1", "model-5_n-5")
exp_col <- c("model-1_n-1", "model-1_n-5", "model-3_n-1", "model-3_n-5")
#exp_col <- c("model-1_n-5", "model-1_n-10", "model-3_n-5", "model-3_n-10", "model-4_n-5", "model-4_n-10", "model-5_n-5", "model-5_n-10")
exp_list <- vector("list", length=length(exp_col))
perm_list <- vector("list", length=length(exp_col))
names(exp_list) <- exp_col
names(perm_list) <- exp_col
print("inital list to save results")
print(exp_list)
for (i in 1:length(file_list)){
  f <- file_list[i]
  print(f)
  dir_path <- file.path(result_dir, f)
  file_name <- list.files(dir_path)[2]
  params <- str_split(f, "_")
  print(params)
  compare_name <- paste(params[[1]][2], params[[1]][4], sep="_")
  print("Name to identify run")
  print(compare_name)
  param_perm <- paste(params[[1]][9], params[[1]][10], sep="_")
  print("current params of interest")
  print(param_perm)
  # Make this list so as to have the permutation names for each subset
  perm_list[[compare_name]] <- c(perm_list[[compare_name]], param_perm)
  if (substr(file_name, str_length(file_name)-3, str_length(file_name))!= "xlsx"){
    print(paste("Incorrect file selected:", file_name))
    print(paste("Available files:", list.files(dir_path)))
    exp_list[[compare_name]] <- c(exp_list[[compare_name]], NA)
  } 
  else{
    #if (substr(file_name, 1, 1) == "r"){
    #  new_file <- substr(file_name, 12, str_length(file_name))
    #  file_path <- file.path(dir_path, new_file)
    #  old_path <- file.path(dir_path, file_name)
    #  file.rename(old_path, file_path)
    #}else{
    #  file_path <- file.path(dir_path, file_name)
    #}
    file_path <- file.path(dir_path, file_name) 
    print("reading file:")
    print(file_path)
    results <- read.xlsx(file_path, sheet="v_stat")
    # Split strings, get params, and then organize table somehow
    exp_list[[compare_name]] <- c(exp_list[[compare_name]], results$v)
  }
}

print("result list")
print(exp_list)
print("Params per model. Make sure each model has params ordered in the same way")
print(perm_list)
for (i in 1:length(exp_col)){
  stopifnot(identical(perm_list[[1]], perm_list[[i]]))
}
v_mat <- matrix(nrow=length(exp_col), ncol=n_f/length(exp_col))
print(dim(v_mat))
for (i in 1:length(exp_col)){
  v_mat[i,] <- exp_list[[exp_col[i]]]
  #v_mat[,i] <- exp_list[[exp_col[i]]]
}
print(v_mat)
v_df <- as.data.frame(v_mat)
rownames(v_df) <- exp_col
colnames(v_df) <- perm_list[[1]]
print(v_df)
caption <- paste("Comparison of models, sampling effort, and parameter permutations")
label <- paste("survey_eval", sep="")
print(kbl(v_df, booktabs = T, escape=T, caption=caption, label=label, 
          align=c('lcccc'), digits=4, format="latex") %>% 
        kable_styling(latex_options = c("HOLD_position")) )
#%>%  
#        add_header_above(c(" " = 1, "w/ station" = 2, "w/o station" = 2)))
