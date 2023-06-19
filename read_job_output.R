library(tidyverse)
library(ggplot2)
library(openxlsx)
library(kableExtra)
library(stringr)


exp_dir <- "C:/Users/msavery/OneDrive - UGent/Documents/ghent_phd_spatial_doe/data/globus_hpc_collection"
exp_name <- "experimental_runs_full_perm"
result_dir <- file.path(exp_dir, exp_name)
result_dir
file_list <- list.files(result_dir)
file_list

n_f <- length(file_list)
n_f
v_mat <- matrix(nrow=2, ncol=n_f)
v_mat
for (i in 1:length(file_list)){
  f <- file_list[i]
  print(f)
  dir_path <- file.path(result_dir, f)
  file_name <- list.files(dir_path)[2]
  file_path <- file.path(dir_path, file_name)
  print(file_path)
  results <- read.xlsx(file_path, sheet="v_stat")
  # Split strings, get params, and then organize table somehow
  #params <- str_split(f, "_")
  if(str_detect(f, "model-1")==T){
  v_mat[1, i] <- results$v
  }else{
    v_mat[2, i] <- results$v
  }
}
v_mat

# Station id will be indexed correctly to account for station order
split_names <- str_split(tf_preds_s$station_id, "_")
# Get just the station name and the obs/pred label
# Move station id to front
new_names <- c("Station", "RMSE", "MAE", "RMSE", "MAE")
colnames(oos_df) <- new_names
caption <- paste("Out-of-sample errors for", toupper(chem), "for horizon", h, 
                 "when comparing effect of station identification variable.")
label <- paste("station_", chem, "_h-", h, sep="")
print(kbl(oos_df, booktabs = T, escape=F, caption=caption, label=label, 
          align=c('lcccc'), digits=4) %>% 
        kable_styling(latex_options = "HOLD_position") %>%  
        add_header_above(c(" " = 1, "w/ station" = 2, "w/o station" = 2))
)