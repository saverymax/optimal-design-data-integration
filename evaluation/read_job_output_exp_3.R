library(tidyverse)
library(ggplot2)
library(openxlsx)
library(kableExtra)
library(stringr)


#exp_dir <- "C:/Users/msavery/OneDrive - UGent/Documents/ghent_phd_spatial_doe/data/globus_hpc_collection"
exp_dir <- "/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/experimental_runs"

exp_name <- "exp_3"
result_dir <- file.path(exp_dir, exp_name)
print(result_dir)
# Look at all files corresponding to this set of experiments
file_list <- list.files(result_dir)
print("available runs")
print(file_list)
n_f <- length(file_list)
print("number of files")
print(n_f)
exp_col <- c("model-1_effort-n-2", "model-1_effort-n-3", "model-1_effort-n-4", "model-1_effort-n-5", "model-1_effort-n-7", "model-1_effort-n-10", "model-3_effort-n-2", "model-3_effort-n-3", "model-3_effort-n-4", "model-3_effort-n-5", "model-3_effort-n-7", "model-3_effort-n-10")
#exp_col <- c("model-1_n-5", "model-1_n-10", "model-3_n-5", "model-3_n-10", "model-4_n-5", "model-4_n-10", "model-5_n-5", "model-5_n-10")
exp_list <- vector("list", length=length(exp_col))
var_list <- vector("list", length=length(exp_col))
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
  param_perm <- paste(params[[1]][7], params[[1]][8], sep="_")
  print("current params of interest")
  print(param_perm)
  # Make this list so as to have the permutation names for each subset
  perm_list[[compare_name]] <- c(perm_list[[compare_name]], param_perm)
  if (substr(file_name, str_length(file_name)-3, str_length(file_name))!= "xlsx"){
    print(paste("Incorrect file selected:", file_name))
    print(paste("Available files:", list.files(dir_path)))
    exp_list[[compare_name]] <- c(exp_list[[compare_name]], NA)
    var_list[[compare_name]] <- c(var_list[[compare_name]], NA)
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
    var_list[[compare_name]] <- c(var_list[[compare_name]], results$v_var)
  }
}

print("result list")
print(exp_list)
print("variance")
print(var_list)
print("Params per model. Make sure each model has params ordered in the same way")
print(perm_list)
for (i in 1:length(exp_col)){
  stopifnot(identical(perm_list[[1]], perm_list[[i]]))
}
v_mat <- matrix(nrow=length(exp_col), ncol=n_f/length(exp_col))
var_mat <- matrix(nrow=length(exp_col), ncol=n_f/length(exp_col))
print(dim(v_mat))
for (i in 1:length(exp_col)){
  v_mat[i,] <- exp_list[[exp_col[i]]]
  var_mat[i,] <- var_list[[exp_col[i]]]
  #v_mat[,i] <- exp_list[[exp_col[i]]]
}
print(v_mat)
print(var_mat)
v_df <- as.data.frame(v_mat)
var_df <- as.data.frame(var_mat)
rownames(v_df) <- exp_col
rownames(var_df) <- exp_col
colnames(v_df) <- perm_list[[1]]
colnames(var_df) <- perm_list[[1]]
print(v_df)
print(var_df)
caption <- paste("Comparison of models, sampling effort, and parameter permutations")
label <- paste("survey_eval", sep="")
print(kbl(v_df, booktabs = T, escape=T, caption=caption, label=label, 
          align=c('lcccc'), digits=4, format="latex") %>% 
        kable_styling(latex_options = c("HOLD_position")) )
#%>%  

# For this we can also make a figure as effort increases
v_df
var_df
v_df$run <- rownames(v_df)
var_df$run <- rownames(var_df)
m1_1 <- v_df$`g-1_d-0.25`[1:6]
m2_1 <- v_df$`g-1_d-0.25`[7:12]
m1_2 <- v_df$`g-1_d-2`[1:6]
m2_2 <- v_df$`g-1_d-2`[7:12]
var_m1_1 <- var_df$`g-1_d-0.25`[1:6]
var_m2_1 <- var_df$`g-1_d-0.25`[7:12]
var_m1_2 <- var_df$`g-1_d-2`[1:6]
var_m2_2 <- var_df$`g-1_d-2`[7:12]
fig_df <- data.frame(m1_g1_d25=m1_1, m2_g1_d25=m2_1, m1_g1_d2=m1_2, m2_g1_d2=m2_2)
var_fig_df <- data.frame(var_m1_g1_d25=var_m1_1, var_m2_g1_d25=var_m2_1, var_m1_g1_d2=var_m1_2, var_m2_g1_d2=var_m2_2)
fig_df$index <- rownames(fig_df)
var_fig_df$index <- rownames(var_fig_df)
fig_df <- pivot_longer(fig_df, cols=1:4, names_to="run", values_to="avg_v")
var_fig_df <- pivot_longer(var_fig_df, cols=1:4, names_to="run", values_to="var")
fig_df$x <- rep(c(2,3,4,5,7,10), each=4)
sqrt(var_fig_df$var)
fig_df$se <- sqrt(var_fig_df$var)/sqrt(10)

# Generally will be run from evaluation directory
fig_name <- file.path(".", "exp_3_comparison.png")
p <- ggplot(data=fig_df, aes(x=x, y=avg_v, colour=run)) +
  geom_line(linewidth=1) +
  #geom_errorbar(aes(ymin=avg_v-se, ymax=avg_v+se)) +
  labs(title="", x="Max visits", y="U(d)") + 
  theme(text=element_text(size=7), axis.title = element_text(size = 7), legend.key.size = unit(0.25, 'cm')) +
  scale_color_discrete(name = "Run", type=c("#c356ea","#ffc100", "#71aef2", "#f7adce"), labels = c("SO, d=0.25", "SO, d=2", "SO+PO, d=0.25", "SO+PO, d=2")) +
  theme_bw()
print(p)
ggsave(fig_name, plot=p, dpi=300, width=10, height=7, units="cm")
