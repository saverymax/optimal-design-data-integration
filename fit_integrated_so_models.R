library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(cmdstanr) 
library(bayesplot)
library(openxlsx)
library(stringr)
library(readr)
library(kableExtra)

#working_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\code\\optimal_design_site_occ"
#result_path <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\globus_hpc_collection\\exp_5_misspec.1"
#working_dir <- "."
working_dir <- "/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration"
result_path <- "/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/experimental_runs/exp_5"

po_path <- file.path(working_dir, "data", "sim_data", "misspec")

source(file.path(working_dir, "experimental_design_functions.R"))
source(file.path(working_dir, "presence_only_functions.R"))
stan_models_path <- file.path(working_dir, "stan_models", "stan_site_occupancy_models.R")
source(stan_models_path)

# Hardcode all run names, sorry :(
experiment_files <- c(
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.01_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.01_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.01_sequential_no_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.05_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.05_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.05_sequential_no_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.1_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.1_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.1_sequential_no_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0_sequential_no_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.01_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.01_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.01_sequential_no_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.05_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.05_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.05_sequential_no_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.1_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.1_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.1_sequential_no_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0_oracle",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0_sequential_int",
  "oe_model-3_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0_sequential_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.01_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.01_pa-only_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.05_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.05_pa-only_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.1_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0.1_pa-only_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-0.25_p-0.2_misspec-0_pa-only_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.01_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.01_pa-only_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.05_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.05_pa-only_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.1_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0.1_pa-only_no_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0_pa-only_int",
  "oe_model-1_sites-5_n-3_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0_pa-only_no_int"
)

k <- 20
sites <- k^2
m <- 5
n <- 3
design_reps <- 2#10
data_reps <- 2#10
n_surveys <- rep(n, m)
sigma <- 1
visits <- 3
p_0 <- 0.2
deviation <- 5
mcmc_iter <- 500
warmup <- 200
n_chains <- 3
centroid <- c(10,4)
centroid_2 <- c(18,18)
link_func <- "cloglog"

# Data gen params are the same here.
alpha <- -2
beta <- 0.5

model_path <- "misspec_fits_with_po.stan"
model_string <- poisson_process_site_occupancy
write(model_string, model_path)
model_so_po <- cmdstan_model(model_path) 
model_path <- "misspec_fits_no_po.stan"
model_string <- cloglog_site_occupancy
write(model_string, model_path)
model_so <- cmdstan_model(model_path) 



all_params <- c("alpha", "beta")
gen_params <- c("occ_gen")
fit_matrix <- matrix(nrow=length(experiment_files), ncol=length(all_params))
v_matrix <- matrix(nrow=length(experiment_files), ncol=data_reps*design_reps)
alpha_matrix <- matrix(nrow=length(experiment_files), ncol=data_reps*design_reps)
beta_matrix <- matrix(nrow=length(experiment_files), ncol=data_reps*design_reps)
#gamma_matrix <- matrix(nrow=length(experiment_files), ncol=data_reps*design_reps)
#delta_matrix <- matrix(nrow=length(experiment_files), ncol=data_reps*design_reps)
chain_array <- array(dim=c(length(experiment_files), design_reps*data_reps*n_chains*mcmc_iter, length(all_params)))
gen_array <- array(dim=c(length(experiment_files), data_reps*design_reps, sites))
text_size <- 10

sampling_surface <- get_sampling_surface_donut(k, deviation)
sampling_surface <- get_bias_surface_exponential(sampling_surface, centroid)

for (e_i in 1:length(experiment_files)){
  # Get particular strength 
  print(paste("Current experiment: ", experiment_files[e_i]))
  print(e_i)
  param_split <- str_split(experiment_files[e_i], "_")
  misspec <- param_split[[1]][11]
  misspec <- as.numeric(str_split(misspec, "-")[[1]][2])
  delta <- param_split[[1]][9]
  delta <- as.numeric(str_split(delta, "-")[[1]][2])
  # Either oracle, sequential, or pa-only. Don't need the int/no int option
  gen_seq <- param_split[[1]][[12]]
  # Don't need to misspecification to fit the model but good to check it
  sampling_surface <- get_bias_surface_misspecified(sampling_surface, centroid_2, misspec)
  if (misspec == 0){
    sampling_surface$aux_e = rep(0, nrow(sampling_surface))
  }
  # Load PO data generated according to a certain level of misspecfication and set params
  if (!file.exists(file.path(po_path, paste("po_gen_ints-donut_peak=5_a=-2_b=", beta, "_g=1_d=", delta, "_e=", misspec, ".Rds", sep="")))){
    stop(paste("Check data ", po_path, "/po_gen_ints-donut_peak=5_a=-2_b=", beta, "_g=1_d=", delta, "_e=", misspec, ".Rds", sep=""))
  }
  po_gen <- read_rds(file.path(po_path, paste("po_gen_ints-donut_peak=5_a=-2_b=", beta, "_g=1_d=", delta, "_e=", misspec, ".Rds", sep="")))
  if (!file.exists(file.path(result_path, experiment_files[e_i], "results.xlsx"))){
    stop(paste("Check run ", experiment_files[e_i], sep=""))
  }
  design_results <- read.xlsx(file.path(result_path, experiment_files[e_i], "results.xlsx"), sheet="best_sites")

  pa_eval_data <- generate_data_so(data_reps, sampling_surface, corr_matrix=NA, p_0, alpha, beta, sigma, visits, sites, link=link_func)
  
  # Check the data
  p <- ggplot(sampling_surface, aes(x, y, fill=aux_e)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, name="E") +
    ggtitle(paste("Initial sampling surface, E covariate with misspec scale of ", misspec, sep="")) +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  print(p)
  # Bias
  p <- ggplot(sampling_surface, aes(x, y, fill=po_gen$bias[data_reps,])) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, "Bias") +
    ggtitle(paste("Generated bias per site, misspec scale of ", misspec, sep="")) +
    theme(text=element_text(size=text_size), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  print(p)
  # Thinned bias
  thinned_intensity <- po_gen$lambda[data_reps,]*po_gen$bias[data_reps,]
  p <- ggplot() +
    geom_tile(sampling_surface, mapping=aes(x, y, fill=thinned_intensity, width=1, height=1), alpha=.6) + 
    scale_fill_viridis(discrete=FALSE, name="L*b") +
    ggtitle(paste("Generated PP thinning per site, \nincluding observed PO individuals, misspec scale=", misspec, sep="")) +
    geom_point(data=po_gen$Y_coords, mapping=aes(x=x, y=y), size=2, col="white") +
    theme(panel.grid.minor = element_line(colour="white")) +
    scale_y_continuous(breaks = seq(0, 20, 1)) +
    scale_x_continuous(breaks = seq(0, 20, 1)) +
    theme(text=element_text(size=text_size), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  print(p)
  # Iterate over designs and data reps
  # There will be d_i designs for the 10 random starts used in the experiments.
  chain_matrix <- matrix(nrow=0, ncol=length(all_params))
  iter_cnt <- 0
  for (d_i in 1:design_reps){
    d_rep_name <- paste("rand-start", d_i, sep="")
    site_idx <- design_results[[d_rep_name]]
    covars_sampled_sites <- sampling_surface[site_idx, ]
    # Then fit over the data reps
    for (data_rep in 1:data_reps){
      iter_cnt <- iter_cnt + 1
      selected_counts <- pa_eval_data$Y[data_rep, site_idx]
      if (gen_seq=="pa-only"){
        data_site_occ = list(n_surveys=n_surveys, n_sites=m, total_sites=sites, X=covars_sampled_sites$aux_x, X_all=sampling_surface$aux_x, Y=selected_counts)
        fit_so_1 <- quiet(model_so$sample(data=data_site_occ, seed=13, chains=n_chains, iter_sampling=mcmc_iter, iter_warmup=warmup, refresh=0, show_messages=F)) 
      }
      else{
        data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, X=covars_sampled_sites$aux_x, Y=selected_counts, PO=po_gen$Y[1, ], X_po=sampling_surface$aux_x, Z_po=sampling_surface$aux_z, model_diag=0)
        fit_so_1 <- quiet(model_so_po$sample(data=data_site_occ, seed=13, chains=n_chains, iter_sampling=mcmc_iter, iter_warmup=warmup, refresh=0, show_messages=F)) 
      }
      fit_draws <- fit_so_1$draws(all_params, format="matrix")
      if (dim(fit_draws)[1]!=n_chains*mcmc_iter){
        warning(paste("Posterior chains for run", e_i, data_rep, "is of unexpected dimension:", dim(fit_draws), sep=" "))
        # Just a little code so it doesn't crash
        fit_draws <- rbind(fit_draws, matrix(0, nrow=(n_chains*mcmc_iter) - dim(fit_draws)[1], ncol=dim(fit_draws)[2]))
      }
      chain_matrix <- rbind(chain_matrix, fit_draws)
      gen_occupancy <- fit_so_1$summary(variables="occ_gen")$mean
      stopifnot(length(gen_occupancy)==length(pa_eval_data$Y[data_rep,]))
      v_est <- design_criteria(criteria="brier", sim_occ=gen_occupancy, obs_occ=pa_eval_data$Y[data_rep,])
      v_matrix[e_i, iter_cnt] <- v_est
      alpha_matrix[e_i, iter_cnt] <- fit_so_1$summary(variables="alpha")$mean
      beta_matrix[e_i, iter_cnt] <- fit_so_1$summary(variables="beta")$mean
      #gamma_matrix[e_i, iter_cnt] <- fit_so_1$summary(variables="gamma")$mean
      #delta_matrix[e_i, iter_cnt] <- fit_so_1$summary(variables="delta")$mean
      gen_array[e_i, iter_cnt, ] <- gen_occupancy
    }
  }
  fit_matrix[e_i, ] <- colMeans(chain_matrix)
  chain_array[e_i, ,] <- chain_matrix
}

# Get the average optimality scores
v_means <- rowMeans(v_matrix)
a_means <- rowMeans(alpha_matrix)
b_means <- rowMeans(beta_matrix)
#g_means <- rowMeans(gamma_matrix)
#d_means <- rowMeans(delta_matrix)
v_var <- apply(v_matrix, MARGIN=1, FUN=var)
a_var <- apply(alpha_matrix, MARGIN=1, FUN=var)
b_var <- apply(beta_matrix, MARGIN=1, FUN=var)
full_exp_name <- c()
exp_list <- list()
a_exp_list <- list()
b_exp_list <- list()
for (exp in 1:length(experiment_files)){
  param_split <- str_split(experiment_files[exp], "_")
  beta <- param_split[[1]][7]
  delta <- param_split[[1]][9]
  misspec_id <- param_split[[1]][11]
  gamma_run <- param_split[[1]][12:length(param_split[[1]])]
  gamma_run_id <- paste(gamma_run, collapse="_")
  exp_name <- paste(beta, "_", delta, "_", misspec_id, sep="")
  full_exp_name <- c(full_exp_name, paste(exp_name, gamma_run_id))
  run_list <- list()
  a_run_list <- list()
  b_run_list <- list()
  run_list[[gamma_run_id]] <- c(v_means[exp], sqrt(v_var[exp]))
  a_run_list[[gamma_run_id]] <- c(a_means[exp], sqrt(a_var[exp]))
  b_run_list[[gamma_run_id]] <- c(b_means[exp], sqrt(b_var[exp]))
  exp_list[[exp_name]] <- c(exp_list[[exp_name]], run_list)
  a_exp_list[[exp_name]] <- c(a_exp_list[[exp_name]], a_run_list)
  b_exp_list[[exp_name]] <- c(b_exp_list[[exp_name]], b_run_list)
}
print(exp_list)
df <- data.frame(exp=full_exp_name, v_means=v_means, v_sd=sqrt(v_var))
print("Long df of all runs")
print(df)
param_df <- data.frame(exp=full_exp_name, alpha=a_means, alpha_sd=sqrt(a_var), beta=b_means, beta_sd=sqrt(b_var))#, gamma=g_means, delta=d_means)
print("Avg parameter estimates over runs")
print(param_df)

result_matrix <- matrix(nrow=(length(exp_list)), ncol=10)
a_result_matrix <- matrix(nrow=(length(exp_list)), ncol=10)
b_result_matrix <- matrix(nrow=(length(exp_list)), ncol=10)
param_names <- c()
for (i in 1:length(exp_list)){
  result_matrix[i, ] <- unlist(exp_list[[i]])
  a_result_matrix[i, ] <- unlist(a_exp_list[[i]])
  b_result_matrix[i, ] <- unlist(b_exp_list[[i]])
  param_names <- c(param_names, names(exp_list)[i])
}

result_df <- as.data.frame(result_matrix)
a_result_df <- as.data.frame(a_result_matrix)
b_result_df <- as.data.frame(b_result_matrix)

rownames(result_df) <- names(exp_list)
rownames(a_result_df) <- names(exp_list)
rownames(b_result_df) <- names(exp_list)
col_names <- c("Oracle mean", "Oracle sd", "Sequential int mean", "Sequential int sd", "Sequential no int mean", "Sequential no int sd", "PA-only mean", "PA-only sd")
colnames(result_df) <- col_names
colnames(a_result_df) <- col_names
colnames(b_result_df) <- col_names

caption <- paste("Comparison of models fit to designs")
label <- paste("survey_eval", sep="")
print(kbl(result_df, booktabs = T, escape=T, caption=caption, label=label, 
          align=c('lcccccccc'), digits=4, format="latex") %>% 
        kable_styling(latex_options = c("HOLD_position")) )

caption <- paste("Comparison of alpha estimates from designs")
label <- paste("alpha_analysis", sep="")
print(kbl(a_result_df, booktabs = T, escape=T, caption=caption, label=label, 
          align=c('lcccccccc'), digits=4, format="latex") %>% 
        kable_styling(latex_options = c("HOLD_position")) )

caption <- paste("Comparison of beta estimtaes from designs")
label <- paste("beta_analysis", sep="")
print(kbl(b_result_df, booktabs = T, escape=T, caption=caption, label=label, 
          align=c('lcccccccc'), digits=4, format="latex") %>% 
        kable_styling(latex_options = c("HOLD_position")) )

#A little too intense to look at all PPDs for $Y$.
#for (i in 1:length(experiment_files)){
#  for (j in 1:data_reps){
#    print(experiment_files[i])
#    generated_occ <- gen_array[i, j, ]
#    ppd_df <- data.frame(x=sampling_surface$x, y=sampling_surface$y, occ_rep=generated_occ)
#    
#    # Plot generated ys and thetas on the sampling grid
#    p <- ggplot(ppd_df, aes(x, y, fill=occ_rep)) + 
#        geom_tile() +
#        scale_fill_viridis(discrete=FALSE) +
#        ggtitle("PPD occupancy per site")
#    print(p)
#  }
#}
