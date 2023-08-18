# Little test script with some of my data and parallel library in R.
library(parallel)
library(spatstat)
library(cmdstanr) 
n_cores <- detectCores()
n_cores
clust <- makeCluster(n_cores)

df <- data.frame(c(1,2,3), c(1,2,3))
results <- parApply(clust, df, 2, FUN=mean)
results

# Then use my data
source("experimental_design_functions.R")
source("presence_only_functions.R")
source("stan_models\\stan_site_occupancy_models.R")
data_reps <- 10
# Area for whole space, which allows us to set area for sites based on number of sites.
area_D <- 100
# For generating data according to GP
gp_bool <- F
# k is one side of grid
k <- 20
sites <- k^2
alpha <- -2 
beta <- 2
gamma <- -2
delta <- 0.5
p_0 <- 0.7
sigma <- 1
# number of surveys at site is equal to n or 0.
n_surveys <- 5
# There will be m sites selected for sampling
# 36/4 was used in paper
m <- 5
# Based on the size of grid get the auxiliary data and coordinates
# Because these don't really depend on any random variables and just location 
# in the grid, there will be one fixed dataset throughout the optimization
sampling_surface <- get_sampling_surface_simple(k)
sampling_surface
# Next, we use this data to generate the rest of the datasets
# The data generating function will sample R occupancy maps|params
# and then R complete datasets|occupancy maps
corr_matrix <- specify_corr(sampling_surface[,1:2])
link_func <- "cloglog"
r_survey_data <- generate_data_so(data_reps, sampling_surface, corr_matrix, p_0, alpha, beta, sigma, n_surveys, sites, link=link_func)
r_survey_data$occupancy[1,]
r_survey_data$Y[1,]
r_survey_data$theta[1,]

# Then generate R presence-only datasets
params <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta)
# Provide occupancy maps and number of data reps, as well as params and sampling surface, to generate pp data.
#r_po_data <- generate_ppp_data_r(sampling_surface, params, sites, r_survey_data$occupancy, data_reps)
r_po_data <- generate_ppp_data_r(sampling_surface, params, sites, data_reps, corr_matrix, gp_bool, area_D)
Y_positive_indices <- which(r_po_data$Y>0)
# This is the data at which there are counts > 0
r_po_data$lambda
r_po_data$Y[Y_positive_indices]


dim(r_survey_data$occupancy)
dim(r_survey_data$Y)
dim(r_po_data$Y)
combined_df <- cbind(r_survey_data$occupancy, r_survey_data$Y, r_po_data$Y)
combined_df[1, ((2*sites)+1):ncol(combined_df)] 
length(combined_df[((2*sites)+1):ncol(combined_df)])

# Each row will be one repetition
combined_df[,1]
dim(combined_df)

# Apply function over each row
result <- parApply(clust, combined_df, 1, FUN=mean)
length(result)
result

#########################
# Now add in a stan model...
#########################

model_path <- "stan_models\\poisson_process_prior_site_occupancy.stan"
# We can experiment with these models in the exchange algorithm
model_strings <- c(cloglog_site_occupancy, site_occupany_detection, poisson_process_site_occupancy, pp_site_occ_no_aux)
model_selection <- 3
write(model_strings[3], model_path)
model <- cmdstan_model(model_path) 
# These are the parameters to report, though this will be model dependent
if (model_selection==3){
  params <- c('p', 'alpha', 'beta', 'gamma', 'delta')
}else if (model_selection==4){
  params <- c('p', 'lambda')
}else{
  params <- c('p', 'alpha', 'beta')
}
generated_vars <- c('g_theta_gen', 'occ_gen')
site_idx <- sample(1:400, m, replace=F)
site_idx
select_sites <- sampling_surface[site_idx, ]
select_sites

combined_df[1, site_idx]
selected_data <- combined_df[1, sites+site_idx]
PO_data <- combined_df[1, (2*sites+1):ncol(combined_df)]

data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, X=select_sites$aux_x, Y=selected_data, PO=PO_data, 
                     X_po=sampling_surface$aux_x, Z_po=sampling_surface$aux_z, model_diag=0)
data_site_occ$X

fit <- model$sample(data=data_site_occ, seed=13, chains=1, 
                    iter_sampling=1000, iter_warmup=100, refresh=0, show_messages=F)

estimate_v_parallel(combined_df[1,], model, n_surveys, m, sites, sampling_surface,
                             site_idx, select_sites, params, generated_vars, model_selection, 1000)

# Export the cluster
clusterExport(clust, varlist=c("design_criteria", "brier_score_stan"), envir=environment())
estimate_mat <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel, model, n_surveys, m, sites, sampling_surface, 
                             site_idx, select_sites, params, generated_vars, model_selection, 1000)
estimate_mat


stopCluster(clust)