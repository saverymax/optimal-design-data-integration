###########################
# Analagous script for optimal design in .qmd file
########################## 

rm(list = ls())

library(ggplot2)
library(viridis)
library(hrbrthemes)
library(cmdstanr) 
library(bayesplot)
library(dplyr)
library(tidyr)
library(reshape2)
library(spatstat)
library(parallel)

source("experimental_design_functions.R")
source("presence_only_functions.R")
source("stan_models\\stan_site_occupancy_models.R")

# R=1000 datasets for monte carlo approx
exp_args <- list(model_selection=3, m=5, data_reps=8, random_starts=3, p_logging=F, 
                 mcmc_iter=1000, intensity_func="simple", v_parallel=T, exch_iter=20)
fig_dir <- paste("figures\\optimal_design_model-", exp_args$model_selection,  "_m-", exp_args$m, "_r-", exp_args$data_reps, "\\", sep="")
dir.create(fig_dir)
data_reps <- exp_args$data_reps
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
#b_0 <- 0
# b_01=2 indicates "good" quality of auxiliary information
#b_1 <- 2
# This assumes spatial variance of 1, which was used in the paper (see supplement)
sigma <- 1
# number of surveys at site is equal to n or 0.
n_surveys <- 5
# There will be m sites selected for sampling
# 36/4 was used in paper
m <- exp_args$m
# Based on the size of grid get the auxiliary data and coordinates
# Because these don't really depend on any random variables and just location 
# in the grid, there will be one fixed dataset throughout the optimization
if (exp_args$intensity_func == "simple"){
  sampling_surface <- get_sampling_surface_simple(k)
}else{
  sampling_surface <- get_sampling_surface(k)
  # Filter for only a quarter of the grid.
  # If we filter, we need to change the total sites as well
  sites <- sites/4
  stopifnot(sites>m)
  sampling_surface <- sampling_surface %>% dplyr::filter(x<11, y<11)
}
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
r_po_data$Y[Y_positive_indices]
# So that gives us R complete datasets for 400 sites.
# These will stay fixed throughout the rest of the procedure

# For the coordinate exchange algorithm, we will need to precompute the nearest neighbors.
# I take a naive approach here of choosing the top l neighbors
l <- 4
nearest_neighbors <- get_neighbors(sampling_surface, l)
dim(nearest_neighbors)
nearest_neighbors

# Use the final generation iteration to look at the presence-absence and presence-only data
survey_data_df <- data.frame(counts=r_survey_data$Y[1,], o=r_survey_data$occupancy[1,], 
                             x=sampling_surface$x, y=sampling_surface$y, theta=r_survey_data$theta[1,], cor_mat=as.vector(corr_matrix[1,]))
head(survey_data_df)
dim(survey_data_df)

# First PA data
p <- ggplot(sampling_surface, aes(x, y, fill=r_survey_data$Y[data_reps,])) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("Generated counts per site")
print(p)

p <- ggplot(sampling_surface, aes(x, y, fill=r_survey_data$occupancy[data_reps,])) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("Generated occupancy per site")
print(p)

p <- ggplot(sampling_surface, aes(x, y, fill=r_survey_data$theta[data_reps,])) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("generated occupancy probability per site")
print(p)

# Then the PO data
p <- ggplot(sampling_surface, aes(x, y, fill=r_po_data$lambda[data_reps,])) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("Generated intensity per site")
print(p)

p <- ggplot(sampling_surface, aes(x, y, fill=r_po_data$bias[data_reps,])) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("Generated bias per site")
print(p)

thinned_intensity <- r_po_data$lambda[data_reps,]*r_po_data$bias[data_reps,]
p <- ggplot() +
  geom_tile(sampling_surface, mapping=aes(x, y, fill=thinned_intensity, width=1, height=1), alpha=.6) + 
  scale_fill_viridis(discrete=FALSE, name="L*b") +
  ggtitle("Generated thinning per site, including generated (or observed) individuals") +
  geom_point(data=r_po_data$Y_coords, mapping=aes(x=x, y=y), size=3, col="white") +
  theme(panel.grid.minor = element_line(colour="white")) +
  scale_y_continuous(breaks = seq(0, 20, 1)) +
  scale_x_continuous(breaks = seq(0, 20, 1)) 
print(p)

model_path <- "stan_models\\poisson_process_prior_site_occupancy.stan"
# We can experiment with these models in the exchange algorithm
model_strings <- c(cloglog_site_occupancy, site_occupany_detection, poisson_process_site_occupancy, pp_site_occ_no_aux)
model_selection <- exp_args$model_selection
write(model_strings[model_selection], model_path)
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

# The reich paper repeats the entire exchange algorithm procedure 10 times,
# and retains solution with lowest V(D)
# Currently I am not repeating the procedure. However, it accounts for uncertainty associated with random m sites 
# selected for sampling
p_logging <- exp_args$p_logging
random_starts <- exp_args$random_starts
v_list <- vector(mode="list", length=random_starts)
names(v_list) <- c(1:random_starts)
v_list
best_v <- vector(mode="numeric", length=random_starts)
best_site_mat <- matrix(nrow=random_starts, ncol=m)
best_site_mat
# Initiate parallel processing if specified
if (exp_args$v_parallel==T){
  n_cores <- detectCores()
  print(paste("Using parallel processing with", n_cores, "cores"))
  clust <- makeCluster(n_cores)
  # Export the environment to the cluster
  clusterExport(clust, varlist=c("design_criteria", "brier_score_stan"), envir=environment())
}
# Start timer
start_time <- Sys.time()
for (r_start in 1:random_starts){
  print(paste("random initialization ", r_start, sep=""))
  # we first randomly select m sites, for each run of the exchange algorithm,
  # from the sampling surface, which includes the auxiliary information.
  # These sites will have n_i = n, the others will have n_i = 0
  # Only sites with n_i=n will contribute to likelihood for the site-occupancy model.
  site_idx <- sample(1:sites, m, replace=F)
  site_idx
  # Initialize for exchange algorithm
  best_neighbor_idx <- site_idx
  select_sites <- sampling_surface[site_idx,]
  select_sites
  # Convergence condition will be met where full iteration through all sampled sites results in no change
  # in sites. ie no changes in sites can improve criterion
  convergence_cond <- FALSE
  exchange_iter <- 0
  v_vec <- c()
  # Compute v for initial design
  if (exp_args$v_parallel==T){
    combined_df <- cbind(r_survey_data$occupancy, r_survey_data$Y, r_po_data$Y)
    estimate_vec <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel, model, n_surveys, m, sites, sampling_surface, 
                                 site_idx, select_sites, params, generated_vars, model_selection, exp_args$mcmc_iter)
    # Compatible format with non-parallel v
    estimate_mat <- matrix(estimate_vec, nrow=data_reps, ncol=1)
  }
  else{
    estimate_mat <- estimate_v(model, n_surveys, data_reps, m, sites, sampling_surface, 
                               site_idx, select_sites, r_survey_data, r_po_data, r,
                               p_logging, params, generated_vars, model_selection, exp_args$mcmc_iter)
  }
  # Once the posterior is computed on each of R datasets, find the average score:
  new_v_est <- sum(estimate_mat) / data_reps
  print("Design score from most initial exchange")
  print(new_v_est)
  v_vec <- c(v_vec, new_v_est)
  current_v_est <- new_v_est
  print("Initial row ids")
  print(site_idx)
  print("Initial coordinates")
  print(sampling_surface[site_idx,1:2])
  print("Initial design score")
  print(new_v_est)
  title <- paste("Inital spatial design: v=", round(new_v_est, 10), sep="")
  p <- plot_sites(sampling_surface, site_idx, title)
  fig_name <- paste(fig_dir, "initial_design.png", sep="")
  ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
  
  while ((convergence_cond==FALSE) & (exchange_iter<exp_args$exch_iter)){
    exchange_iter <- exchange_iter + 1
    print(paste("New exchange iteration: ", exchange_iter))
    # Data structure for each score estimate
    # Then compute the posterior based on those sites and generated data, for each r dataset
    # We iterate through the sites, computing the estimate of $V(D)$ for each so that we explore the effect of each site on the design
    for (s in 1:length(site_idx)){
      print(paste("ex iter: ", exchange_iter, ", current site index: ", s, sep=""))
      current_site <- site_idx[s]
      print(paste("current site: ", current_site, sep=""))
      neighbor_set <- nearest_neighbors[current_site,]
      n_count <- 0
      # We first iterate through the neighbors of each site, and compute V for each exchange. 
      # The inital estimate will be for our initial design.
      # We set the current best set of indices
      #best_neighbor_idx <- site_idx
      for(nn in neighbor_set){
        n_count <- n_count + 1
        # Select new data using the neighbors (switch out local points)
        neighbor_idx <- exchange_coordinates_deterministic(best_neighbor_idx, s, nn)
        select_sites <- sampling_surface[neighbor_idx,] 
        if (exp_args$v_parallel==T){
          combined_df <- cbind(r_survey_data$occupancy, r_survey_data$Y, r_po_data$Y)
          estimate_vec <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel, model, n_surveys, m, sites, sampling_surface, 
                                       neighbor_idx, select_sites, params, generated_vars, model_selection, exp_args$mcmc_iter)
          # Compatible format with non-parallel v
          estimate_mat <- matrix(estimate_vec, nrow=data_reps, ncol=1)
        }
        else{
          estimate_mat <- estimate_v(model, n_surveys, data_reps, m, sites, sampling_surface, 
                                     neighbor_idx, select_sites, r_survey_data, r_po_data, r,
                                     p_logging, params, generated_vars, model_selection, exp_args$mcmc_iter)
        }
        # Once the posterior is computed on each of R datasets, find the average score:
        new_v_est <- sum(estimate_mat) / data_reps
        # TODO: Write this to file
        #print("Design score from most recent exchange")
        #print(new_v_est)
        v_vec <- c(v_vec, new_v_est)
        
        if (new_v_est < current_v_est){
           current_v_est <- new_v_est
           title <- paste("New optimal spatial design: v=", round(new_v_est, 10), sep="")
           # Plot the new best site compared to the previous selection, but need to reverse arguments to function
           p <- plot_sites_vs_best(sampling_surface, current_site, best_neighbor_idx, neighbor_idx, title)
           # Then set new best indices
           #site_idx <- select_idx
           best_neighbor_idx <- neighbor_idx
           print("New optimal row ids")
           print(best_neighbor_idx) 
           print("New coordinates")
           print(sampling_surface[best_neighbor_idx,1:2])
           print("New optimal design score")
           print(current_v_est)
        }
        else{
          title <- paste("Non-optimal spatial design: v=", round(new_v_est, 10), 
                         "\nvs current optimal design: v=", round(current_v_est, 10), sep="")
          #print("No change in optimal design")
          p <- plot_sites_vs_best(sampling_surface, current_site, neighbor_idx, best_neighbor_idx, title)
        }
        fig_name <- paste(fig_dir, "site_locs_rand-start-", r_start, "_ex-iter_", 
                          exchange_iter, "_site-iter-", s, "_nn-iter", n_count,".png", sep="")
        # TODO: Decrease legend size
        ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
        # Then go to the next neighbor or site
      }
      # Set design to best from iteration through neighbors of one site
      # If there is no change from any neighbors, site_idx will not change
      site_idx <- best_neighbor_idx
    }
    # If after a complete iteration through all the sites, the best sites haven't changed
    # then we can call that convergence
    # If it's the first iteration we need to initialize the best sites
    if(exchange_iter==1){
      best_iter_idx <- site_idx
    }
    else if(all(best_iter_idx==site_idx)){
      convergence_cond <- TRUE
    }
    else{
      print("Previous best sites and current best sites")
      print(best_iter_idx)
      print(site_idx)
      best_iter_idx <- site_idx
    }
  }
  v_list[[r_start]] <- v_vec
  best_site_mat[r_start,] <- best_iter_idx
  best_v[r_start] <- current_v_est
}
end_time <- Sys.time()
run_time <- end_time - start_time
run_time

v_list
best_site_mat

# Plot the convergence of V(D)
v_concat <- c()
x_concat <- c()
rep_labels <- c()
for(i in 1:random_starts){
  v_concat <- c(v_concat, v_list[[i]])
  l_vec <- rep(paste("y", i, sep=""), length(v_list[[i]]))
  print(length(l_vec))
  x_concat <- c(x_concat, 1:length(l_vec))
  rep_labels <- c(rep_labels, l_vec)
}
  length(rep_labels)
length(v_concat)
v_df <- data.frame(x=x_concat, v=v_concat, r=rep_labels)
fig_name <- paste(fig_dir, "exchange_convergence.png", sep="")
p <- ggplot(data=v_df, aes(x=x, y=v, colour=r)) +
  geom_line() +
  theme_bw()
print(p)
ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")

for(rs in 1:random_starts){
  title <- paste("Optimal sites from random init ", rs, " with V(D)=", best_v[rs], sep="")
  p <- plot_sites(sampling_surface, best_site_mat[rs, ], title) 
  fig_name <- paste(fig_dir, "optimal_sites_random_start-", rs, ".png", sep="")
  ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
}

print("Avg V(D)")
print(sum(best_v) / random_starts)

# End cluster
stopCluster(clust)

