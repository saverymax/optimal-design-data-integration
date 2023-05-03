##########################
# Sandbox for running parallel stan models
# TODO: Remove complexity and fit simple stan model.
##########################

library(ggplot2)
library(viridis)
library(hrbrthemes)
library(cmdstanr) 
library(bayesplot)
library(dplyr)
library(tidyr)
library(reshape2)
library(parallel)

source("experimental_design_functions.R")
source("presence_only_functions.R")
source("stan_models\\stan_site_occupancy_models.R")

n_cores <- detectCores()
# R=1000 datasets for monte carlo approx
exp_args <- list(model_selection=3, m=10, data_reps=5, random_starts=5, p_logging=F)
fig_dir <- paste("figures\\optimal_design_model-", exp_args$model_selection,  "_m-", exp_args$m, "_r-", exp_args$data_reps, "\\", sep="")
dir.create(fig_dir)
data_reps <- exp_args$data_reps
# k is one side of grid
k <- 20
sites <- k^2
alpha <- -2 
beta <- 2
gamma <- -4
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
sampling_surface <- get_sampling_surface(k)
# Filter for only a quarter of the grid.
# If we filter, we need to change the total sites as well
sites <- sites/4
stopifnot(sites>m)
sampling_surface <- sampling_surface %>% dplyr::filter(x<11, y<11)
sampling_surface
# Next, we use this data to generate the rest of the datasets
# The data generating function will sample R occupancy maps|params
# and then R complete datasets|occupancy maps
corr_matrix <- specify_corr(sampling_surface[,1:2])
link_func <- "cloglog"
r_survey_data <- generate_data(data_reps, sampling_surface, corr_matrix, p_0, alpha, beta, sigma, n_surveys, sites, link=link_func)
r_survey_data$occupancy[1,]
r_survey_data$Y[1,]
r_survey_data$theta[1,]

# Then generate R presence-only datasets
params <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta)
# Provide occupancy maps and number of data reps, as well as params and sampling surface, to generate pp data.
#r_po_data <- generate_ppp_data_r(sampling_surface, params, sites, r_survey_data$occupancy, data_reps)
r_po_data <- generate_ppp_data_r(sampling_surface, params, sites, data_reps)
Y_positive_indices <- which(r_po_data$Y>0)
# This is the data at which there are counts > 0
r_po_data$Y[Y_positive_indices]
# So that gives us R complete datasets for 400 sites.
# These will stay fixed throughout the rest of the procedure

# For the coordinate exchange algorithm, we will need to precompute the nearest neighbors.
# I take a naive approach here of choosing the top l neighbors
l <- 8
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
model_strings <- c(cloglog_site_occupancy, site_occupany_detection, poisson_process_site_occupancy)
model_selection <- exp_args$model_selection
write(model_strings[model_selection], model_path)
model <- cmdstan_model(model_path) 
# These are the parameters to report, though this will be model dependent
if (model_selection==3){
  params <- c('p', 'alpha', 'beta', 'gamma', 'delta')
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
# Start timer
start_time <- Sys.time()
for (r_start in 1:random_starts){
  print(paste("random initialization ", r_start, sep=""))
  # we first randomly select m sites, for each run of the exchange algorithm,
  # from the sampling surface, which includes the auxiliary information.
  # These sites will have n_i = n, the others will have n_i = 0
  # Only sites with n_i=n will contribute to likelihood for the site-occupancy model.
  select_idx <- sample(1:sites, m, replace=F)
  select_idx
  # Initialize for exchange algorithm
  best_select_idx <- select_idx
  select_sites <- sampling_surface[select_idx,]
  select_sites
  # Convergence condition will be met where full iteration through all sampled sites results in no change
  # in sites. ie no changes in sites can improve criterion
  convergence_cond <- FALSE
  exchange_iter <- 0
  v_vec <- c()
  while (convergence_cond==FALSE){
    exchange_iter <- exchange_iter + 1
    print(paste("New exchange iteration: ", exchange_iter))
    # Data structure for each score estimate
    # Then compute the posterior based on those sites and generated data, for each r dataset
    # We iterate through the sites, computing the estimate of $V(D)$ for each so that we explore the effect of each site on the design
    # It is also possible to iterate through all neighbors of the site as well
    for (s in 1:length(select_idx)){
      print(paste("ex iter: ", exchange_iter, ", current site index: ", s, sep=""))
      current_site <- select_idx[s]
      print(paste("current site: ", current_site, sep=""))
      clust <- makeCluster(n_cores)
      # Generate estimate matrix, with V estimate for each data rep
      combined_df <- cbind(r_survey_data$occupancy, r_survey_data$Y, r_po_data$Y)
      dim(combined_df)
      parApply(clust, r_survey_data$occupancy, 1, FUN=mean)
      # TODO: Getting stan to work here?
      estimate_mat <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel, model, n_surveys, data_reps, m, sites, sampling_surface, 
                       select_idx, select_sites, params, generated_vars, model_selection)
      # Nice functional version of this, but I need to put in form for apply?
      estimate_mat <- estimate_v(model, n_surveys, data_reps, m, sites, sampling_surface, 
                       select_idx, select_sites, r_survey_data, r_po_data, 
                       p_logging, params, generated_vars, model_selection)
      # Once the posterior is computed on each of R datasets, find the average score:
      new_v_est <- sum(estimate_mat) / data_reps
      print("Design score from most recent exchange")
      print(new_v_est)
      v_vec <- c(v_vec, new_v_est)
      
      # In the first exchange iteration and for the first site, save the first V(D) as the best
      if ((exchange_iter==1)&(s==1)){
        current_v_est <- new_v_est
        print("Initial row ids")
        print(select_idx)
        print("Initial coordinates")
        print(sampling_surface[best_select_idx,1:2])
        print("Initial design score")
        print(new_v_est)
        title <- paste("Inital spatial design: v=", round(new_v_est, 10), sep="")
        p <- plot_sites(sampling_surface, best_select_idx, title)
      } 
      else if (new_v_est < current_v_est){
        current_v_est <- new_v_est
        title <- paste("New optimal spatial design: v=", round(new_v_est, 10), sep="")
        # Plot the new best site compared to the previous selection, but need to reverse arguments to function
        # Since the best one is now old and the select will is now the best
        p <- plot_sites_vs_best(sampling_surface, best_select_idx, select_idx, title)
        # Then set new best indices
        best_select_idx <- select_idx
        print("New optimal row ids")
        print(best_select_idx)
        print("New coordinates")
        print(sampling_surface[best_select_idx,1:2])
        print("New optimal design score")
        print(current_v_est)
      }
      else{
        title <- paste("Non-optimal spatial design: v=", round(new_v_est, 10), 
                       ", vs current optimal design: v=", round(current_v_est, 10), sep="")
        print("No change in optimal design")
        p <- plot_sites_vs_best(sampling_surface, select_idx, best_select_idx, title)
      }
      fig_name <- paste(fig_dir, "site_locs_rand-start-", r_start, "_ex-iter_", exchange_iter, "_site-iter-", s, ".png", sep="")
      ggsave(fig_name, plot=p, dpi=300)
      # Select new data using the best coordinates (switch out local points)
      # TODO: Compare to deterministic exchange
      select_idx <- exchange_coordinates(best_select_idx, current_site, s, nearest_neighbors, l)
      select_sites <- sampling_surface[select_idx,] 
      # Then go to the next site
    }
    # If after a complete iteration through all the sites, the best sites haven't changed
    # then we can call that convergence
    # If it's the first iteration we need to initialize the best sites
    if(exchange_iter==1){
      best_iter_idx <- best_select_idx
    }
    else if(all(best_iter_idx==best_select_idx)){
      convergence_cond <- TRUE
    }
    else{
      print("Previous best sites and current best sites")
      print(best_iter_idx)
      print(best_select_idx)
      best_iter_idx <- best_select_idx
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
ggsave(fig_name, plot=p, dpi = 300)

for(rs in 1:random_starts){
  title <- paste("Optimal sites from random init ", rs, " with V(D)=", best_v[rs], sep="")
  p <- plot_sites(sampling_surface, best_site_mat[rs, ], title) 
  fig_name <- paste(fig_dir, "optimal_sites_random_start-", rs, ".png", sep="")
  ggsave(fig_name, plot=p, dpi = 300)
}

print("Avg V(D)")
sum(best_v) / random_starts

