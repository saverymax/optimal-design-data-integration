library(ggplot2)
library(viridis)
library(hrbrthemes)

#user_dir <- "C:\\Users\\saver\\Documents\\ghent_phd\\code\\point-processes\\"
user_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\code\\point_processes\\"
source(paste(user_dir, "mcmc_functions.R", sep=""))
source(paste(user_dir, "experimental_design_functions.R", sep=""))


plot_sites <- function(sampling_surface, select_idx, title){
  p <- ggplot(sampling_surface, aes(x, y, fill=aux)) + 
        geom_tile() +
        geom_point(data=sampling_surface[select_idx,], aes(x=x, y=y), colour = "white", size = 3) +
        scale_fill_viridis(discrete=FALSE) +
        ggtitle(title)
  print(p)
}
 
# R=1000 datasets for monte carlo approx
data_reps <- 10
# k is one side of grid
k <- 20
sites <- k^2
p_0 <- 0.7
b_0 <- 0
# b_01=2 indicates "good" quality of auxiliary information
b_1 <- 2
# This assumes spatial variance of 1, which was used in the paper (see supplement)
sigma <- 1
# number of surveys at site is equal to n or 0.
n_surveys <- 5
# There will be m sites selected for sampling
# 36 was used in paper
m <- 10
# Based on the size of grid get the auxiliary data and coordinates
# Because these don't really depend on any random variables and just location 
# in the grid, there will be one fixed dataset throughout the optimization
sampling_surface <- get_sampling_surface(k)
sampling_surface
# Next, we use this data to generate the rest of the datasets
# The data generating function will sample R occupancy maps|params
# and then R complete datasets|occupancy maps
corr_matrix <- specify_corr(sampling_surface[,1:2])
r_survey_data <- generate_data(data_reps, sampling_surface, corr_matrix, p_0, b_0, b_1, sigma, n_surveys, sites)
# The first replicate dataset
r_survey_data$occupancy[1,]
r_survey_data$Y[1,]
# So that gives us R complete datasets for 400 sites.
# These will stay fixed throughout the rest of the procedure

# We then begin to prepare for the exchange algorithm
# we first randomly select m sites
# from the sampling surface, which includes the auxiliary information
select_idx <- sample(1:sites, m, replace=F)
select_idx
select_sites <- sampling_surface[select_idx,]
select_sites
    
# For the coordinate exchange algorithm, we will need to precompute the nearest neighbors.
# I take a naive approach here of choosing the top l neighbors
# 9 is nice because it includes exactly the adjacent cells
l <- 15
nearest_neighbors <- get_neighbors(sampling_surface, l)
nearest_neighbors

# Choose number of iterations for metropolis-hastings
mh_iter <- 1000
# Replicate n surveys for data to be passed to mh, for m sites 
total_rep <- rep(n_surveys, m)
total_exchanges <- 100
v_mat <- matrix(nrow=total_exchanges, ncol=1)
p_logging <- T
# Start timer 
start_time <- Sys.time()
for (ex in 1:total_exchanges){
  print(paste("Exchange iteration: ", ex))
  # Data structure for each score estimate
  estimate_mat <- matrix(nrow=data_reps, ncol=1)
  # Then compute the posterior based on those sites and generated data, for each r dataset
  for (r in 1:data_reps){
    # For each rth dataset get the m randomly chosen sites
    selected_data <- r_survey_data$Y[r, select_idx]
    selected_data
    selected_occ <- r_survey_data$occupancy[r, select_idx]
    selected_occ

    # Format into df in the way I initially wrote mh function.
    # Center the aux data first
    observed_df <- data.frame(Y=selected_data, N=total_rep, aux=select_sites$aux-mean(select_sites$aux))
    posterior <- mh(observed_df, mh_iter, corr_matrix, logging=F)
    if (p_logging==T){
      print("acceptance rate:")
      print(sum(posterior$acceptance) / nrow(posterior$acceptance))
      # Posterior estimate for betas
      print("parameter estimates, beta0, beta 1, p")
      print(mean(posterior$sim_data[,1]))
      print(mean(posterior$sim_data[,2]))
      print(mean(posterior$sim_data[,3]))
      plot(type="l",1:mh_iter, posterior$sim_data[,1], main="Traceplot of beta_0")
      plot(type="l",1:mh_iter, posterior$sim_data[,2], main="Traceplot of beta_1")
      plot(type="l",1:mh_iter, posterior$sim_data[,3], main="Traceplot of detection probability p")
    }
    # print("criterion")
    #print(design_criteria(criteria="brier", posterior$occupancy, selected_occ))
    # Average estimate after R iterations through the datasets.
    estimate_mat[r,1] <- design_criteria(criteria="brier", posterior$occupancy, selected_occ)
  }
  # Once the posterior is computed on each of R datasets, find the average score:
  new_v_est <- sum(estimate_mat) / data_reps
  v_mat[ex, 1] <- new_v_est
  print("Design score from most recent exchange")
  print(new_v_est)
  
  if (ex==1){
    current_v_est <- new_v_est
    best_select_idx <- select_idx
    print("Initial row ids")
    print(select_idx)
    print("Initial coordinates")
    print(sampling_surface[best_select_idx,1:2])
    print("Initial design score")
    print(new_v_est)
    title <- paste("Inital spatial design: v=", round(new_v_est, 5), sep="")
    plot_sites(sampling_surface, best_select_idx, title)
  } else if (new_v_est < current_v_est){
    # Set new best indices
    best_select_idx <- select_idx
    print("New optimal row ids")
    print(best_select_idx)
    print("New coordinates")
    print(sampling_surface[best_select_idx,1:2])
    current_v_est <- new_v_est
    print(current_v_est)
    title <- paste("Optimal spatial design: v=", round(new_v_est, 5), sep="")
    plot_sites(sampling_surface, best_select_idx, title)
  }
  else{
    title <- paste("Non-optimal spatial design: v=", round(new_v_est, 5), sep="")
    plot_sites(sampling_surface, select_idx, title)
  }
  # Select new data (switch out local points)
  select_idx <- exchange_coordinates(best_select_idx, nearest_neighbors, l)
  select_sites <- sampling_surface[best_select_idx,] 
  selected_occ <- r_survey_data$occupancy[best_select_idx]
  selected_data <- r_survey_data$Y[best_select_idx]
  # Then iterate again.
}
end_time <- Sys.time()
run_time <- end_time - start_time
run_time


