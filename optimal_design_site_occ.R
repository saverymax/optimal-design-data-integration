###########################
# Sandbox script for optimal design in .qmd file
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

fig_dir <- "figures\\optimal_design\\"
source("experimental_design_functions.R")
source("presence_only_functions.R")
source("stan_models\\stan_site_occupancy_models.R")

# R=1000 datasets for monte carlo approx
exp_args <- list(model_selection=3, data_reps=10, random_starts=1, p_logging=F)
data_reps <- exp_args$data_reps
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
# 36/4 was used in paper
m <- 10
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
r_survey_data <- generate_data(data_reps, sampling_surface, corr_matrix, p_0, b_0, b_1, sigma, n_surveys, sites)

# Then generate R presence-only datasets
alpha <- -2 
beta <- 2
gamma <- -4
delta <- 0.5
params <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta)
# Provide occupancy maps and number of data reps, as well as params and sampling surface, to generate pp data.
#r_po_data <- generate_ppp_data_r(sampling_surface, params, sites, r_survey_data$occupancy, data_reps)
r_po_data <- generate_ppp_data_r(sampling_surface, params, sites, data_reps)
Y_positive_indices <- which(r_po_data$Y>0)
# This is the data at which there are counts > 0
r_po_data$Y[Y_positive_indices]
# The first replicate dataset
r_survey_data$occupancy[1,]
r_survey_data$Y[1,]
# So that gives us R complete datasets for 400 sites.
# These will stay fixed throughout the rest of the procedure

# For the coordinate exchange algorithm, we will need to precompute the nearest neighbors.
# I take a naive approach here of choosing the top l neighbors
l <- 8
nearest_neighbors <- get_neighbors(sampling_surface, l)
dim(nearest_neighbors)
nearest_neighbors

# Use the final generation iteration
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
c(1:random_starts)
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
      estimate_mat <- matrix(nrow=data_reps, ncol=1)
      # TODO: Parallelize this.
      for (r in 1:data_reps){
        # For each rth dataset get the m randomly chosen sites
        selected_data <- r_survey_data$Y[r, select_idx]
        selected_occ <- r_survey_data$occupancy[r, select_idx]
        # Here I use the same X covariate for the presence-only data as used for the survey data. The difference is that
        # the survey data here is a subset (select_sites) of the sites, whereas the presence only data 
        # needs covariates for the whole grid to approximate the expected count in the entire region.
        if (model_selection==3){
          data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, X=select_sites$aux_x, Y=selected_data, PO=r_po_data$Y[r,], 
                             X_po=sampling_surface$aux_x, Z_po=sampling_surface$aux_z, model_diag=0)
        }
        else if(model_selection==1){
          data_site_occ = list(n_surveys=n_surveys, n_sites=m, total_sites=sites, X=select_sites$aux_x, X_all=sampling_surface$aux_x, Y=selected_data)
        }
        else{
          stop("No model selected")
        }
        # refresh=0 turns off messages except errors from stan
        # quiet function silences stan output 
        fit <- quiet(model$sample(data=data_site_occ, seed=13, chains=1, 
                                  iter_sampling=1000, iter_warmup=100, refresh=0, show_messages=F))
        #fit <- model$sample(data=data_site_occ, seed=13, chains=1, 
        #                          iter_sampling=1000, iter_warmup=100)
        if (p_logging==T){
          print("Logging posterior")
          posterior <- fit$draws()
          print(fit$summary(variables=params))
          
          color_scheme_set("mix-blue-pink")
          p_trace <- mcmc_trace(posterior,  pars = params,
                                facet_args = list(nrow = 2, labeller = label_parsed))
          print(p_trace + facet_text(size = 15))
  
          plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
          p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = c("alpha", "beta")) + plot_title
          print(p_post)
          
          if(model_selection==3){
            plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
            p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = c("gamma", "delta")) + plot_title
            print(p_post)
          }
          
          plot_title <- ggtitle(paste("Posterior distribution of detection probability, mean and 90% interval"))
          p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", pars = c("p")) + plot_title
          print(p_post)
          
          # Check issues with divergences
          color_scheme_set("darkgray")
          nuts_fit <- nuts_params(fit)
          diverge_p <- mcmc_parcoord(posterior, pars = params, np = nuts_fit, alpha=.1, 
                                     np_style=parcoord_style_np(div_alpha=1, div_size=.5))
          print(diverge_p)
          
          ppd_count_df <- data.frame(x=sampling_surface$x, y=sampling_surface$y, 
                                     occ_prob=fit$summary(variables=generated_vars[1])$mean)
          
          p <- ggplot(ppd_count_df, aes(x, y, fill=occ_prob)) + 
            geom_tile() +
            scale_fill_viridis(discrete=FALSE) +
            ggtitle("Generated occupancy probability per site")
          print(p)
          
        }
        # Average estimate after R iterations through the datasets.
        # Need to select_idx the sites since it generates for all of them.
        # The 2 index is the ppd for Z
        gen_occupancy <- fit$summary(variables=generated_vars[2])$mean[select_idx]
        # Take the mean of the generated quantity for the posterior estimate
        # TODO: Check if the magnitude of V makes sense
        estimate_mat[r,1] <- design_criteria(criteria="brier", sim_occ=gen_occupancy, obs_occ=selected_occ)
      }
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
      fig_name <- paste(fig_dir, "site_locs_iter_", exchange_iter, ".png", sep="")
      ggsave(fig_name, plot=p, dpi = 300)
      # Select new data using the best coordinates (switch out local points)
      # TODO: Compare to deterministic exchange
      select_idx <- exchange_coordinates(best_select_idx, current_site, s, nearest_neighbors, l)
      select_sites <- sampling_surface[select_idx,] 
      # Then go to the next site
    }
    # If after a complete iteration through all the sites the best sites haven't changed
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
ggplot(data=v_df, aes(x=x, y=v, colour=r)) +
  geom_line()

for(rs in 1:random_starts){
  title <- paste("Optimal sites from random init ", rs, sep="")
  plot_sites(sampling_surface, best_site_mat[rs, ], title) 
}
