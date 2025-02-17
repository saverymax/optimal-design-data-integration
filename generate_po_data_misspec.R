# This script is intended to run additional experiments generating PO data which will lead to misspecification in the working model.
# But first the appropriate set of PO datasets must be generated using a data generating model with an additional covariate.
library(optparse)
library(dplyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)
library(cmdstanr)

parser <- OptionParser()
parser <- add_option(parser, "--working_dir", type="character", default=".", help="Path to the directory containing script")
parser <- add_option(parser, "--exp_name", type="character", default="po_gen_misspec", help="Base name to save data")
parser <- add_option(parser, "--alpha", type="double", default=-2, help="Intercept for intensity")
parser <- add_option(parser, "--beta", type="double", default=0.5, help="Slope for intensity")
parser <- add_option(parser, "--gamma", type="double", default=1, help="Intercept for bias")
parser <- add_option(parser, "--delta", type="double", default=0.25, help="Slope for bias")
parser <- add_option(parser, "--epsilon", type="double", default=1.5, help="Slope for covariate inducing misspec")
parser <- add_option(parser, "--intensity_func", type="character", default="donut", help="Intensity function for sampling surface")
parser <- add_option(parser, "--sd", type="double", default=5, help="Standard deviation for donut intensity surface")
parser <- add_option(parser, "--bias_func", type="character", default="exponential", help="Bias function for sampling surface")
parser <- add_option(parser, "--area", type="integer", default=100, help="Area of region D")
parser <- add_option(parser, "--k", type="integer", default=20, help="Number of sites along one side of grid")
parser <- add_option(parser, "--data_reps", type="integer", default=4, help="Number of dataset reps for criterion estimation")
parser <- add_option(parser, "--gamma_reps", type="integer", default=10, help="Number of dataset reps for criterion estimation")
exp_args <- parse_args(parser)
print(exp_args)

#exp_args <- parse_args(parser, args = c("--gamma_reps=100", "--alpha=-2", "--beta=0.5", "--gamma=1", "--delta=0.5", 
#"--epsilon=1.5", "--intensity_func=donut", "--sd=5", "--bias_func=exponential", "--area=100", "--k=20", "--data_reps=4", "--working_dir=.")) 
#print(exp_args)

source(file.path(exp_args$working_dir, "experimental_design_functions.R"))
source(file.path(exp_args$working_dir, "presence_only_functions.R"))
set.seed(13)

data_dir <- file.path(exp_args$working_dir, "data", "sim_data", "misspec")
dir.create(data_dir)

poisson_process <- '
    data{
      int<lower = 1> N;
      int<lower=0> k_i; // number of predictors for intensity
      int<lower=0> k_b; // number of predictors for bias
      matrix[N, k_i] X; // predictor matrix for intensity
      matrix[N, k_b] Z; // predictor matrix for bias
      array[N] int y;
    }
    parameters{
      real alpha;
      real gamma;
      vector[k_i] beta;       // vector of params for intensity
      vector[k_b] delta;      // vector of params for bias
    }
    model{
      //priors
      target += normal_lpdf(alpha | 0,10);
      target += normal_lpdf(beta | 0,10);
      target += normal_lpdf(gamma | 0,10);
      target += normal_lpdf(delta | 0,10);

      // likelihood
      target += poisson_log_lpmf(y | alpha + X * beta + gamma + Z * delta);
    }
    generated quantities{
      array[N] int y_rep;
      y_rep = poisson_log_rng(alpha + X * beta + gamma + Z * delta);
    }
'

# Will need this model to do the gamma integration
poisson_process_gamma_constant <- '
    data{
      int<lower = 1> N;
      real gamma;
      int<lower=0> k_i; // number of predictors for intensity
      int<lower=0> k_b; // number of predictors for bias
      matrix[N, k_i] X; // predictor matrix for intensity
      matrix[N, k_b] Z; // predictor matrix for bias
      array[N] int y;
    }
    parameters{
      real alpha;
      vector[k_i] beta;       // vector of params for intensity
      vector[k_b] delta;      // vector of params for bias
    }
    model{
      //priors
      target += normal_lpdf(alpha | 0,10);
      target += normal_lpdf(beta | 0,10);
      target += normal_lpdf(delta | 0,10);

      // likelihood
      target += poisson_log_lpmf(y | alpha + X * beta + gamma + Z * delta);
    }
    generated quantities{
      array[N] int y_rep;
      y_rep = poisson_log_rng(alpha + X * beta + gamma + Z * delta);
    }
'

model_path <- "poisson_process_sim_misspec.stan"
model_string <- poisson_process_gamma_no_int
write(model_string, model_path)
model_pp_gamma_no_int <- cmdstan_model(model_path) 
model_path <- "poisson_process_sim_gamma_misspec.stan"
model_string <- poisson_process_gamma_constant
write(model_string, model_path)
model_pp_gamma_int <- cmdstan_model(model_path) 

# Data generating functions for the misspecification that are distinct from the orginal
# PO data generation functions
get_bias_surface_misspecified <- function(sampling_grid, centroid, strength){
  x <- seq(1:k)
  y <- seq(1:k)
  # Generate all possible coordinate points
  grid_points <- expand.grid(x, y)
  # Then compute distance from centroid to every other location
  r <- apply(grid_points, 1, function(x, center_coord){sqrt((center_coord[1] - x[1])^2 + (center_coord[2] - x[2])^2)}, center_coord=centroid)
  # Then create z covariate
  x <- exp(-strength*((r)/5))
  x_1 <- qnorm(0.98*x + .01)
  sampling_grid$aux_e <- x_1
  return(sampling_grid)
}

generate_ppp_data_misspecified <- function(surface_data, params, n_sites, data_reps, corr_matrix, gp_bool, area_D){
  # Generate some random covariate data per site on a grid that will be used to model lambda and b
  # There will be r (data_reps) datasets replicated, for use with the exchange 
  # algorithm/monte carlo integration and to maintain efficiency in the parallel 
  # function computation of V
  lambda <- (area_D / n_sites) * exp(params$alpha + params$beta*surface_data$aux_x)
  eta <- exp(params$gamma + params$delta*surface_data$aux_z + params$epsilon*surface_data$aux_e)
  b <- eta / (1 + eta)
  if (gp_bool==T){
    gp <- exp(rnorm(n_sites, 0, corr_matrix))
    Y_counts <- rpois(n_sites, lambda * b * gp)
  }
  else{
    Y_counts <- rpois(n_sites, lambda * b)
  }
  # Generate coordinate points within each site for each observation.
  Y_coords <- generate_coords(surface_data, Y_counts)
  # Create matrix of data reps, though for po data we just copy everything,
  # as opposed to random generation for each rep, for later convenience in optimal design algorithm
  Y_po <- matrix(rep(Y_counts, data_reps), nrow=data_reps, ncol=n_sites, byrow=T)
  lambdas <- matrix(rep(lambda, data_reps), nrow=data_reps, ncol=n_sites, byrow=T)
  biases <- matrix(rep(b, data_reps), nrow=data_reps, ncol=n_sites, byrow=T)
  return(list(Y=Y_po, lambda=lambdas, bias=biases, Y_coords=Y_coords))
}

# Area for whole space, which allows us to set area for sites based on number of sites.
area_D <- exp_args$area
# For generating data according to GP
gp_bool <- F
data_reps <- exp_args$data_reps
# k is one side of grid
k <- exp_args$k
sites <- k^2
gamma_reps <- exp_args$gamma_reps
centroid <- c(10,4)
centroid_2 <- c(18,18)
n_chains <- 2
mcmc_iter <- 1500

# alpha=2 indicates "good" quality of auxiliary information
alpha <- exp_args$alpha
beta <- exp_args$beta
gamma <- exp_args$gamma
delta <- exp_args$delta
# Fixed value for 3rd covar param. Only strength of covariate itself is changed.
epsilon <- exp_args$epsilon
k_i <- 1
k_b <- 1
po_param_vec <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta, epsilon=epsilon)
# Gamma sample for integration
gamma_sample <- runif(gamma_reps, -1, 1)

# standard deviation of surface
deviation <- exp_args$sd
# Different strengths of the extra data generating covariate
misspec_strength <- c(0, 0.01, 0.05, 0.1)

sampling_surface <- get_sampling_surface_donut(k, deviation)
sampling_surface <- get_bias_surface_exponential(sampling_surface, centroid)

for (si in 1:length(misspec_strength)){
  print(paste("Misspecification effect: ", misspec_strength[si]))
  strength <- misspec_strength[si]
  # here the covariate is forced to be 0, as a strength of 0 doesn't equal 0 via the function 
  if (strength == 0){
    sampling_surface$aux_e = rep(0, nrow(sampling_surface))
  }else{
    sampling_surface <- get_bias_surface_misspecified(sampling_surface, centroid_2, strength)
  }
  # Generate PO data for that particular setting
  r_po_data <- generate_ppp_data_misspecified(sampling_surface, po_param_vec, sites, data_reps, corr_matrix=NA, gp_bool=FALSE, area_D)
  # Save params for later reference
  params <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta)
  
  # Save the PO data sets
  param_setting <- paste(exp_args$exp_name, "_", "ints-", exp_args$intensity_func, "_peak=", 
  deviation ,"_a=", alpha, "_b=", beta, "_g=", gamma, "_d=", delta, "_e=", strength, sep="")
  
  p <- ggplot(sampling_surface, aes(x, y, fill=aux_e)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, name="E") +
    ggtitle(paste("Initial sampling surface, E covariate with misspec scale of ", strength, sep="")) +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  fig_name <- file.path(data_dir, paste("misspecification_strength_", param_setting, ".png", sep=""))
  save_basic_plots(fig_name, p)

  # Bias
  p <- ggplot(sampling_surface, aes(x, y, fill=r_po_data$lambda[data_reps,])) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, name="Lambda") +
    ggtitle("Generated intensity per site") +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  fig_name <- file.path(data_dir, paste("intensity-site_", param_setting, ".png", sep=""))
  save_basic_plots(fig_name, p)
  
  p <- ggplot(sampling_surface, aes(x, y, fill=r_po_data$bias[data_reps,])) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE, "Bias") +
    ggtitle("Generated bias per site") +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  fig_name <- file.path(data_dir, paste("bias-site_", param_setting,".png", sep=""))
  save_basic_plots(fig_name, p)
  
  thinned_intensity <- r_po_data$lambda[data_reps,]*r_po_data$bias[data_reps,]
  p <- ggplot() +
    geom_tile(sampling_surface, mapping=aes(x, y, fill=thinned_intensity, width=1, height=1), alpha=.6) + 
    scale_fill_viridis(discrete=FALSE, name="L*b") +
    ggtitle(paste("Generated PP thinning per site, \nincluding observed PO individuals, misspecification of ", strength, sep="")) +
    geom_point(data=r_po_data$Y_coords, mapping=aes(x=x, y=y), size=2, col="orange") +
    theme(panel.grid.minor = element_line(colour="white")) +
    scale_y_continuous(breaks = seq(0, 20, 1)) +
    scale_x_continuous(breaks = seq(0, 20, 1)) +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed()
  fig_name <- file.path(data_dir, paste("po-thinning-per-site_", param_setting, ".png", sep=""))
  save_basic_plots(fig_name, p)
  
  saveRDS(r_po_data, file=file.path(data_dir, paste(param_setting, ".Rds", sep="")))
  
  po_params <- c("alpha", "beta[1]", "delta[1]")
  # 3 cols for alpha, beta, and delta (gamma is integrated out)
  pp_rep <- matrix(nrow=gamma_reps, ncol=length(po_params))
  # First save the posterior with gamma integration
  for (i in 1:gamma_reps){
    data_site_occ = list(gamma=gamma_sample[i], N=sites, X=as.matrix(sampling_surface$aux_x), y=r_po_data$Y[1,], 
                         Z=as.matrix(sampling_surface$aux_z), k_i=k_i, k_b=k_b)
    fit_pp_gamma_const <- quiet(model_pp_gamma_int$sample(data=data_site_occ, seed=13, chains=n_chains, iter_sampling=mcmc_iter, iter_warmup=500, show_messages=F, refresh=0))
    fit_means <- fit_pp_gamma_const$summary(po_params)$mean
    pp_rep[i,] <- fit_means
  }
  # Get mean of each parameter from the iterations
  pp_posterior <- matrix(rep(colMeans(pp_rep[,1:2]), data_reps), nrow=data_reps, byrow = T)
  saveRDS(pp_posterior, file.path(data_dir, paste("pp_posterior_", param_setting, ".Rds", sep="")))

  # And then save it without gamma integration
  po_params <- c("alpha", "beta[1]", "gamma", "delta[1]")
  data_site_occ = list(gamma=gamma_sample[i], N=sites, X=as.matrix(sampling_surface$aux_x), y=r_po_data$Y[1,], 
                       Z=as.matrix(sampling_surface$aux_z), k_i=k_i, k_b=k_b)
  fit_pp_gamma_no_int <- quiet(model_pp_gamma_no_int$sample(data=data_site_occ, seed=13, chains=n_chains, iter_sampling=mcmc_iter, iter_warmup=500, show_messages=F, refresh=0))
  fit_means <- fit_pp_gamma_no_int$summary(po_params)$mean
  pp_posterior <- matrix(rep(fit_means, data_reps), nrow=data_reps, byrow = T)
  saveRDS(pp_posterior, file.path(data_dir, paste("pp_posterior_", param_setting, ".Rds", sep="")))
}
