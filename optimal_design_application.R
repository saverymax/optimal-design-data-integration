######################################################
# Script for optimal design of surveys applied to Brown-headed nuthatch 
######################################################

library(tidyverse)
library(ggplot2)
library(viridis)
#library(hrbrthemes)
library(cmdstanr) 
library(bayesplot)
library(reshape2)
library(spatstat)
library(parallel)
library(optparse)
library(openxlsx)
# Spatial libs
library(sf)
library(terra)
library(tidyterra)

set.seed(13)

# R=1000 datasets for monte carlo approx
# Create command line arguments
parser <- OptionParser()
parser <- add_option(parser, "--working_dir", type="character", default=".", help="Path to the directory containing code to source for the main script")
parser <- add_option(parser, "--base_data_dir", type="character", default="./data", help="Path to the directory containing covariate data")
parser <- add_option(parser, "--ebird_data_dir", type="character", default="ebird", help="Name of directory containing processed ebird data, within basedir")
parser <- add_option(parser, "--data_save_dir", type="character", default="data/ebird", help="Directory to in which preporcessed covariate data is saved")
parser <- add_option(parser, "--map_file", type="character", help="Name of file containing processed US geopackage fil")
parser <- add_option(parser, "--landcover_file", type="character", help="Name of landcover tif")
parser <- add_option(parser, "--modis_file", type="character", help="Name of modis EVI tif")
parser <- add_option(parser, "--elevation_file", type="character", help="Name of elevation tif")
parser <- add_option(parser, "--exp_name", type="character", default="oe_application", help="Name of current experiment, which is used for dir to save output")
parser <- add_option(parser, "--m", type="integer", default=5, help="Number of sites to survey")
parser <- add_option(parser, "--max_visits", type="integer", default=5, help="Maximum number of time to visit each site")
parser <- add_option(parser, "--min_visits", type="integer", default=1, help="Minimum number of time to visit each site")
parser <- add_option(parser, "--vary_visits", action="store_true", default=F, help="Allow varying survey effort between sites")
parser <- add_option(parser, "--model_selection", type="integer", default=1, help="Occupancy model to use: (1) with or (2) without PO data")
parser <- add_option(parser, "--data_reps", type="integer", default=1, help="number of dataset reps for criterion estimation")
parser <- add_option(parser, "--random_starts", type="integer", default=1, help="Number of random starts to run the exchange")
parser <- add_option(parser, "--exch_iter", type="integer", default=2, 
                     help="Number of iterations of exchange before ending optimization. Recommended is 20 but default is set low for test runs.")
parser <- add_option(parser, "--mcmc_iter", type="integer", default=1000, help="Number of MCMC iterations in Stan")
parser <- add_option(parser, "--p_logging", action="store_true", default=F, help="Boolean for logging information about posterior estimates")
parser <- add_option(parser, "--v_parallel", action="store_true", default=F, help="Boolean for parallel computation of V criterion")
parser <- add_option(parser, "--cores", type="integer", default=4, help="Number of cores to use for parallel processing")
parser <- add_option(parser, "--p", type="double", default=0.2, help="Probability of detection")
parser <- add_option(parser, "--cell_size", type="integer", default=10000, help="Size of one side of cell in point process grid")

exp_args <- parse_args(parser)
print(exp_args)
stopifnot(exp_args$p_logging==F)
select <- dplyr::select
sort <- base::sort

# Set important global variables if we're just running within Rstudio. hacky :)
#base_data_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\"
#ebd_download_dir <- "ebd_US_bnhnut_201901_201912_smp_relJul-2023"
#lc_path <- "copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
#modis_path <- "modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif"
#elev_path <- "elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif"
#map_path <- "us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg"
#data_save_dir <- file.path(".", "data/ebird")

# Source modules
source(file.path(exp_args$working_dir, "data_utils", "load_ebird_data.R"))
source(file.path(exp_args$working_dir, "experimental_design_functions.R"))
source(file.path(exp_args$working_dir, "presence_only_functions.R"))
stan_models_path <- file.path(exp_args$working_dir, "stan_models", "stan_nuthatch_models.R")
source(stan_models_path)
# Set data paths
exp_name <- exp_args$exp_name
exp_dir <- file.path(exp_args$working_dir, "experimental_runs", exp_name)
base_data_dir <- file.path(exp_args$base_data_dir)
ebd_download_dir <- file.path(exp_args$ebird_data_dir)
map_path <- exp_args$map_file
modis_path <- exp_args$modis_file
lc_path <- exp_args$landcover_file
elev_path <- exp_args$elevation_file
data_save_dir <- file.path(exp_args$working_dir, exp_args$data_save_dir)
# Create dir for figures and stan files
fig_dir <- file.path(exp_dir, "figures")
stan_dir <- file.path(exp_dir, "stan")
dir.create(exp_dir)
dir.create(fig_dir)
dir.create(stan_dir)
# Set up the rest of the parameters
# Most of the parameters from the simulation code we don't need. Some we keep, such as m and the prob of detection
data_reps <- exp_args$data_reps
p_0 <- exp_args$p
# There will be m sites selected for sampling
m <- exp_args$m
cell_size <- exp_args$cell_size

# That set the intial experimental environment up. Now we can focus on loading our covariates and data
# This involves a lot of data processing which we pack intothe main_data_handling function
# It is necessary to run the process_ebird script first to generate the rds file. See documentation about data generation.
rast_surface_path <- file.path(data_save_dir, "rast_surface.tif")
data_pack <- read_rds(file.path(data_save_dir, "data_pack.RDS"))
# Unpack
bias_covars <- data_pack$bias_covars
intensity_covars <- data_pack$intensity_covars
# Standardize covars
standardize <- function(x){ 
  z <- (x - mean(x)) / sd(x) 
  return( z)
}
intensity_covars <- apply(intensity_covars, 2, standardize)
k_param_intn <- data_pack$k_param_intn
k_param_bias <- data_pack$k_param_bias
site_centroids <- data_pack$site_centroids
subgrid <- data_pack$subgrid
state_po_prj <- data_pack$state_po_prj
rast_surface <- rast(rast_surface_path) 
sites <- nrow(intensity_covars)
# Generate po data to match format of simulation script. 
Y_po <- matrix(rep(data_pack$counts, data_reps), nrow=data_reps, ncol=sites, byrow=T)
r_po_data <- list(Y=Y_po)

stan_path <- file.path(exp_args$working_dir, "stan_models")
# Cell size is in meters but let's work with our parameters in kilometer scale
area_a <- (cell_size/1000)^2
# Set to one for practical purposes
#area_a <- 1
# The point process model needs to be fit in the process_ebird script.Please see documentation regarding that before using the optimal design code 
if(!file.exists(file.path(data_save_dir, "pp_posterior_ebird.RDS"))){
  stop("Please run process_ebird.R with the relevant CLI arguments before running the optimal design! See the documentation for more details")
}else{
  pp_posterior <- read_rds(file.path(data_save_dir, "pp_posterior_ebird.RDS"))
  pp_posterior <- matrix(rep(colMeans(pp_posterior), data_reps), nrow=data_reps, byrow = T)
  print("Intensity fit from point process")
  print(pp_posterior)
}
# There is a question of whether to use draws from the posterior or just the expectation. This will have to be resolved
# later when discussing Bayesian optimal  design

# Then generate PA data.
# It turns out that if we want to compare different survey efforts it is convenient to have pre-generated datasets for
# each number of visits, though it is coded a bit awkwardly here.
# The number of surveys can differ between sites so it is a vector
# We will select m sites to have these visits, otherwise the sites will have 0 visits
# We will compare the number of visits during the optimization
link_func <- "cloglog"
if (exp_args$vary_visits == TRUE){
  visits <- c(exp_args$min_visits, exp_args$max_visits)
  print("Creating datasets for varying survey effort between sites")
  r_survey_data_n1 <- generate_ebird_pa(data_reps, area_a, intensity_covars, p_0, pp_posterior, visits[1], sites, link=link_func)
  r_survey_data_n5 <- generate_ebird_pa(data_reps, area_a, intensity_covars, p_0, pp_posterior, visits[2], sites, link=link_func)
} else{
  visits <- c(exp_args$max_visits)
  print("Creating datasets for fixed survey effort across sites")
  r_survey_data_n1 <- generate_ebird_pa(data_reps, area_a, intensity_covars, p_0, pp_posterior, visits[1], sites, link=link_func)
  r_survey_data_n5 <- generate_ebird_pa(data_reps, area_a, intensity_covars, p_0, pp_posterior, visits[1], sites, link=link_func)
}
print(paste("Current visit options:", paste(visits, collapse=" ")))

# For the coordinate exchange algorithm, we will need to precompute the nearest neighbors.
# Need matrix of size: matrix(nrow=nrow(surface), ncol=l). In the simulate we computed distances between 
# cells. st_distance allows us to do the same here, between grid cells
l <- 4
nearest_neighbors <- get_neighbors_sf_grid(site_centroids, l)
print("Dim of NNs")
print(dim(nearest_neighbors))

# Have to use subgrid_covars for this since it still has geometry and not intensity_/bias_covars
# Create plots of the occupancy and probability maps
# Plot the data for n=5 so that when vary_visits=T we can observe counts above 1.
d_examine <- ifelse(data_reps<3, data_reps, 3)
for(d_i in 1:d_examine){
  occ_map <- st_geometry(subgrid) %>% st_sf()
  occ_map$occ <- r_survey_data_n5$occupancy[d_i,]
  occ_map$prob <- r_survey_data_n5$theta[d_i,]
  occ_map$counts <- r_survey_data_n5$Y[d_i,]
  
  # Then plot the predictions
  occ_rast <- terra::rasterize(vect(occ_map), rast_surface, field="occ")
  p <- ggplot() + 
    geom_spatraster(data=occ_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle("Occupancy maps generated over Prior Predictive Distribtion")
  fig_name=file.path(exp_dir, paste("gen_pa_occ", d_i, ".png", sep=""))
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  prob_rast <- terra::rasterize(vect(occ_map), rast_surface, field="prob")
  p <- ggplot() + 
    geom_spatraster(data=prob_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle("Probability maps generated over Prior Predictive Distribtion")
  fig_name=file.path(exp_dir, paste("gen_pa_prob", d_i, ".png", sep=""))
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  count_rast <- terra::rasterize(vect(occ_map), rast_surface, field="counts")
  p <- ggplot() + 
    geom_spatraster(data=count_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle("Counts at sites generated over Prior Predictive Distribtion")
  fig_name=file.path(exp_dir, paste("gen_pa_y", d_i, ".png", sep=""))
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
}

# We can experiment with these models in the exchange algorithm
# Model 1 uses just basic priors over params. Model 3 uses PO data as prior. Model 2
# uses different link function but is otherwise the same as model 1. Model has constant
# intensity, which doesn't make that much sense to use in this case.
model_strings <- list(
	"nuthatch_site_occ"=nuthatch_poisson_process_site_occupancy,
	"nuthatch_site_occ_no_po"=nuthatch_site_occ_no_po
)
model_selection <- exp_args$model_selection
model_name <- names(model_strings)[model_selection]
print(paste("Using model", model_name, model_selection))
model_path <- file.path(stan_dir, paste(model_name, ".stan", sep=""))
print(paste("Cmdstan model", model_strings[[model_selection]]))
write(model_strings[[model_selection]], model_path)
model <- cmdstan_model(model_path) 
# These are the parameters to report
# Only the PO prior model needs gamma and delta. 
# If I look at models with and without PO data then I will need to add more options
if (model_selection==1){
  params <- c("alpha", "gamma", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "delta[1]")
}else if (model_selection==2){
  params <- c("alpha", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]")
}else{
  stop("No additional models implemented")
}
print(paste("Using params", paste(params, collapse=" "), "in model", model_selection))
generated_vars <- c('g_theta_gen')

# The reich paper repeats the entire exchange algorithm procedure 10 times,
# and retains solution with lowest V(D)
# Currently I am not repeating the procedure. However, it accounts for uncertainty associated with random m sites 
# selected for sampling
p_logging <- exp_args$p_logging
random_starts <- exp_args$random_starts
v_list <- vector(mode="list", length=random_starts)
names(v_list) <- c(1:random_starts)
best_v <- vector(mode="numeric", length=random_starts)
best_site_mat <- matrix(nrow=random_starts, ncol=m)
optimal_visit_mat <- matrix(nrow=random_starts, ncol=m)

# Initiate parallel processing if specified
if (exp_args$v_parallel==T){
  n_cores <- detectCores()
  print(paste("Using parallel processing. Cores detected: ", n_cores, ". Currently using ", exp_args$cores, " cores."))
  clust <- makeCluster(exp_args$cores)
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
  print("Initial row ids")
  print(site_idx)
  # Initialize for exchange algorithm
  # I don't need best_neighbor_idx but it allows me to not modify the 
  # the vector that is looped over during the exchange. Even though this concurrent looping should be ok, as the sites are independently 
  # exchanged, but for organization purposes they are separate variables. 
  best_neighbor_idx <- site_idx
  # Vector to hold potential sampling effort at each site
  # optimal_visits will be initiated in algorithm
  possible_visits <- rep(visits[length(visits)], m)
  # Initial sites
  # Use the covariates for just the intensity since these correspond to the sites in the PA likelihood
  select_sites <- intensity_covars[site_idx,]
  # Convergence condition will be met where full iteration through all sampled sites results in no change
  # in sites. ie no changes in sites can improve criterion
  convergence_cond <- FALSE
  exchange_iter <- 0
  v_vec <- c()
  # Compute v for initial design
  # Not comparing sampling effort here
  if (exp_args$v_parallel==T){
    combined_df <- cbind(r_survey_data_n5$occupancy, r_survey_data_n5$Y, r_po_data$Y)
    estimate_vec <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel_nuthatch, model, possible_visits, m, sites, area_a, 
                             intensity_covars, bias_covars, site_idx, select_sites, generated_vars, k_param_intn, k_param_bias, 
                             model_selection, exp_args$mcmc_iter)
    if (class(estimate_vec)!="numeric"){
      print("Error in parallel computation")
      stop(print(estimate_vec))
    }
    # Compatible format with non-parallel v
    estimate_mat <- matrix(estimate_vec, nrow=data_reps, ncol=1)
  }
  else{
    estimate_mat <- estimate_v_nuthatch(model, possible_visits, data_reps, m, sites, area_a, intensity_covars, bias_covars, 
                               site_idx, select_sites, r_survey_data_n5, r_po_data,
                               p_logging, params, generated_vars, k_param_intn, k_param_bias, model_selection, exp_args$mcmc_iter)
  }
  # Once the posterior is computed on each of R datasets, find the average score:
  new_v_est <- sum(estimate_mat) / data_reps
  print("Design score from initial exchange")
  print(new_v_est)
  v_vec <- c(v_vec, new_v_est)
  current_v_est <- new_v_est
  print("Initial design score")
  print(new_v_est)
  title <- paste("Inital spatial design ", r_start, ": v=", round(new_v_est, 10), sep="")
  # Use the evi as the background for these plots
  p <- plot_sites_ebird(site_centroids, rast_surface, site_idx, title)
  fig_name <- file.path(fig_dir, paste("initial_design_", r_start, ".png", sep=""))
  ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm", bg='white', device="png", type="cairo")
  
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
      print("neighbor set")
      print(neighbor_set)
      # Set the current sites to the best from iteration over sites neighbors
      current_visits <- possible_visits
      # We first iterate through the neighbors of each site, and compute V for each exchange. 
      # The initial estimate will be for our initial design.
      for (visit in visits){
        n_count <- 0
        current_visits[s] <- visit
        for(nn in neighbor_set){
          print("current neighbor")
          print(nn)
          # Handle duplicates in sites (if two nearly adjacent sites have the same neighbor) 
          # because the visit optimization is not compatible with 
          # the same site occurring multiple times in site_idx
          # Check that the neighbor isn't in the main set of sites. This is the only case that should trigger the if statement
          if (nn %in% site_idx){
            print("skipping site")
            print(site_idx)
            print(nn)
            next
          }
          n_count <- n_count + 1
          # Select new data using the neighbors (switch out local points)
          neighbor_idx <- exchange_coordinates_deterministic(best_neighbor_idx, s, nn)
          select_sites <- intensity_covars[neighbor_idx,] 
          # Then set the survey effort for current id set
          # This code chunk will "double index" if the current site s already has been selected 
          # to have 1 visit, but that is ok since we just need the data at nn to correspond to 1 visit.
          # For example if possible_visits == c(1,5,5) and then we are the neighbor of the 1st site 
          # so that current_visits == c(1, 5,5) we will select that index==1 neighbor here 
          # If visits is fixed between sites this vector will always be of 0 length.
          visit_idx <- c(neighbor_idx[which(current_visits==exp_args$min_visits)])
          print("current s and site")
          print(s)
          print(current_site)
          print("current visit")
          print(visit)
          print("possible visits")
          print(possible_visits)
          print("current visits")
          print(current_visits)
          print("ids with visit 1")
          print(visit_idx)
          r_survey_data <- r_survey_data_n5
          r_survey_data$occupancy[, visit_idx] <- r_survey_data_n1$occupancy[, visit_idx]
          r_survey_data$Y[, visit_idx] <- r_survey_data_n1$Y[, visit_idx]
          stopifnot(all(r_survey_data$Y[, visit_idx]<=exp_args$min_visits))
          # Check that we're selecting right sites
          if (visit==exp_args$min_visits){
            stopifnot(nn%in%visit_idx)
          }
          # Check that vector is empty if not varying visits
          if(exp_args$vary_visits==F){
            stopifnot(length(visit_idx)==0)
          }
          # Don't really need theta as it's only for data generation purposes
          r_survey_data$theta[, visit_idx] <- r_survey_data_n1$theta[, visit_idx]
          print("site set")
          print(site_idx)
          print("best neighbor idx")
          print(best_neighbor_idx)
          print("temp neighbor idx")
          print(neighbor_idx)
          #print("current data selction after visits altered")
          #print(r_survey_data$Y[,neighbor_idx])
          if (exp_args$v_parallel==T){
            combined_df <- cbind(r_survey_data$occupancy, r_survey_data$Y, r_po_data$Y)
            estimate_vec <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel_nuthatch, model, current_visits, m, sites, area_a, 
                                     intensity_covars, bias_covars, neighbor_idx, select_sites, generated_vars, k_param_intn, k_param_bias, 
                                     model_selection, exp_args$mcmc_iter)
            # Error handling
            if (class(estimate_vec)!="numeric"){
              print("Error in parallel computation")
              stop(print(estimate_vec))
            }
            # Compatible format with non-parallel v
            estimate_mat <- matrix(estimate_vec, nrow=data_reps, ncol=1)
          }
          else{
            estimate_mat <- estimate_v_nuthatch(model, current_visits, data_reps, m, sites, area_a, intensity_covars, bias_covars, 
                                                neighbor_idx, select_sites, r_survey_data, r_po_data,
                                                p_logging, params, generated_vars, k_param_intn, k_param_bias, model_selection, exp_args$mcmc_iter)
          }
          # Once the posterior is computed on each of R datasets, find the average score:
          new_v_est <- sum(estimate_mat) / data_reps
          v_vec <- c(v_vec, new_v_est)
          
          if (new_v_est < current_v_est){
             current_v_est <- new_v_est
             title <- paste("New optimal spatial design: v=", round(new_v_est, 10), sep="")
             # Plot the new best site compared to the previous selection, but need to reverse arguments to function, so
             # use current_visits as optimal visits.
             print(possible_visits)
             print(current_visits)
             print(best_neighbor_idx)
             print(neighbor_idx)
             print(current_site)
             p <- plot_sites_vs_best_ebird(site_centroids, rast_surface, current_site, best_neighbor_idx, neighbor_idx, possible_visits, current_visits, title)
             # Then set new best indices
             # best_neighbor and possible_visits will hold the optimal for the local round of iteration
             best_neighbor_idx <- neighbor_idx
             possible_visits <- current_visits
             print("New optimal row ids")
             print(best_neighbor_idx) 
             print("Optimum survey effort")
             print(possible_visits)
             print("New optimal design score")
             print(current_v_est)
          }
          else{
            title <- paste("Non-optimal spatial design: v=", round(new_v_est, 10), 
                           "\nvs current optimal design: v=", round(current_v_est, 10), sep="")
            #print("No change in optimal design")
            p <- plot_sites_vs_best_ebird(site_centroids, rast_surface, current_site, neighbor_idx, best_neighbor_idx, current_visits, possible_visits, title)
          }
          fig_name <- file.path(fig_dir, paste("site_locs_rand-start-", r_start, "_ex-iter_", 
                            exchange_iter, "_site-iter-", s, "_effort_", visit, "_nn-iter", n_count,".png", sep=""))
          ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm", bg="white", device="png", type="cairo")
          # Then go to the next neighbor or site
        }
      }
      # Set design to best from iteration through visits AND neighbors for one site
      # If there is no change from any neighbors, site_idx will not change
      site_idx <- best_neighbor_idx
    }
    # If after a complete iteration through all the sites, the best sites haven't changed
    # then we can call that convergence
    # If it's the first iteration we need to initialize the best sites
    if(exchange_iter==1){
      best_iter_idx <- site_idx
      optimal_visits <- possible_visits
    }
    else if(all(best_iter_idx==site_idx)&all(optimal_visits==possible_visits)){
      convergence_cond <- TRUE
    }
    else{
      print("Previous best sites and current best sites")
      print(best_iter_idx)
      print(site_idx)
      best_iter_idx <- site_idx
      print("Previous optimal survey effort and new optimal effort")
      print(optimal_visits)
      print(possible_visits)
      # Need to track if there is a change in visits over the course of full iteration
      optimal_visits <- possible_visits
    }
    # Print run time per exchange
    cur_time <- Sys.time()
    run_time <- cur_time - start_time
    print(paste("Current run time is", run_time))
    
    # Write running results for current start after each exchange
    v_df <- data.frame(x=1:length(v_vec), v=v_vec)
    fig_name <- file.path(fig_dir, paste("running_exchange_convergence_random_start-", r_start, ".png", sep=""))
    p <- ggplot(data=v_df, aes(x=x, y=v)) +
      geom_line() +
      theme_bw()
    ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
    
    # Plot best sites
    title <- paste("PO data and Optimal sites from random init ", r_start, "\n with V(D)=", current_v_est, sep="")
    p <- plot_po_optimal_sites_ebird(site_centroids, rast_surface, state_po_prj, best_iter_idx, optimal_visits, title)
    fig_name <- file.path(fig_dir, paste("running_optimal_sites_random_start-", r_start, ".png", sep=""))
    ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm", bg="white", device="png", type="cairo")
  }
  
  v_list[[r_start]] <- v_vec
  best_site_mat[r_start,] <- best_iter_idx
  optimal_visit_mat[r_start,] <- optimal_visits
  best_v[r_start] <- current_v_est
}
end_time <- Sys.time()
run_time <- end_time - start_time
print("Total run time:")
print(run_time)


# Plot the convergence of V(D)
v_concat <- c()
x_concat <- c()
rep_labels <- c()
# Format v for data frame
for(i in 1:random_starts){
  v_concat <- c(v_concat, v_list[[i]])
  l_vec <- rep(paste("y", i, sep=""), length(v_list[[i]]))
  x_concat <- c(x_concat, 1:length(l_vec))
  rep_labels <- c(rep_labels, l_vec)
}

v_df <- data.frame(x=x_concat, v=v_concat, r=rep_labels)
fig_name <- file.path(fig_dir, "final_exchange_convergence.png")
p <- ggplot(data=v_df, aes(x=x, y=v, colour=r)) +
  geom_line() +
  theme_bw()
print(p)
ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")

# Plot best sites
for(rs in 1:random_starts){
  title <- paste("PO data and Optimal sites from random init ", rs, "\n with V(D)=", best_v[rs], sep="")
  p <- plot_po_optimal_sites_ebird(site_centroids, rast_surface, state_po_prj, best_site_mat[rs,], optimal_visit_mat[rs,], title)
  fig_name <- file.path(fig_dir, paste("optimal_sites_random_start-", rs, ".png", sep=""))
  ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm", bg="white", device="png", type="cairo")
}

print("V list")
print(v_list)
print("Best sites")
print(best_site_mat)
print("Optimal visits")
print(optimal_visit_mat)
print("best V")
print(best_v)
print("Avg V(D)")
print(sum(best_v) / random_starts)

write_results(random_starts, best_v, best_site_mat, optimal_visit_mat, v_df, exp_dir, exp_name)

# End cluster
stopCluster(clust)

