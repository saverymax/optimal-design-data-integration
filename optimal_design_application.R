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
library(auk)
library(lubridate)
library(exactextractr)

set.seed(13)

# R=1000 datasets for monte carlo approx
# Create command line arguments
parser <- OptionParser()
parser <- add_option(parser, "--working_dir", type="character", default=".", help="Path to the directory containing code to source for the main script")
parser <- add_option(parser, "--base_data_dir", type="character", default="./data", help="Path to the directory containing covariate data")
parser <- add_option(parser, "--base_data_dir", type="character", default="./data", help="Path to the directory containing covariate data")
parser <- add_option(parser, "--exp_name", type="character", default="oe_run", help="Name of current experiment, which is used for dir to save output")
parser <- add_option(parser, "--data_reps", type="integer", default=4, help="Number of dataset reps for criterion estimation")
parser <- add_option(parser, "--m", type="integer", default=5, help="Number of sites to survey")
parser <- add_option(parser, "--max_visits", type="integer", default=5, help="Maximum number of time to visit each site")
parser <- add_option(parser, "--min_visits", type="integer", default=1, help="Minimum number of time to visit each site")
parser <- add_option(parser, "--vary_visits", action="store_true", default=F, help="Allow varying survey effort between sites")
parser <- add_option(parser, "--model_selection", type="integer", default=1, help="Occupancy model to use. For application, currently using only 1 model")
parser <- add_option(parser, "--random_starts", type="integer", default=3, help="Number of random starts to run the exchange")
parser <- add_option(parser, "--exch_iter", type="integer", default=3, 
                     help="Number of iterations of exchange before ending optimization. Recommended is 20 but default is set low for test runs.")
parser <- add_option(parser, "--mcmc_iter", type="integer", default=1000, help="Number of MCMC iterations in Stan")
parser <- add_option(parser, "--p_logging", action="store_true", default=F, help="Boolean for logging information about posterior estimates")
parser <- add_option(parser, "--pp_diagnostic", action="store_true", default=F, help="Boolean for printing diagnostics for point process model")
parser <- add_option(parser, "--v_parallel", action="store_true", default=F, help="Boolean for parallel computation of V criterion")
parser <- add_option(parser, "--cores", type="integer", default=4, help="Number of cores to use for parallel processing")
parser <- add_option(parser, "--p", type="double", default=0.7, help="Probability of detection")
parser <- add_option(parser, "--area", type="integer", default=100, help="Area of region D")

exp_args <- parse_args(parser)
print(exp_args)
stopifnot(exp_args$p_logging==F)
select <- dplyr::select
sort <- base::sort

# Set important global variables 
base_data_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\"
ebd_download_dir <- "ebd_US_bnhnut_201901_201912_smp_relJul-2023"
lc_path <- "copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
modis_path <- "modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif"
elev_path <- "elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif"
map_prj <- st_crs("ESRI:102003")

# Source modules
source(file.path(exp_args$working_dir, "data_utils", "load_ebird_data.R"))
source(file.path(exp_args$working_dir, "experimental_design_functions.R"))
source(file.path(exp_args$working_dir, "presence_only_functions.R"))
stan_models_path <- file.path(exp_args$working_dir, "stan_models", "stan_nuthatch_models.R")
source(stan_models_path)
map_path <- file.path(base_data_dir, "us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg")

exp_name <- exp_args$exp_name
exp_dir <- file.path(exp_args$working_dir, "experimental_runs", exp_name)
fig_dir <- file.path(exp_dir, "figures")
stan_dir <- file.path(exp_dir, "stan")
dir.create(exp_dir)
dir.create(fig_dir)
dir.create(stan_dir)
data_reps <- exp_args$data_reps
# Most of the parameters from the simulation code we don't need. Some we keep, such as m and the prob of detection
p_0 <- exp_args$p
# There will be m sites selected for sampling
m <- exp_args$m
# That set the intial experimental environment up. Now we can focus on loading our covariates and data
# This involves a lot of data processing
state_bound <- load_map(map_prj, map_path)

# TODO: Include preprocessing file for creating the file "nuthatch_filtered_for_occ.csv"
nuthatch_list <- load_ebird(base_data_dir, ebd_download_dir, state_bound, map_prj)
state_pa <- nuthatch_list$state_pa
state_po <- nuthatch_list$state_po

# TODO: Save these babies
# Plot the data
p <- ggplot() + 
  #geom_spatvector(data=state_grid, fill = 'transparent', colour="lightblue") +
  geom_sf(data = state_bound, color=alpha("white",0.9)) + 
  geom_sf(data=state_pa, color=alpha("darkgrey",0.2), size=0.5)+
  geom_sf(data=state_po, color=alpha("#FFC81C",0.5), size=0.5)+
  ggtitle("Brown-headed Nuthatch checklists in Tennessee") +
  theme_minimal() +
  theme(text=element_text(size=10)) +
  coord_sf()
print(p)

crop_lc_rast <- load_landcover(base_data_dir, lc_path, state_bound)
# Convert crs to the landcover
state_pa_prj <- terra::project(vect(state_pa), crop_lc_rast)
state_po_prj <- terra::project(vect(state_po), crop_lc_rast)
prj_state <- terra::project(vect(state_bound), crop_lc_rast)
state_pa_prj <- st_transform(state_pa, crs(crop_lc_rast))
state_po_prj <- st_transform(state_po, crs(crop_lc_rast))
prj_state <- st_transform(state_bound, crs(crop_lc_rast))

# Plot landcover
p <- ggplot() + 
  geom_spatraster(data=crop_lc_rast) +
  #geom_sf(data=state_pp_prj, color=alpha("darkgrey", 0.3))+
  geom_sf(data=state_po_prj, color=alpha("#FFC81C",0.5), size=0.5)+
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.4) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch and landcover in Tennessee")
print(p)
#fig_name="data/ebird/landcover.png"
#ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")

# Get modis data
crop_evi_rast <- load_evi(base_data_dir, modis_path, prj_state, crop_lc_rast)
p <- ggplot() + 
  geom_spatraster(data=crop_evi_rast) +
  #geom_sf(data=state_pp_prj, color=alpha("darkgrey", 0.3))+
  geom_sf(data=state_po_prj, color=alpha("#FFC81C",0.5), size=0.5)+
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.4) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch and EVI in Tennessee")
print(p)

crop_elev_rast <- load_elevation(base_data_dir, elev_path, prj_state)
# Shouldn't need to project either source of rast for elevation rasters
stopifnot(crs(crop_elev_rast)==crs(crop_lc_rast))
stopifnot(st_crs(crop_elev_rast)==st_crs(prj_state))
p <- ggplot() + 
  geom_spatraster(data=crop_elev_rast) +
  #geom_sf(data=state_pp_prj, color=alpha("darkgrey", 0.3))+
  geom_sf(data=state_po, color=alpha("#FFC81C",0.5), size=0.5)+
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.4) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch and elevation in Tennessee")
print(p)

subgrid <- create_pp_grid(state_bound, crs(crop_lc_rast))
stopifnot(st_crs(subgrid)==st_crs(crop_lc_rast))

p <- ggplot() + 
  geom_spatvector(data=subgrid, color=alpha("#B6B6B6", 0.5)) +
  #geom_spatvector(data=subgrid, aes(fill=counts), color=alpha("#B6B6B6", 0.5)) +
  geom_sf(data=state_po, color=alpha("#FFC81C",0.5), size=0.5)+
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("2.5km by 2.5km grid over Tennessee")
print(p)
#fig_name="data/ebird/state_grid.png"
#ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")

pp_counts <- aggregate_point_counts(subgrid, state_po_prj)
pp_counts[which(pp_counts>1)]
hist(pp_counts[which(pp_counts>1)])
subgrid$counts <- pp_counts

# We can plot the subgrid (either as sfc or sf) with spatvector
# make_grid produces sfc but st_sf turns this to sf, to which 
# we can add attributes such as point counts
p <- ggplot() + 
  geom_spatvector(data=subgrid, aes(fill=counts), color=alpha("#B6B6B6", 0.5)) +
  scale_fill_viridis_c(begin=0.2, end=1, option="magma",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch aggregated intensity in Tennessee")
print(p)

# Next calculate covariate aggregation based on the grid.
evi_buff <- exact_extract(crop_evi_rast, subgrid, "mean", progress=T)
plot(evi_buff)
# TODO: Add ebird covariates from the buuffer: number_observers and duration_minutes
agg_covars <- generate_gridded_covariates(subgrid, crop_lc_rast, crop_evi_rast, crop_elev_rast, state_po_prj)
agg_covars$lc_frac
subgrid$evi <- agg_covars$evi_agg
subgrid$elev <- agg_covars$elev_agg
subgrid_covars <- bind_cols(subgrid, agg_covars$lc_frac)

p <- ggplot() + 
  geom_spatvector(data=subgrid_covars, aes(fill=counts), color=alpha("#B6B6B6", 0.5)) +
  scale_fill_viridis_c(begin=0.2, end=1, option="magma",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch aggregated intensity in Tennessee")
print(p)

# rasterized version of site evi 
evi_rast <- terra::rasterize(subgrid_covars, crop_evi_rast, field="evi")
#evi_rast[is.na(evi_rast)] <- 0
plot(evi_rast)

# rasterized version of site elev 
elev_site_rast <- terra::rasterize(subgrid_covars, crop_evi_rast, field="elev")
elev_site_rast[is.na(elev_site_rast)] <- 0
plot(elev_site_rast)

p <- ggplot() + 
  #  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatraster(data=evi_rast) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, na.value="white") +
  theme_minimal()+
  ggtitle("Nuthatch EVI in Tennessee at Nuthatch sites")
print(p)

p <- ggplot() + 
  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatraster(data=elev_site_rast) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, na.value="white") +
  theme_minimal()+
  ggtitle("Elevation in Tennessee at Nuthatch sites")
print(p)

colnames(subgrid_covars)
covar_names <- c(colnames(subgrid_covars)[5:(length(colnames(subgrid_covars))-1)])
covar_names
length(covar_names)
# List of all covar labels
covar_labels <- c("shrubs", "herbaceous_vegetation", "agriculture", "urban", "bare/sparse", "water", "wetland", "closed_forest_evergreen_needle", 
                  "close_forest_decid_broad", "closed_forest_mixed", "closed_forest_unknown", "open_forest_evergreen_needle", 
                  "open_forest_decid_broad leaf", "open_forest_unknown")
covar_labels
names(covar_names) <- covar_labels
covar_names

# Take a look at all the fractional covariates by plotting rasterized versions of them, which is faster
for(i in 1:length(covar_names[1:length(covar_names)])){
  covar <- covar_names[i]
  covar_name <- covar_labels[i]
  covar_rast <- terra::rasterize(vect(subgrid_covars), crop_evi_rast, field=covar)
  p <- ggplot() + 
    geom_spatraster(data=covar_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle(paste(covar, "; ", covar_name, " in Tennessee at Nuthatch sites", sep=""))
  print(p)
  #fig_name <- paste("data/ebird/", covar, "_covar_tenn.png", sep="")
  #ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
}

covar_df <- st_drop_geometry(subgrid_covars)
colnames(covar_df)
#intensity_covars <- covar_df[c("evi", "elev", "frac_111", )]
intensity_covars <- covar_df[c("evi", "elev", "frac_111", "frac_115", "frac_124")]
bias_covars <- covar_df[c("frac_50")]#"duration_minutes", "number_observers",  "effort_distance_km", "time_observations_started", 
k_param_intn <- length(colnames(intensity_covars)) 
k_param_bias <- length(colnames(bias_covars)) 
dim(as.matrix(covar_df))
# Number of sites depends upon how the grid is created
sites <- nrow(covar_df)
sites
sum(covar_df$counts)
length(covar_df$counts)
# Generate po data to match format of simulation script. 
# TODO: Put this all in function or preprocessing script
Y_po <- matrix(rep(covar_df$counts, data_reps), nrow=data_reps, ncol=sites, byrow=T)
r_po_data <- list(Y=Y_po)

# Fit the point process model to the PO data
model_path <- "nuthatch_poisson_process.stan"
model_string <- nuthatch_poisson_process
write(model_string, model_path)
data_site_occ = list(N=sites, X=intensity_covars, y=r_po_data$Y[1,], Z=bias_covars, k_i=k_param_intn, k_b=k_param_bias)
model <- cmdstan_model(model_path) 
fit <- model$sample(data=data_site_occ, seed=13, chains=3, iter_sampling=2000, iter_warmup=500) 

print(fit$summary())
all_params <- c("alpha", "gamma", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "delta[1]")
#all_params <- c("alpha", "gamma", "beta[1]", "beta[2]", "beta[3]", "delta[1]")
#intensity_params <- c("alpha", "beta[1]", "beta[2]", "beta[3]")
intensity_params <- c("alpha", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]")
params_intercept <- c("alpha", "gamma")

if (exp_args$pp_diagnostic == T){
  posterior <- fit$draws(all_params)
  color_scheme_set("mix-blue-pink")
  p_trace <- mcmc_trace(posterior,
                        facet_args = list(nrow = 2, labeller = label_parsed))
  print(p_trace + facet_text(size = 15))
  
  plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
  mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[1]")) + plot_title
  mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[2]")) + plot_title
  mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[3]")) + plot_title
  mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[4]")) + plot_title
  mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[5]")) + plot_title
  
  plot_title <- ggtitle(paste("Posterior distributions of detection probability, mean and 90% interval"))
  p_post <- mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars = params_intercept) + plot_title
  print(p_post)
  
  plot_title <- ggtitle(paste("Posterior distributions of detection probability, mean and 90% interval"))
  p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", pars = c("delta[1]")) + plot_title
  print(p_post)
  
  mcmc_intervals(fit$draws(), pars=all_params)
  mcmc_hist(fit$draws(), pars = all_params)
  mcmc_pairs(fit$draws(), pars=all_params)
  mcmc_scatter(fit$draws(), pars=c('alpha', 'gamma'))
  mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[2]'))
  mcmc_scatter(fit$draws(), pars=c('beta[2]', 'beta[3]'))
  mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[3]'))
  mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[4]'))
  mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[5]'))
}
  
# Examine the ppd for y
generated_yrep <- fit$draws("y_rep", format="matrix")
yrep_means <- colMeans(generated_yrep)
# Create new sf for the preds
subgrid_geo <- st_geometry(subgrid_covars) %>% st_sf()
subgrid_geo$preds <- yrep_means
# Then plot the predictions
pp_pred_rast <- terra::rasterize(vect(subgrid_geo), crop_evi_rast, field="preds")
p <- ggplot() + 
  geom_spatraster(data=pp_pred_rast) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
  #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
  theme_minimal()+
  ggtitle("Predicted counts at nuthatch sites")
print(p)

# Take expectation from the PP params
pp_posterior <- fit$draws(intensity_params, format="matrix")
pp_posterior <- matrix(rep(colMeans(pp_posterior), data_reps), nrow=data_reps, byrow = T)
print(dim(pp_posterior))

# We fit the PP model to the PO data and use the expectation of this posterior to generate the PA data
# There is a question of whether to use draws from the posterior or just the expectation. This will have to be resolved
# later when discussing Bayesian optimal  design
# It turns out that if we want to compare different survey efforts it is convenient to have pre-generated datasets for
# each number of visits
# The number of surveys can differ between sites so it is a vector
# We will select m sites to have these visits, otherwise the sites will have 0 visits
# We will compare the number of visits during the optimization
link_func <- "cloglog"
if (exp_args$vary_visits == TRUE){
  visits <- c(exp_args$min_visits, exp_args$max_visits)
  print("Creating datasets for varying survey effort between sites")
  r_survey_data_n1  <- generate_ebird_pa(data_reps, intensity_covars, p_0, pp_posterior, visits[1], sites, link=link_func)
  r_survey_data_n5  <- generate_ebird_pa(data_reps, intensity_covars, p_0, pp_posterior, visits[2], sites, link=link_func)
} else{
  visits <- c(exp_args$max_visits)
  print("Creating datasets for fixed survey effort across sites")
  r_survey_data_n1  <- generate_ebird_pa(data_reps, intensity_covars, p_0, pp_posterior, visits[1], sites, link=link_func)
  r_survey_data_n5  <- generate_ebird_pa(data_reps, intensity_covars, p_0, pp_posterior, visits[1], sites, link=link_func)
}
print(paste("Current visit options:", paste(visits, collapse=" ")))

# For the coordinate exchange algorithm, we will need to precompute the nearest neighbors.
# Need matrix of size: matrix(nrow=nrow(surface), ncol=l). In the simulate we computed distances between 
# cells. st_distance allows us to do the same here, between grid cells
l <- 4
nearest_neighbors <- get_neighbors_sf_grid(subgrid, l)
dim(nearest_neighbors)

# Create plots of the occupancy and probability maps
d_examine <- ifelse(data_reps<10, data_reps, 10)
for(d_i in 1:d_examine){
  occ_map <- st_geometry(subgrid_covars) %>% st_sf()
  occ_map$occ <- r_survey_data_n1$occupancy[d_i,]
  occ_map$prob <- r_survey_data_n1$theta[d_i,]
  occ_map$counts <- r_survey_data_n1$Y[d_i,]
  
  # Then plot the predictions
  occ_rast <- terra::rasterize(vect(occ_map), crop_evi_rast, field="occ")
  p <- ggplot() + 
    geom_spatraster(data=occ_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle("Occupancy maps generated over Prior Predictive Distribtion")
  print(p)
  
  prob_rast <- terra::rasterize(vect(occ_map), crop_evi_rast, field="prob")
  p <- ggplot() + 
    geom_spatraster(data=prob_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle("Probability maps generated over Prior Predictive Distribtion")
  print(p)
  
  count_rast <- terra::rasterize(vect(occ_map), crop_evi_rast, field="counts")
  p <- ggplot() + 
    geom_spatraster(data=count_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle("Counts at sites generated over Prior Predictive Distribtion")
  print(p)
}

# We can experiment with these models in the exchange algorithm
# Model 1 uses just basic priors over params. Model 3 uses PO data as prior. Model 2
# uses different link function but is otherwise the same as model 1. Model has constant
# intensity, which doesn't make that much sense to use in this case.
model_strings <- list(
	"nuthatch_site_occ"=nuthatch_poisson_process_site_occupancy
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
params <- c("alpha", "gamma", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "delta[1]")
print(paste("Using params", paste(params, collapse=" ")))
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
    estimate_vec <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel_nuthatch, model, possible_visits, m, sites, 
                             intensity_covars, bias_covars, site_idx, select_sites, generated_vars, k_param_intn, k_param_bias, 
                             model_selection, exp_args$mcmc_iter)
    # Compatible format with non-parallel v
    estimate_mat <- matrix(estimate_vec, nrow=data_reps, ncol=1)
  }
  else{
    # TODO: fix arguments, create function
    estimate_mat <- estimate_v_nuthatch(model, possible_visits, data_reps, m, sites, sampling_surface, 
                               site_idx, select_sites, r_survey_data_n5, r_po_data,
                               p_logging, params, generated_vars, model_selection, exp_args$mcmc_iter)
  }
  # Once the posterior is computed on each of R datasets, find the average score:
  new_v_est <- sum(estimate_mat) / data_reps
  print("Design score from initial exchange")
  print(new_v_est)
  v_vec <- c(v_vec, new_v_est)
  current_v_est <- new_v_est
  print("Initial row ids")
  print(site_idx)
  print("Initial design score")
  print(new_v_est)
  title <- paste("Inital spatial design: v=", round(new_v_est, 10), sep="")
  p <- plot_sites(sampling_surface, site_idx, title)
  fig_name <- file.path(fig_dir, paste("initial_design_", r_start, ".png", sep=""))
  ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
  
  while ((convergence_cond==FALSE) & (exchange_iter<=exp_args$exch_iter)){
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
      # Set the current sites to the best from iteration over sites neighbors
      current_visits <- possible_visits
      # We first iterate through the neighbors of each site, and compute V for each exchange. 
      # The initial estimate will be for our initial design.
      for (visit in visits){
        n_count <- 0
        current_visits[s] <- visit
        for(nn in neighbor_set){
          # Handle duplicates in sites because the visit optimization is not compatible with 
          # the same site occurring multiple times in site_idx
          # Check that the neighbor isn't in the main set of sites
          if (nn %in% site_idx){
            print("skipping site")
            print(site_idx)
            print(nn)
            next
          }
          n_count <- n_count + 1
          # Select new data using the neighbors (switch out local points)
          neighbor_idx <- exchange_coordinates_deterministic(best_neighbor_idx, s, nn)
          select_sites <- sampling_surface[neighbor_idx,] 
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
          print("posible visits")
          print(possible_visits)
          print("current vistits")
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
          print("neighbors")
          print(neighbor_idx)
          print("current neighbor")
          print(nn)
          #print("current data selction after visits altered")
          #print(r_survey_data$Y[,neighbor_idx])
          if (exp_args$v_parallel==T){
            combined_df <- cbind(r_survey_data$occupancy, r_survey_data$Y, r_po_data$Y)
            estimate_vec <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel_nuthatch, model, current_visits, m, sites, 
                                     intensity_covars, bias_covars, neighbor_idx, select_sites, generated_vars, k_param_intn, k_param_bias, 
                                     model_selection, exp_args$mcmc_iter)
            # Compatible format with non-parallel v
            estimate_mat <- matrix(estimate_vec, nrow=data_reps, ncol=1)
          }
          else{
            estimate_mat <- estimate_v(model, current_visits, data_reps, m, sites, sampling_surface, 
                                       neighbor_idx, select_sites, r_survey_data, r_po_data,
                                       p_logging, params, generated_vars, model_selection, exp_args$mcmc_iter)
          }
          # Once the posterior is computed on each of R datasets, find the average score:
          new_v_est <- sum(estimate_mat) / data_reps
          v_vec <- c(v_vec, new_v_est)
          
          if (new_v_est < current_v_est){
             current_v_est <- new_v_est
             title <- paste("New optimal spatial design: v=", round(new_v_est, 10), sep="")
             # Plot the new best site compared to the previous selection, but need to reverse arguments to function
             # Use current_visits as optimal visits.
             p <- plot_sites_vs_best_nuthatch(intensity_covars, current_site, best_neighbor_idx, neighbor_idx, possible_visits, current_visits, title)
             # Then set new best indices
             best_neighbor_idx <- neighbor_idx
             possible_visits <- current_visits
             print("New optimal row ids")
             print(best_neighbor_idx) 
             print("Optimum survey effort")
             print(possible_visits)
             print("New coordinates")
             print(sampling_surface[best_neighbor_idx,1:2])
             print("New optimal design score")
             print(current_v_est)
          }
          else{
            title <- paste("Non-optimal spatial design: v=", round(new_v_est, 10), 
                           "\nvs current optimal design: v=", round(current_v_est, 10), sep="")
            #print("No change in optimal design")
            p <- plot_sites_vs_best_nuthatch(sampling_surface, current_site, neighbor_idx, best_neighbor_idx, current_visits, possible_visits, title)
          }
          fig_name <- file.path(fig_dir, paste("site_locs_rand-start-", r_start, "_ex-iter_", 
                            exchange_iter, "_site-iter-", s, "_effort_", visit, "_nn-iter", n_count,".png", sep=""))
          ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
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
fig_name <- file.path(fig_dir, "exchange_convergence.png")
p <- ggplot(data=v_df, aes(x=x, y=v, colour=r)) +
  geom_line() +
  theme_bw()
print(p)
ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")

# Plot best sites
for(rs in 1:random_starts){
  title <- paste("PO data and Optimal sites from random init ", rs, "\n with V(D)=", best_v[rs], sep="")
  #p <- plot_sites(sampling_surface, best_site_mat[rs, ], title) 
  p <- plot_po_optimal_sites(sampling_surface, r_po_data, best_site_mat[rs,], optimal_visit_mat[rs,], title)
  fig_name <- file.path(fig_dir, paste("optimal_sites_random_start-", rs, ".png", sep=""))
  ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
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

