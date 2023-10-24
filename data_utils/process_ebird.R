################################################
# Module for processing the initial ebird download
# and saving csv as "nuthatch_filtered_for_occ.csv"
# and "nuthatch_filtered_po_2019.txt"

# We will also have functions to fit/save the PP posterior 
# to the ebird data and create the dataframe of the covariates 
# for the optimal design

# Example usage:
#"C:/PROGRA~1/R/R-42~1.2/bin/x64/Rscript.exe" data_utils/process_ebird.R 
#  --working_dir=. 
#  --base_data_dir="C:/Users/msavery/OneDrive - UGent/Documents/ghent_phd_spatial_doe/data/" 
#  --ebird_data_dir="ebd_US_bnhnut_201901_201912_smp_relJul-2023" 
#  --map_file="us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg" 
#  --landcover_file="copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif" 
#  --modis_file="modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif" 
#  --elevation_file="elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif" 
#  --save_dir="data/ebird" --pp_fit --pp_diagnostic
################################################

library(tidyverse)
library(ggplot2)
library(viridis)
library(cmdstanr) 
library(bayesplot)
library(reshape2)
library(spatstat)
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
parser <- add_option(parser, "--ebird_data_dir", type="character", default="ebird", help="Name of directory containing processed ebird data csv's, within basedir")
parser <- add_option(parser, "--save_dir", type="character", default="data", help="Directory to save data within basedir")
parser <- add_option(parser, "--map_file", type="character", help="Nmae of file containing processed US geopackage fil")
parser <- add_option(parser, "--landcover_file", type="character", help="Name of landcover tif")
parser <- add_option(parser, "--modis_file", type="character", help="Name of modis EVI tif")
parser <- add_option(parser, "--elevation_file", type="character", help="Name of elevation tif")
parser <- add_option(parser, "--data_reps", type="integer", default=4, help="Number of dataset reps for criterion estimation")
parser <- add_option(parser, "--cell_size", type="integer", default=10000, help="Size of one side of cell in point process grid")
parser <- add_option(parser, "--auk_process", action="store_true", default=F, help="Boolean for running the initial auk filtering steps to generate smallers csv's")
parser <- add_option(parser, "--pp_fit", action="store_true", default=T, help="Boolean to fit the Point Process posterior after data saving steps")
parser <- add_option(parser, "--pp_diagnostic", action="store_true", default=T, help="Boolean for printing diagnostics for point process model")

exp_args <- parse_args(parser)
print(exp_args)

# Source modules
source(file.path(exp_args$working_dir, "data_utils", "load_ebird_data.R"))
source(file.path(exp_args$working_dir, "experimental_design_functions.R"))
source(file.path(exp_args$working_dir, "presence_only_functions.R"))
stan_models_path <- file.path(exp_args$working_dir, "stan_models", "stan_nuthatch_models.R")
source(stan_models_path)
# Set data paths
exp_name <- exp_args$exp_name
save_dir <- file.path(exp_args$working_dir, exp_args$save_dir)
base_data_dir <- file.path(exp_args$base_data_dir)
ebd_download_dir <- file.path(exp_args$ebird_data_dir)
map_path <- exp_args$map_file
modis_path <- exp_args$modis_file
lc_path <- exp_args$landcover_file
elev_path <- exp_args$elevation_file
# Create dir for saving figs and data
dir.create(save_dir)
# Set up the rest of the parameters
# Most of the parameters from the simulation code we don't need. Some we keep, such as m and the prob of detection
data_reps <- exp_args$data_reps
cell_size <- exp_args$cell_size

# This will save the PO data and PA data to the specified dirs in hardcoded file names.
if (exp_args$auk_process==TRUE){
  initial_auk_processing(base_data_dir, ebd_download_dir)
}

rast_surface_path <- file.path(save_dir, "rast_surface.tif")
# Hardcode this for now
map_prj <- st_crs("ESRI:102003")
main_data_handling_oe(map_prj, map_path, base_data_dir, ebd_download_dir, save_dir, lc_path, modis_path, elev_path, rast_surface_path, cell_size)
# Then read in the data and fit the PP model if specified
data_pack <- read_rds(file.path(save_dir, "data_pack.RDS"))
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
print(paste("Total sites: ", sites, sep=""))
# Generate po data to match format of simulation script. 
Y_po <- matrix(rep(data_pack$counts, data_reps), nrow=data_reps, ncol=sites, byrow=T)
r_po_data <- list(Y=Y_po)

stan_path <- file.path(exp_args$working_dir, "stan_models")
# Cell size is in meters but let's work with our parameters in kilometer scale
# We can also set area to 1 for convenience
area_a <- (cell_size/1000)^2
#area_a <- 1
if (exp_args$pp_fit == TRUE){
  fit_point_process_ebird(stan_path, save_dir, sites, area_a, intensity_covars, bias_covars, r_po_data, k_param_intn, k_param_bias, pp_diagnostic,
                          subgrid, rast_surface)
  # Check the posterior
  # The function should save the matrix of MCMC draws, so we can later sample from them or take the means.
  pp_posterior <- read_rds(file.path(save_dir, "pp_posterior_ebird.RDS"))
  print(colMeans(pp_posterior))
  pp_posterior <- matrix(rep(colMeans(pp_posterior), data_reps), nrow=data_reps, byrow = T)
  print(pp_posterior)
}
