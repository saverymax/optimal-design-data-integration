################################################
# Module for processing the initial ebird downloda 
# and saving csv as "nuthatch_filtered_for_occ.csv"
# and "nuthatch_filtered_po_2019.txt"

# We will also have functions to fit/save the PP posterior 
# to the ebird data and create the dataframe of the covariates 
# for the optimal design
################################################

library(tidyverse)
library(ggplot2)
library(viridis)
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
parser <- add_option(parser, "--ebird_data_dir", type="character", default="ebird", help="Name of directory containing processed ebird data, within basedir")
parser <- add_option(parser, "--map_file", type="character", help="Nmae of file containing processed US geopackage fil")
parser <- add_option(parser, "--landcover_file", type="character", help="Name of landcover tif")
parser <- add_option(parser, "--modis_file", type="character", help="Name of modis EVI tif")
parser <- add_option(parser, "--elevation_file", type="character", help="Name of elevation tif")
parser <- add_option(parser, "--exp_name", type="character", default="oe_run", help="Name of current experiment, which is used for dir to save output")
parser <- add_option(parser, "--data_reps", type="integer", default=4, help="Number of dataset reps for criterion estimation")
parser <- add_option(parser, "--m", type="integer", default=5, help="Number of sites to survey")
parser <- add_option(parser, "--max_visits", type="integer", default=5, help="Maximum number of time to visit each site")
parser <- add_option(parser, "--min_visits", type="integer", default=1, help="Minimum number of time to visit each site")
parser <- add_option(parser, "--vary_visits", action="store_true", default=F, help="Allow varying survey effort between sites")
parser <- add_option(parser, "--mcmc_iter", type="integer", default=1000, help="Number of MCMC iterations in Stan")
parser <- add_option(parser, "--pp_diagnostic", action="store_true", default=F, help="Boolean for printing diagnostics for point process model")

exp_args <- parse_args(parser)
print(exp_args)