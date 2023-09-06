######################################################## 
# Load the ebird data in the way it will be used in the sampling procedure, and fit it to a Bayesian model.
library(sf)
library(terra)
library(tidyterra)
library(sf)
library(auk)
library(tidyverse)
library(lubridate)


select <- dplyr::select
source(file.path("data_utils", "load_ebird_data.R"))

# Set important global variables 
base_data_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\"
ebd_download_dir <- "ebd_US_bnhnut_201901_201912_smp_relJul-2023"
lc_path <- "copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif"
modis_path <- "modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif"
elev_path <- "elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif"
map_prj <- st_crs("ESRI:102003")

map_path <- file.path(base_data_dir, "us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg")
state_bound <- load_map(map_prj, map_path)

# TODO: Include preprocessing for creating the file "nuthatch_filtered_for_occ.csv"
nuthatch_list <- load_ebird(base_data_dir, ebd_download_dir, state_bound, map_prj)
state_pa <- nuthatch_list$state_pa
state_po <- nuthatch_list$state_po

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
stopifnot(crs(crop_elev_rast)==crs(prj_state))
p <- ggplot() + 
  geom_spatraster(data=crop_elev_rast) +
  #geom_sf(data=state_pp_prj, color=alpha("darkgrey", 0.3))+
  geom_sf(data=state_po, color=alpha("#FFC81C",0.5), size=0.5)+
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.4) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch and elevation in Tennessee")
print(p)

# TODO: 
# Create grid and counts: create_pp_grid
# Based on grid, generate covariate summaries: generate_covariates
