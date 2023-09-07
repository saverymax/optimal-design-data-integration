######################################################## 
# Load the ebird data in the way it will be used in the sampling procedure, and fit it to a Bayesian model.
library(sf)
library(terra)
library(tidyterra)
library(sf)
library(auk)
library(tidyverse)
library(lubridate)
library(exactextractr)


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

# TODO: 
# Based on grid, generate covariate summaries: generate_covariates
# Need to use crs for grid with units meters but then transform to the crs of raster, so we also us the pp data in crs in meters
# TODO: Consider transforming raster (nn for categorical;bilinear for continuous)
subgrid <- create_pp_grid(state_bound, state_pa, crs(crop_lc_rast))
stopifnot(st_crs(subgrid)==st_crs(crop_lc_rast))

p <- ggplot() + 
  geom_spatvector(data=subgrid, color=alpha("#B6B6B6", 0.5)) +
  #geom_spatvector(data=subgrid, aes(fill=counts), color=alpha("#B6B6B6", 0.5)) +
  geom_sf(data=state_po, color=alpha("#FFC81C",0.5), size=0.5)+
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("2.5km by 2.5km grid over Tennessee")
print(p)
fig_name="data/ebird/state_grid.png"
ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")

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
evi_buff
plot(evi_buff)
agg_covars <- generate_gridded_covariates(subgrid, crop_lc_rast, crop_evi_rast, crop_elev_rast)
# TODO: Add frac
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

# Everything looks good which means we can proceed to modelling





