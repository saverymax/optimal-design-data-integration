# Script for using Ebird data for occupancy modelling. The big thing is learning how to process ebird data for occupancy modelling
# From https://strimas.com/ebp-workshop-au/advanced.html#advanced-unmarked: 
# First, we need to extract a subset of observations that are suitable for occupancy modeling. 
#   In particular, occupancy models typically require data from repeated visits to a single site during a time frame over 
#   which the population can be considered closed. The auk function filter_repeat_visits() is designed to extract subsets 
#   of eBird data that meet these criteria.
library(auk)
library(terra)
library(tidyterra)
library(tidyverse)
library(rnaturalearth)
library(sf)
library(lubridate)
library(dggridR)
library(unmarked)
library(exactextractr)
library(cartography)
library(raster)

select <- dplyr::select
map <- purrr::map
projection <- raster::projection


base_data_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\"

#################### 
# Modelling occupancy
# Following https://cornelllabofornithology.github.io/ebird-best-practices/occupancy.html
ebird <- read_csv(file.path(base_data_dir, "eb_tutorial_data/data/ebd_woothr_june_bcr27_zf.csv")) %>% 
  mutate(year = year(observation_date), species_observed = as.integer(species_observed))
head(ebird)
habitat <- read_csv(file.path(base_data_dir, "eb_tutorial_data/data/pland-elev_location-year.csv")) %>% 
   mutate(year = as.integer(year))
head(habitat)

# combine ebird and modis data
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))

pred_surface <- read_csv(file.path(base_data_dir, "eb_tutorial_data/data/pland-elev_prediction-surface.csv"))
# latest year of landcover data
max_lc_year <- pred_surface$year[1]
r <- raster(file.path(base_data_dir, "eb_tutorial_data/data/prediction-surface.tif"))

map_proj <- st_crs("ESRI:102003")
map_proj
ne_land <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_land
bcr <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_country_lines <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
ne_state_lines <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
# st_geometry is at the end so that the data is converted to sf format.
ne_state_lines
         
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5,
                         year == max(year))
ebird_filtered

occ <- filter_repeat_visits(ebird_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))

occ
occ$pland_13_urban

occ_wide <- format_unmarked_occu(occ, 
                                 site_id = "site", 
                                 response = "species_observed",
                                 site_covs = c("n_observations", 
                                               "latitude", "longitude", 
                                               "pland_04_deciduous_broadleaf", 
                                               "pland_05_mixed_forest",
                                               "pland_12_cropland",
                                               "pland_13_urban"),
                                 obs_covs = c("time_observations_started", 
                                              "duration_minutes", 
                                              "effort_distance_km", 
                                              "number_observers", 
                                              "protocol_type",
                                              "pland_04_deciduous_broadleaf", 
                                              "pland_05_mixed_forest"))
occ_wide
occ_wide$time_observations_started.1

dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id for each site
occ_wide_cell <- occ_wide %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
occ_wide_cell

dim(occ_wide_cell)
names(occ_wide_cell)
point_counts <- occ_wide_cell %>% group_by(cell) %>% summarise(count=n())
point_counts
dim(point_counts)
colnames(point_counts)


# sample one site per grid cell
occ_ss <- occ_wide_cell %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  dplyr::select(-cell)
# calculate the percent decrease in the number of sites
1 - nrow(occ_ss) / nrow(occ_wide)

occ_um <- formatWide(occ_ss, type = "unmarkedFrameOccu")
summary(occ_um)

occ_model <- occu(~ time_observations_started + 
                    duration_minutes + 
                    effort_distance_km + 
                    number_observers + 
                    protocol_type +
                    pland_04_deciduous_broadleaf + 
                    pland_05_mixed_forest
                  ~ pland_04_deciduous_broadleaf + 
                    pland_05_mixed_forest + 
                    pland_12_cropland + 
                    pland_13_urban, 
                  data = occ_um)

pred_surface
occ_pred <- predict(occ_model, 
                    newdata = as.data.frame(pred_surface), 
                    type = "state")

occ_pred

pred_occ <- bind_cols(pred_surface, 
                      occ_prob = occ_pred$fit, 
                      occ_se = occ_pred$se.fit) %>% 
  select(latitude, longitude, occ_prob, occ_se)
pred_occ

