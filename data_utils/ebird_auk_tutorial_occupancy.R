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
View(ebird_filtered)

occ <- filter_repeat_visits(ebird_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))

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

# Working to make a grid
dgcellstogrid(point_counts)

ggplot() +
  geom_polygon(data=occ_wide_cell,  aes(x=longitude, y=latitude, group=group), fill=NA, color="black")   +
  scale_fill_gradient(low="blue", high="red")+
  geom_sf(data=hgrids[[1]], fill=NA, color="#1B9E77")+
  geom_sf(data=hgrids[[2]], fill=NA, color="#D95F02")+
  geom_sf(data=hgrids[[3]], fill=NA, color="#7570B3")+
  # coord_sf(crs="+proj=ortho +lat_0=0 +lon_0=90")+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())

ggplot() +
  geom_polygon(data=occ_wide_cell,  aes(x=longitude, y=latitude))   +
  scale_fill_gradient(low="blue", high="red")+
  geom_sf(data=hgrids[[1]], fill=NA, color="#1B9E77")+
  xlab('')+ylab('')+
  theme(axis.ticks.x=element_blank())+
  theme(axis.ticks.y=element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.text.y=element_blank())


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

##############################
# Now start working on occupancy with my own data.
# I need to get the covariate data (2019 ideally) and the SO data in the same year (MAKE SURE THERE IS NUTHATCH DATA FOR THE YEAR)
# Also need to work with single season.
# Where is that species_observed variable
map_proj <- st_crs("ESRI:102003")
ebd_nh <- auk_ebd(file.path(base_data_dir, "ebd_bnhnut_smp_relJun-2023/ebd_bnhnut_smp_relJun-2023.txt"))
ebd_nh %>% auk_date(date = c("2019-01-01", "2019-12-31")) %>% 
  auk_complete() -> ebd_nh_filtered

auk_filter(ebd_nh_filtered, file = file.path(base_data_dir, "ebd_bnhnut_smp_relJun-2023/nuthatch_filtered_2019.txt"), overwrite=T)
nuthatch <- read_ebd(file.path(base_data_dir, "ebd_bnhnut_smp_relJun-2023/nuthatch_filtered_2019.txt"))

read_sf(file.path(base_data_dir, "us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg"))
state_bound <- read_sf(file.path(base_data_dir, "us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg")) %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
attributes(state_bound)
state_bound
plot(state_bound)
st_crs(state_bound)

nuthatch_sf <- nuthatch %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(crs="ESRI:102003") %>% st_geometry()
st_crs(nuthatch_sf) == st_crs(state_bound) 

par(mar = c(0.25, 0.25, 0.25, 0.25))
plot(state_bound)
plot(nuthatch_sf, pch = 19, cex = 0.1, col = alpha("#555555", 0.25), add = TRUE)

# Then plot over the state
# Set up a grid first
state_grid <- state_bound %>% st_make_grid(cellsize=10000, what = "polygons", crs = "ESRI:102003")
state_grid
# get points just in state
state_pp <- st_intersection(nuthatch_sf, state_bound)
# Using [] to select the cells is the way to go
subgrid <- state_grid[state_bound]
plot(state_bound)
plot(state_pp, pch = 19, cex = 0.5, col = alpha("orange", 0.5), add = TRUE)
plot(subgrid, col = alpha("black", 0.0001), add=T)
#plot(state_grid, col = alpha("black", 0.0001), add=T)

# More advanced using tidyterra and ggplot
p <- ggplot() + 
  #geom_spatvector(data=state_grid, fill = 'transparent', colour="lightblue") +
  geom_sf(data = state_bound, color=alpha("white",0.9)) + 
  geom_sf(data=state_pp, color=alpha("#FFC81C",0.2), size=0.5)+
  ggtitle("Brown-headed Nuthatch observations in Tennessee") +
  theme_minimal() +
  theme(text=element_text(size=10)) +
  coord_sf()
print(p)
fig_name="data/ebird/ebd_po_nuthatch.png"
# To save the degree symbol, we use cairo: 
# https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
# This is also used for embedding custom fonts in pngs/svgs in general
# To show it in R studio, change graphics backend to cairo
ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")


# Then create an intensity map of the points per cell
point_counts <- st_intersects(subgrid, state_pp, sparse=F)
dim(point_counts)
# Count the number of points in each cell
pp_counts <- apply(point_counts, MARGIN=1, FUN=sum)
length(pp_counts)
which(pp_counts>0)

# Let's convert this to a data frame which is how we will require the data for stan
pp_df <- data.frame(counts=pp_counts)
subgrid <- st_sf(subgrid)
subgrid$counts <- pp_df$counts
class(st_geometry(subgrid))
class(st_geometry(subgrid)[[1]])
subgrid$counts <- pp_counts

# We can plot the subgrid (either as sfc or sf) with spatvector
# make_grid produces sfc but st_sf turns this to sf, to which 
# we can add attributes such as point counts
ggplot() + 
  geom_spatvector(data=subgrid, aes(fill=counts), color=alpha("#B6B6B6", 0.5)) +
  geom_sf(data=state_pp, color=alpha("lightgrey",0.2))+
  #geom_spatvector(data=state_grid, fill = 'transparent', colour="lightblue") +
  geom_sf(data = state_bound, color=alpha("white",0.9), fill='transparent') + 
  scale_fill_viridis_c(begin=0.2, end=1, option="magma",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

subgrid$counts
st_geometry(subgrid)
# That sets up our grid and PO data properly
# Now prepare the landcover covariate
landcover_filename <- file.path(base_data_dir, "copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif")
file.exists(landcover_filename)
lc_se_us <- rast(landcover_filename) 
# Project state to raster crs
# Then crop raster in crs of raster, and convert raster back to original state crs
max(lc_se_us)
unique_rasts_val <- unique(lc_se_us)
unique_rasts_val
dim(unique_rasts_val)
prj_state <- terra::project(vect(state_bound), lc_se_us)
plot(prj_state)
crs(prj_state)
crs(lc_se_us)
crop_lc_rast <- crop(lc_se_us, prj_state)
plot(crop_lc_rast)
# projected_lc <- terra::project(crop_lc_rast, crs(vect(state_bound)))
# unique(projected_lc)
# max(projected_lc)
state_pp_prj <- terra::project(vect(state_pp), lc_se_us)


#ggplot() + 
#  geom_spatraster(data=projected_lc) +
#  geom_sf(data=state_pp_prj, color=alpha("orange", 0.9))+
#  #geom_spatvector(data=state_grid, fill = 'transparent', colour="lightblue") +
#  geom_sf(data = state_bound, color=alpha("white",0.9), fill='transparent') + 
#  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
#  theme_minimal()+
#  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

ggplot() + 
  geom_spatraster(data=crop_lc_rast) +
  geom_sf(data=state_pp_prj, color=alpha("orange", 0.9))+
  #geom_sf(data = state_bound, color=alpha("white",0.9), fill='transparent') + 
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.7) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

# Now we need to do some processing of the MODIS data to make it a bit more suitable for occupancy modelling
# "approximately 2.5 km by 2.5 km neighborhood (5 by 5 MODIS cells) centered on the checklist location is 
# sufficient to account for the spatial precision in the data when the maximum distance of travelling counts has been limited to 5 km"
neighborhood_radius <- 5 * ceiling(max(res(crop_lc_rast))) / 2
neighborhood_radius

nuthatch_unique <- nuthatch %>% distinct(locality_id, latitude, longitude) 

nuthatch_buffer <- nuthatch_unique %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326) %>% 
  # transform to state crs projection
  st_transform(crs = crs(state_pp)) %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius) 
  
class(st_as_sf(prj_state))
sf_prj_state <- st_as_sf(prj_state)
state_border <- getBorders(sf_prj_state)

# Do the intersection in the crs of the state
nh_state_buff <- st_intersection(nuthatch_buffer, state_bound)
nh_state_buff
plot(nh_state_buff)

# TODO: Working on calculating pland from the buffers
# Then go back to crs of raster
buff_prj <- st_transform(nh_state_buff, crs=crs(crop_lc_rast))
lc_ext <- exact_extract(crop_lc_rast, buff_prj)
lc_ext
lc_ext[1:2] %>% lapply(function(x) head(x))

require(reshape2)
mapped_lc <- map(lc_ext, dplyr::count, value)
mapped_lc
rbind(mapped_lc)
melt(mapped_lc)
mapped_lc
tibble(nh_state_buff, mapped_lc)

lc_group <- lc_ext %>% group_by(ID)
 
pland <- lc_ext %>% 
# calculate proporiton
group_by(locality_id) %>% 
mutate(pland = n / sum(n)) %>% 
ungroup() %>% 
select(-n) %>% 
# remove NAs after tallying so pland is relative to total number of cells
filter(!is.na(landcover))                            

# Then prepare occupancy data for modelling
# combine ebird and modis data
# We will first need to join the covariate data to the ebird data using some sort of site id.
# It is also possible to do this just by indexing since I am aggregating the counts to grid cells anyway
# TODO: Create site id somehow
ebird_habitat <- inner_join(ebird, habitat, by = c("locality_id", "year"))
# Then filter based on a few criteria I need to check
# so that the data are closed PO data
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5,
                         year == max(year))

occ <- filter_repeat_visits(ebird_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))

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
