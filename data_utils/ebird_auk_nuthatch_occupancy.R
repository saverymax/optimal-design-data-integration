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
library(AICcmodavg)
library(stringr)

select <- dplyr::select

base_data_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\"
ebd_download_dir <- "ebd_US_bnhnut_201901_201912_smp_relJul-2023"
##############################
# Now start working on occupancy with my own data.
# I need to get the covariate data (2019 ideally) and the SO data in the same year (MAKE SURE THERE IS NUTHATCH DATA FOR THE YEAR)
# Also need to work with single season.
# Where is that species_observed variable
map_proj <- st_crs("ESRI:102003")
ebd_nh <- auk_ebd(file.path(base_data_dir, ebd_download_dir, "ebd_US_bnhnut_201901_201912_smp_relJul-2023.txt"),
                  file_sampling = file.path(base_data_dir, ebd_download_dir, "ebd_US_bnhnut_201901_201912_smp_relJul-2023_sampling.txt"))
ebd_nh %>% auk_date(date = c("2019-01-01", "2019-12-31")) %>% 
  auk_complete() -> ebd_nh_filtered
auk_filter(ebd_nh_filtered, file = file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_2019.txt"), 
                            file_sampling=file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_2019_sampling.txt"), overwrite=T) 
nuthatch_obs <- read_ebd(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_2019.txt"))
# TODO: Make plots only using the PO data and not the merged sampling + PO data
nuthatch_sampling <- read_sampling(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_2019_sampling.txt"))
nuthatch_obs
nuthatch_sampling

# This allows us to combine the 2 sets. 
nuthatch_zf <- auk_zerofill(nuthatch_obs, nuthatch_sampling, collapse = TRUE)
# Some 4 million rows
nuthatch_zf

# Do some preprocessing of the data so that it is easier for modelling.
# This follows https://cornelllabofornithology.github.io/ebird-best-practices/ebird.html exactly
time_to_decimal <- function(x) {
  x <- hms(x, quiet = TRUE)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

# clean up variables
nuthatch_zf_filter1 <- nuthatch_zf %>% 
  mutate(
    # convert X to NA
    observation_count = if_else(observation_count == "X", 
                                NA_character_, observation_count),
    observation_count = as.integer(observation_count),
    # effort_distance_km to 0 for non-travelling counts
    effort_distance_km = if_else(protocol_type != "Traveling", 
                                 0, effort_distance_km),
    # convert time to decimal hours since midnight
    time_observations_started = time_to_decimal(time_observations_started),
    # split date into year and day of year
    year = year(observation_date),
    day_of_year = yday(observation_date)
  )

# additional filtering
nuthatch_zf_filter2 <- nuthatch_zf_filter1 %>% 
  filter(
    # effort filters
    duration_minutes <= 5 * 60,
    effort_distance_km <= 5,
    # last 10 years of data
    year >= 2010,
    # 10 or fewer observers
    number_observers <= 10)

nuthatch <- nuthatch_zf_filter2 %>% 
  select(checklist_id, observer_id, sampling_event_identifier,
         scientific_name,
         observation_count, species_observed, 
         state_code, locality_id, latitude, longitude,
         protocol_type, all_species_reported,
         observation_date, year, day_of_year,
         time_observations_started, 
         duration_minutes, effort_distance_km,
         number_observers)
# Save that csv for later use
write_csv(nuthatch, file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_for_occ.csv"), na = "")
# Read that data (useful to skip prev steps)
nuthatch <- read_csv(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_for_occ.csv")) %>% 
  mutate(year = year(observation_date),
         # occupancy modeling requires an integer response
         species_observed = as.integer(species_observed))

state_bound <- read_sf(file.path(base_data_dir, "us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg")) %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
attributes(state_bound)
state_bound
plot(state_bound)
st_crs(state_bound)

# TODO: There's a problem because now that I have the PA data all the processing is much more intensive than with the PO only.
# Convert to sf for spatial processing
nuthatch_sf <- nuthatch %>% 
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(crs="ESRI:102003") %>% st_geometry()
st_crs(nuthatch_sf) == st_crs(state_bound) 

par(mar = c(0.25, 0.25, 0.25, 0.25))
plot(state_bound)
plot(nuthatch_sf, pch = 19, cex = 0.1, col = alpha("#555555", 0.25), add = TRUE)

# Then plot over the state
# Set up a grid first
# Shoudl be 2500 x 2500 meters
crs(state_bound)
state_grid <- state_bound %>% st_make_grid(cellsize=c(2500,2500), what = "polygons", crs = "ESRI:102003")
state_grid
# get points in state
# One way to do it
# state_pp_within <- st_within(nuthatch_sf, state_bound, prepared = T, sparse=F)
# state_pp <- nuthatch_sf[state_pp_within]
state_pp <- st_intersection(nuthatch_sf, state_bound)
class(state_pp)
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
  ggtitle("Brown-headed Nuthatch checklists in Tennessee") +
  theme_minimal() +
  theme(text=element_text(size=10)) +
  coord_sf()
print(p)
fig_name="data/ebird/ebd_checklist_nuthatch.png"
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
# See https://land.copernicus.eu/global/products/lc and go the the viewer.
# The documentation is also on this page ^
# docs: https://land.copernicus.eu/global/sites/cgls.vito.be/files/products/CGLOPS1_PUM_LC100m-V3_I3.4.pdf
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
lc_mask <- terra::mask(lc_se_us, prj_state)
plot(lc_mask)
# projected_lc <- terra::project(crop_lc_rast, crs(vect(state_bound)))
# unique(projected_lc)
# max(projected_lc)
state_pp_prj <- terra::project(vect(state_pp), lc_se_us)

# I am also gonig to add some additional information: The Vegetation index from  MODIS
# Different layers from this download are described here: https://lpdaac.usgs.gov/products/mcd12q2v061/
# See data readme for description of how files are obtained
# See https://www.usgs.gov/landsat-missions/landsat-enhanced-vegetation-index#:~:text=EVI%20is%20similar%20to%20Normalized,in%20areas%20with%20dense%20vegetation.
# for info about evi
modis_evi_filename <- file.path(base_data_dir, "modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif")
file.exists(modis_evi_filename)
modis_evi_rast <- rast(modis_evi_filename) 
modis_evi_rast
plot(modis_evi_rast)
# Project full raster then crop.
evi_rast_prj <- terra::project(modis_evi_rast, lc_se_us)
crop_evi_rast <- crop(evi_rast_prj, prj_state)

# Plot evi
p <- ggplot() + 
  geom_spatraster(data=crop_evi_rast) +
  geom_sf(data=state_pp_prj, color=alpha("black", 0.9))+
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.4) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch and EVI in Tennessee")
print(p)
fig_name="data/ebird/evi.png"
ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")


# Plot landcover
p <- ggplot() + 
  geom_spatraster(data=crop_lc_rast) +
  geom_sf(data=state_pp_prj, color=alpha("orange", 0.9))+
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.4) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch and landcover in Tennessee")
print(p)
fig_name="data/ebird/landcover.png"
ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")


# Now we need to do some processing of the MODIS data to make it a bit more suitable for occupancy modelling
# "approximately 2.5 km by 2.5 km neighborhood (5 by 5 MODIS cells) centered on the checklist location is 
# sufficient to account for the spatial precision in the data when the maximum distance of travelling counts has been limited to 5 km"
neighborhood_radius <- 5000 * ceiling(max(res(crop_lc_rast))) / 2
neighborhood_radius

nuthatch_unique <- nuthatch %>% distinct(locality_id, latitude, longitude) 
nuthatch_unique

nuthatch_unique_state <- nuthatch_unique %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326) %>% 
  # transform to state crs projection
  st_transform(crs = crs(state_bound)) %>% 
  st_intersection(state_bound)
nuthatch_unique_state

nuthatch_buffer <- nuthatch_unique_state %>% 
  # buffer to create neighborhood around each point
  st_buffer(dist = neighborhood_radius)
nuthatch_buffer

#class(prj_state)
#sf_prj_state <- st_as_sf(prj_state)
#state_border <- getBorders(sf_prj_state)
#getBorders(sf_prj_state)
# :(

# Take a look at the buff on -headed Nuthatch  in Tennessee the map
# The buffer will have all sampled sites, not just observed PO sites
buff_prj <- st_transform(nuthatch_buffer, crs=crs(crop_lc_rast))
crs(buff_prj)
p <- ggplot() + 
  geom_spatraster(data=crop_lc_rast) +
  geom_sf(data=buff_prj, color=alpha("black", 1), linewidth=0.5)+
  geom_sf(data=state_pp_prj, color=alpha("orange", 0.5))+
  #geom_sf(data = state_bound, color=alpha("white",0.9), fill='transparent') + 
  geom_sf(data = prj_state, color=alpha("white",0.9), fill='transparent', linewidth=0.4) + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Buffer on observed sites")
print(p)
fig_name="data/ebird/site_buffer.png"
ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")

# Create the fractional covariates
lc_ext_frac <- exact_extract(crop_lc_rast, buff_prj, "frac")
lc_ext_frac
dim(lc_ext_frac)
max_cols <- max.col(lc_ext_frac)
max_col_names <- colnames(lc_ext_frac)[max_cols]
mode_sites <- as.numeric(str_sub(max_col_names, 6, -1))
# Take mode
lc_ext_mode <- exact_extract(crop_lc_rast, buff_prj, "mode")
# Manual mode from frac should equal mode using function .
stopifnot(lc_ext_mode == mode_sites)
#lc_ext_terra <- terra::extract(crop_lc_rast, vect(buff_prj))

dim(crop_evi_rast)
# Let's do the same for the EVI data
evi_buff <- exact_extract(crop_evi_rast, buff_prj, "mean")
evi_buff
class(evi_buff)
evi_df <- data.frame(evi=evi_buff, locality_id=buff_prj$locality_id)

# Mapping unique sites to covariates using the locality_id
dim(nuthatch_unique_state)
lc_ext_frac$locality_id <- buff_prj$locality_id
names(lc_ext_frac)
#nuthatch_unique_state$landcover <- lc_ext_mode 
nuthatch_unique_covars <- inner_join(buff_prj, lc_ext_frac, by=c("locality_id"))
nuthatch_unique_covars <- inner_join(nuthatch_unique_covars, evi_df, by=c("locality_id"))
dim(nuthatch_unique_covars)
nuthatch_unique_covars
names(nuthatch_unique_covars)
class(nuthatch_unique_covars)
#plot(nuthatch_unique_covars, max.plot09=19)


evi_rast <- terra::rasterize(vect(nuthatch_unique_covars), crop_evi_rast, field="evi")
evi_rast
evi_rast[is.na(evi_rast)] <- 0
plot(evi_rast)

p <- ggplot() + 
  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatvector(data=nuthatch_unique_covars, aes(col=frac_126)) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, name="open_forest_unknown") +
  theme_minimal()+
    ggtitle("Brown-headed Nuthatch Intensity in Tennessee")
print(p)

p <- ggplot() + 
  geom_sf(data=buff_prj, color=alpha("black", 1), linewidth=0.5)+
  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatvector(data=nuthatch_unique_covars, aes(col=frac_40)) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

p <- ggplot() + 
  geom_sf(data=buff_prj, color=alpha("black", 1), linewidth=0.5)+
  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatvector(data=nuthatch_unique_covars, aes(col=frac_30)) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

p <- ggplot() + 
  geom_sf(data=buff_prj, color=alpha("black", 1), linewidth=0.5)+
  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatvector(data=nuthatch_unique_covars, aes(col=frac_114)) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

p <- ggplot() + 
  geom_sf(data=buff_prj, color=alpha("black", 1), linewidth=0.5)+
  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatvector(data=nuthatch_unique_covars, aes(col=frac_116)) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

p <- ggplot() + 
  geom_sf(data=buff_prj, color=alpha("black", 1), linewidth=0.5)+
  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
  geom_spatraster(data=evi_rast) +
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Nuthatch EVI in Tennessee at Nuthatch sites")
print(p)
fig_name="data/ebird/site_evi.png"
ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")




# Then final step to prepare occupancy data for modelling
# combine ebird and modis data
# This will reduce nuthatch to only locality_id in the state
ebird_habitat <- inner_join(nuthatch, nuthatch_unique_covars, by = c("locality_id"))
ebird_habitat
# Then filter based on a few criteria I need to check
# so that the data are closed PO data
ebird_habitat$number_observers
ebird_filtered <- filter(ebird_habitat, 
                         number_observers <= 5)

occ <- filter_repeat_visits(ebird_filtered, 
                            min_obs = 2, max_obs = 10,
                            annual_closure = TRUE,
                            date_var = "observation_date",
                            site_vars = c("locality_id", "observer_id"))

class(occ)
covar_names <- c(colnames(occ)[26:length(colnames(occ))], "latitude", "longitude")
covar_names
select_covar <- covar_names[c(1,2, 7, 8, 9, 10, 11, 13, 14, 15)]
select_covar
covar_labels <- c("herbaceuos_vegetation", "agriculture", "closed_forest_evergreen_needle", 
                  "close_forest_decid_broad", "closed_forest_mixed", 
                  "closed_forest_unknown", "open_forest_evergreen_needle", 
                  "open_forest_decid_needle", "open_forest_unknown")

# Use closed forests for detection prob
occ$protocol_type <- as.factor(occ$protocol_type)
occ_wide <- format_unmarked_occu(occ, 
                                 site_id = "site", 
                                 response = "species_observed",
                                 site_covs = select_covar,
                                 obs_covs = c("duration_minutes", "number_observers",  "effort_distance_km",
                                              "time_observations_started", "protocol_type", 
                                              "frac_111", "frac_114", "frac_116"))

# Do some spatial subsampling to reduce bias
dggs <- dgconstruct(spacing = 5)
# get hexagonal cell id for each site
occ_wide_cell <- occ_wide %>% 
  mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum)
# sample one site per grid cell
occ_ss <- occ_wide_cell %>% 
  group_by(cell) %>% 
  sample_n(size = 1) %>% 
  ungroup() %>% 
  select(-cell)
# calculate the percent decrease in the number of sites
1 - nrow(occ_ss) / nrow(occ_wide)
occ_ss

# Unmarked formatting
occ_um <- formatWide(occ_ss, type = "unmarkedFrameOccu")
summary(occ_um)

occ_model <- occu(~ duration_minutes + frac_111 + frac_114 + frac_116 +
                    number_observers + effort_distance_km + time_observations_started + protocol_type
                  ~ frac_30 + frac_40 + frac_111 + frac_114 + frac_115 + frac_116 + frac_121 +
                    frac_126, 
                  data = occ_um)
summary(occ_model)
occ_gof <- mb.gof.test(occ_model, nsim = 10, plot.hist = T)
occ_gof
occ_gof$chisq.table <- NULL
print(occ_gof)
# Just hacking this together based on already existing data
# Using the proportions for the buffers already calculated
occ_pred <- unmarked::predict(occ_model, 
                    newdata = lc_ext_frac,
                    type = "state")
occ_pred$Predicted
nuthatch_unique_covars
nuthatch_unique_covars$preds <- occ_pred$Predicted

# Prepare to convert to raster
preds_select <- nuthatch_unique_covars %>% select(geometry, preds)
preds_select
class(preds_select)
# We can just plot as vector
plot(preds_select)
preds_select$preds
preds_select$geometry

r_pred <- preds_select %>% 
  # convert to spatial features
    #st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = crs(crop_lc_rast))

# TODO
# This isn't using the prediction probabilities for some reason
r_pred <- terra::rasterize(r_pred, crop_lc_rast)
r_pred
crs(r_pred)
plot(r_pred)
plot(prj_state, add=T)
