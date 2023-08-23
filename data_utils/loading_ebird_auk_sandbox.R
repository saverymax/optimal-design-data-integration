library(auk)
library(raster)
library(tidyverse)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(sf)
library(lubridate)
library(dggridR)
library(unmarked)

base_data_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\"
eb_data <- file.path(base_data_dir, "ebp-workshop-data/ebd/ebd_2014-2015_yucatan.txt")
ebd_top <- read_tsv(eb_data, n_max = 10)
head(ebd_top)
dim(ebd_top)

# Set auk data path
# auk_set_ebd_path(file.path(base_data_dir, "ebp-workshop-data/ebd"))
eb_file <- "ebd_2014-2015_yucatan.txt"
file.exists(eb_data)
auk_ebd(eb_file) %>% 
  auk_species("Resplendent Quetzal") %>% 
  auk_country("Guatemala") %>% 
  auk_date(c("2015-06-01", "2015-06-30"))

auk_ebd(eb_data) %>% 
  auk_species("Magnolia Warbler") %>% 
  auk_country("BZ") %>% 
  auk_protocol(c("Traveling", "Stationary")) %>% 
  auk_time(c("5:00", "9:00"))

nuthatch <- auk_ebd(eb_data) %>% 
  auk_species("Resplendent Quetzal") %>% 
  auk_complete() %>% 
  auk_filter(file=file.path(base_data_dir, "ebp-workshop-data/ebd/quetzal.csv"), overwrite=T)

eb_quetzal <- read_ebd(file.path(base_data_dir, "ebp-workshop-data/ebd/quetzal.csv"))
class(eb_quetzal)

# Let's plot the nutchatch data
gt_map <- ne_countries(country = "Guatemala", returnclass = "sf")
plot(st_geometry(gt_map))
gt_grid <- gt_map %>% st_make_grid(cellsize = 1, what = "polygons") %>% st_intersection(gt_map)
plot(gt_grid)
ggplot() + 
geom_sf(data = gt_map) + 
geom_sf(data = gt_grid)

eb_quetzal <- eb_quetzal %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

par(mar = c(0.25, 0.25, 0.25, 0.25))
quetzal_plot <- eb_quetzal %>% st_geometry()
plot(st_geometry(gt_map) , col = "grey80", border = "white")
plot(quetzal_plot, col = alpha("grey20", 0.3), pch = 19, add=T)
gt_grid <- gt_map %>% st_make_grid(cellsize = 1, what = "polygons") #%>% st_intersection(gt_map)
plot(gt_grid, add=T)

# This doesn't work if we intersect grid and map
point_counts <- st_intersects(gt_grid, st_geometry(eb_quetzal), sparse=F)
dim(point_counts)
# Count the number of points in each cell
pp_counts <- apply(point_counts, MARGIN=1, FUN=sum)
pp_counts
which(pp_counts>0)
# It's nice that these easily correspond to the cells starting from the origin

# ggplot way
ggplot() + 
  geom_sf(data = gt_map) + 
  geom_sf(data=quetzal_plot)+
  geom_spatvector(data=gt_grid, aes(fill=pp_counts), alpha=0.2)+
  theme_minimal()


##############
# Following https://cornelllabofornithology.github.io/ebird-best-practices/ebird.html, chapter 2
# https://cornelllabofornithology.github.io/ebird-best-practices/ebird.html#ebird-explore

# Load their csv data directly, since they already did some processing on that.
ebird <- read_csv(file.path(base_data_dir, "eb_tutorial_data/data/ebd_woothr_june_bcr27_zf.csv"))
# Census projection
map_proj <- st_crs("ESRI:102003")
map_proj
ne_land <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "ne_land") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
plot(ne_land)
bcr <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "bcr") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
plot(bcr)
ne_country_lines <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "ne_country_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
plot(ne_country_lines)
ne_state_lines <- read_sf(file.path(base_data_dir, "eb_tutorial_data/data/gis-data.gpkg"), "ne_state_lines") %>% 
  st_transform(crs = map_proj) %>% 
  st_geometry()
# st_geometry is at the end so that the data is converted to sf format.
plot(ne_state_lines)

ebird_sf <- ebird %>% 
  # convert to spatial points
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% 
  st_transform(crs = map_proj) %>% 
  select(species_observed)

# map
par(mar = c(0.25, 0.25, 0.25, 0.25))
# set up plot area
plot(st_geometry(ebird_sf), col = NA)
# contextual gis data
plot(ne_land, col = "#dddddd", border = "#888888", lwd = 0.5, add = TRUE)
plot(bcr, col = "#cccccc", border = NA, add = TRUE)
plot(ne_state_lines, col = "#ffffff", lwd = 0.75, add = TRUE)
plot(ne_country_lines, col = "#ffffff", lwd = 1.5, add = TRUE)
# ebird observations
# not observed
plot(st_geometry(ebird_sf),
     pch = 19, cex = 0.1, col = alpha("#555555", 0.25),
     add = TRUE)
# observed
plot(filter(ebird_sf, species_observed) %>% st_geometry(),
     pch = 19, cex = 0.3, col = alpha("#4daf4a", 1),
     add = TRUE)
# legend
legend("bottomright", bty = "n",
       col = c("#555555", "#4daf4a"),
       legend = c("eBird checklists", "Wood Thrush sightings"),
       pch = 19)
box()
par(new = TRUE, mar = c(0, 0, 3, 0))
title("Wood Thrush eBird Observations\nJune 2010-2019, BCR 27")

crs(ne_land)
# We choose the cell sized based on the units of the projection :)
# This is pretty slow to run
se_us_grid <- ne_land %>% st_make_grid(cellsize=10000,what = "polygons", crs = "ESRI:102003") #%>% st_intersection(ne_land)
plot(se_us_grid, add=T)
se_us_grid[1]
us_cnts <- st_intersects(se_us_grid, st_geometry(ebird_sf), sparse=F)


######## 
# Processing my own downloaded data: Brown-Headed Nuthatch
# file name: ebd_bnhnut_smp_relJun-2023
ebd_nh <- auk_ebd(file.path(base_data_dir, "ebd_bnhnut_smp_relJun-2023/ebd_bnhnut_smp_relJun-2023.txt"))
ebd_nh %>% auk_date(date = c("2019-01-01", "2019-12-31")) %>% 
  auk_complete() -> ebd_nh_filtered

auk_filter(ebd_nh_filtered, file = file.path(base_data_dir, "ebd_bnhnut_smp_relJun-2023/nuthatch_filtered_2019.txt"), overwrite=T)
nuthatch <- read_ebd(file.path(base_data_dir, "ebd_bnhnut_smp_relJun-2023/nuthatch_filtered_2019.txt"))

# Full map of us
map_proj <- st_crs("ESRI:102003")
us_map <- ne_countries(country = "united states of america", returnclass = "sf") %>% st_transform(crs=map_proj)
class(us_map)
crs(us_map)
us_vect <- vect(us_map)
writeVector(us_vect, "us_filetype.shp")
# Or load just one state downloaded from https://apps.nationalmap.gov/downloader/
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
st_crs(state_bound) == st_crs(us_map)
st_crs(us_map)

par(mar = c(0.25, 0.25, 0.25, 0.25))
plot(st_geometry(us_map))
plot(state_bound, add=TRUE)
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

# ggplot way
ggplot() + 
  geom_sf(data = state_bound) + 
  geom_sf(data=state_pp)+
  geom_sf(data=subgrid, alpha=0.1)+
  theme_minimal()

# More advanced using tidyterra
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
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

subgrid$counts
st_geometry(subgrid)

# Let's also add some covariate data. Here is copernicus landcover 
# This was downloaded from 
# https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2019/W100N40/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif
# using the landcover viewer: https://lcviewer.vito.be/2015
landcover_filename <- file.path(base_data_dir, "copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif")
file.exists(landcover_filename)
lc_se_us <- rast(landcover_filename) 
# Project raster to crs we're using here
# Easier to go from vector to the crs of raster but in this case we need the crs in units of m
projected_lc <- terra::project(lc_se_us, state_bound)
class(projected_lc)
ggplot() + 
  geom_spatraster(data=projected_lc) +
  geom_sf(data=state_pp, color=alpha("orange", 0.9))+
  #geom_spatvector(data=state_grid, fill = 'transparent', colour="lightblue") +
  geom_sf(data = state_bound, color=alpha("white",0.9), fill='transparent') + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

state_sf <- st_sf(state_bound)
class(state_sf)
crop_lc_rast <- crop(projected_lc, state_sf)
plot(crop_lc_rast)

ggplot() + 
  geom_spatraster(data=crop_lc_rast) +
  geom_sf(data=state_pp, color=alpha("orange", 0.9))+
  #geom_spatvector(data=state_grid, fill = 'transparent', colour="lightblue") +
  geom_sf(data = state_bound, color=alpha("white"), fill='transparent', linewidth=0.7) + 
  scale_fill_viridis_c(name="landcover", option="viridis",alpha=0.6) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")

# Can we also load some modis data?
# This was downloaded from https://appeears.earthdatacloud.nasa.gov/
# by selecting the landcover type prodcut for 2019.
# We can see more information about the layers from this download here:
# https://lpdaac.usgs.gov/products/mcd12q1v061/
# Different layers can be loaded.
# Which to use for occupancy modelling?
#modis_landcover_filename <- file.path(base_data_dir, "modis_landcover_types/MCD12Q1.061_LC_Prop1_Assessment_doy2019001_aid0001.tif")
modis_landcover_filename <- file.path(base_data_dir, "modis_landcover_types/MCD12Q1.061_LC_Type1_doy2019001_aid0001.tif")
file.exists(modis_landcover_filename)
modis_lc_rast <- rast(modis_landcover_filename) 
#modis_lc_prj <- terra::project(modis_lc_rast, state_bound)
# Do cropping in the crs of raster
state_bound_proj <- terra::project(vect(state_bound), modis_lc_rast)
plot(state_bound_proj)
crop_lc_rast <- crop(modis_lc_rast, state_bound_proj)
plot(crop_lc_rast)
# Project raster back to state
modis_lc_prj <- terra::project(crop_lc_rast, crs(vect(state_bound)))
plot(modis_lc_prj)

ggplot() + 
  geom_spatraster(data=modis_lc_prj) +
  geom_sf(data=state_pp, color=alpha("orange", 0.9))+
  geom_sf(data = state_bound, color=alpha("white",0.9), fill='transparent') + 
  scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
  theme_minimal()+
  ggtitle("Brown-headed Nuthatch Intensity in Tennessee")
