# Learning how to plot points in R and put a grid over 
# This is nice to do, but it's only the initial start to using real data. 
# What I have to do is have some presence only data, such as that from hans (if it can be considered PO)
# and then also have some covariates that correspond to this data. Then I can use the data as my prior.
# What I am doing now is just plotting some data on a map. I also have a 1km grid that
# goes along with the map, which is nice, but not what I need for the end product.
library(sf)
library(terra)
library(tidyterra)
library(dplyr)
library(ggplot2)

base_data_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\"

shp_file <- file.path(base_data_dir, "us_regions/regions", "GAP_Regions.shp")
us_shp <- terra::vect(shp_file)
crs(us_shp)

# From Reich: 
# We use eBird data (Sullivan et al., 2009) to construct an initial estimate of the species distribution map of the brown-headed nuthatch (BHNU).
ebird_data_path <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\code\\supplement_code_from_papers\\kpacifici-ecodesign-df5b292\\ebdata.csv"
ebird_data <- read.csv(ebird_data_path)
ebird_data %>% dplyr::filter(EB_Count > 0) -> ebird_data
ebird_data$EB_Count
coords <- data.frame(long=ebird_data$Long, lat=ebird_data$Lat, count=ebird_data$EB_Count)
coords
ebird_vec <- vect(coords, geom=c("long", "lat"), crs="EPSG:4055")
ebird_vec
class(ebird_vec)

us_prj <- project(us_shp, "EPSG:4055")
plot(us_prj)
plot(ebird_vec, add=T)

# TODO: Plot with fill based on counts
ggplot() +
  geom_spatvector(data=us_prj) +
  geom_spatvector(data=ebird_vec, aes(fill=count)) + 
  scale_fill_viridis_b() +
  theme_minimal()



# Then look at Belgium
shp_file <- file.path(base_data_dir, "Belgium_shapefile", "be_1km.shp")
shp_file
rast_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\be_map\\"
rast_file <- paste(rast_dir, "be_topo_map.tif", sep="")
# Landcover from lifewatch
rast_file <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\copernicus_land_service\\lifewatch.tif"
rast_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\Belgium_shapefile\\"
#rast_dir <- "C:\\Users\\msavery\\OneDrive - UGent\\Documents\\ghent_phd_spatial_doe\\data\\be_map\\9490a6cf-808e-11e9-9847-482ae30f98d9_tiff_3812\\GeoTiff_L08\\"
#rast_file <- paste(rast_dir, "13map.tif", sep="")
rast_file
p <- terra::vect(shp_file)
p
plot(p)
r <- terra::rast(rast_file)
r
plot(r)
lines(p)
terra::crs(r)
terra::crs(p)

# project p to CRS of r?
# EPSG:3812, raster
# EPSG:3035, polygons
prj_p <- project(p, "EPSG:3812")
plot(r)
lines(prj_p, col="gray", alpha=.5)
pp <- spatSample(prj_p, 10, "random")
pp$be_topo_map
points(pp, cex=1, col="blue")

# or r to p, but it's slower
#prj_r <- project(r, "EPSG:3035")
#prj_r
#plot(prj_r)
#lines(p)

# That's all very interesting.
# Now try to map the butterfly data from https://zenodo.org/record/4004609#.ZF4XZHZBy3B to the sites
living_planet_index <- read.csv(paste(base_data_dir, "butterflies_observations_finalgrid.csv", sep=""))
living_planet_index
living_planet_index$site_id
living_planet_index$count
length(unique(living_planet_index$site_id))

prj_p$CELLCODE[2]
length(prj_p$CELLCODE[2])

sites <- prj_p$CELLCODE
sites





