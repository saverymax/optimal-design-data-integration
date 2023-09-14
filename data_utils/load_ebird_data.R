#############################################################
# Module for processing PO ebird data and covariates for 
# optimal design usage
#############################################################

load_ebird <- function(base_data_dir, ebd_download_dir, state_bound, map_prj){
  # Load pre-processed PO data
  nuthatch_obs <- read_ebd(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_po_2019.txt"))
  nuthatch <- read_csv(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_for_occ.csv")) %>% 
    mutate(year = year(observation_date),
           # occupancy modeling requires an integer response
           species_observed = as.integer(species_observed))
  
  # Crs is 4326 via ebird
  # Convert to sf for spatial processing
  # Separate files for the pure PO data and the data merged with checklists so that absences are included.
  nuthatch_sf <- nuthatch %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(crs="ESRI:102003") #%>% st_geometry()
  st_crs(nuthatch_sf) == st_crs(state_bound) 
 
  nuthatch_po_sf <- nuthatch_obs %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(crs="ESRI:102003") #%>% st_geometry()
  st_crs(nuthatch_po_sf) == st_crs(state_bound) 
  
  # Get the points within the state
  state_pa <- st_intersection(nuthatch_sf, state_bound)
  state_po <- st_intersection(nuthatch_po_sf, state_bound)
  
  return(list(state_pa=state_pa, state_po=state_po))
}

load_map <- function(map_proj, map_path){
  # Load region for optimal design to work on.
  state_bound <- read_sf(map_path) %>% 
    st_transform(crs = map_proj) %>% 
    st_geometry()
  return(state_bound)
  
}

load_landcover <- function(base_data_dir, lc_path, state_bound){
  # Load preprocessed landcover raster
  landcover_filename <- file.path(base_data_dir, lc_path)
  stopifnot(file.exists(landcover_filename))
  lc_se_us <- rast(landcover_filename) 
  # Project state to raster crs
  # Then crop raster in crs of raster, and convert raster back to original state crs
  prj_state <- terra::project(vect(state_bound), lc_se_us)
  crop_lc_filename <- file.path(base_data_dir, "copernicus_landcover/cropped_lc.tif")
  if (!file.exists(crop_lc_filename)){
    crop_lc_rast <- crop(lc_se_us, prj_state)
    writeRaster(crop_lc_rast, crop_lc_filename)
  }else{
    crop_lc_rast <- rast(crop_lc_filename)
  }
  return(crop_lc_rast)
}

load_evi <- function(base_data_dir, modis_path, prj_state, crop_lc_rast){
  # Load preprocessed EVI
  modis_evi_filename <- file.path(base_data_dir, modis_path)
  stopifnot(file.exists(modis_evi_filename))
  modis_evi_rast <- rast(modis_evi_filename) 
  # Project full raster then crop.
  crop_evi_filename <- file.path(base_data_dir, "modis_landcover_dynamics/cropped_evi.tif")
  if(!file.exists(crop_evi_filename)){
    evi_rast_prj <- terra::project(modis_evi_rast, crop_lc_rast)
    crop_evi_rast <- crop(evi_rast_prj, prj_state)# Plot evi
    writeRaster(crop_evi_rast, crop_evi_filename)
  }else{
    crop_evi_rast <- rast(crop_evi_filename)
  }
  return(crop_evi_rast)
}

load_elevation <- function(base_data_dir, elev_path, prj_state){
  #elev_filename <- file.path(base_data_dir, "elevation_earth_env/elevation_1KMmd_GMTEDmd.tif")
  elev_filename <- file.path(base_data_dir, elev_path)
  stopifnot(file.exists(elev_filename))
  elev_rast <- rast(elev_filename)
  crop_elev_filename <- file.path(base_data_dir, "elevation_aster/cropped_elevation.tif")
  if(!file.exists(crop_elev_filename)){
    crop_elev_rast <- crop(elev_rast, prj_state)# Plot evi
    writeRaster(crop_elev_rast, crop_elev_filename)
  }else{
    crop_elev_rast <- rast(crop_elev_filename)
  }
  return(crop_elev_rast)
}

create_pp_grid <- function(state_bound, rast_crs){
  cell_size <- 10000
  # WOrk in the CRS with units of meters that state_bound was put in
  state_grid <- state_bound %>% st_make_grid(cellsize=c(cell_size,cell_size), what = "polygons", crs = "ESRI:102003")
  # Using [] to select the cells is the way to go
  subgrid <- state_grid[state_bound]
  # Then project back to the raster crs and convert to sf object from sfc
  subgrid_prj <- st_transform(subgrid, rast_crs) %>% st_sf()
  return(subgrid_prj)
  
}

aggregate_point_counts <- function(subgrid, state_po_prj){
  # At the moment, the aggregating is done in the crs of the raster
  # Then create an intensity map of the points per cell
  point_counts <- st_intersects(subgrid, state_po_prj, sparse=F)
  # Count the number of points in each cell
  pp_counts <- apply(point_counts, MARGIN=1, FUN=sum)
  return(pp_counts)
}

generate_gridded_covariates <- function(subgrid, crop_lc_rast, crop_evi_rast, crop_elev_rast, state_po_prj){
  # Function to create covariates based on the grid.
  # Create the fractional covariates
  lc_ext_frac <- exact_extract(crop_lc_rast, subgrid, "frac", progress=F)
  colnames(lc_ext_frac)
  dim(lc_ext_frac)
  max_cols <- max.col(lc_ext_frac)
  max_col_names <- colnames(lc_ext_frac)[max_cols]
  mode_sites <- as.numeric(str_sub(max_col_names, 6, -1))
  # Take mode
  lc_ext_mode <- exact_extract(crop_lc_rast, subgrid, "mode", progress=F)
  # Manual mode from frac should equal mode using function .
  stopifnot(lc_ext_mode == mode_sites)
  # Let's do the same for the EVI data
  evi_buff <- exact_extract(crop_evi_rast, subgrid, "mean", progress=F)
  # Then extract from the elevation
  elev_buff <- exact_extract(crop_elev_rast, subgrid, "mean", progress=F)
  # Get ebird covariates buffered; didn't want to work.
  #observers_buff <- subgrid %>% st_intersection(state_po_prj) %>% group_by(id) %>% summarise(avg_observers=mean(number_observers))
  
  return(list(lc_frac=lc_ext_frac, evi_agg=evi_buff, elev_agg=elev_buff))
}
