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
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(crs="ESRI:102003") %>% st_geometry()
  st_crs(nuthatch_sf) == st_crs(state_bound) 
 
  nuthatch_po_sf <- nuthatch_obs %>% 
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>% st_transform(crs="ESRI:102003") %>% st_geometry()
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

create_pp_grid <- function(){
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
  
  # Then create an intensity map of the points per cell
  point_counts <- st_intersects(subgrid, state_pp, sparse=F)
  dim(point_counts)
  # Count the number of points in each cell
  pp_counts <- apply(point_counts, MARGIN=1, FUN=sum)
  length(pp_counts)
  which(pp_counts>0)
}

generate_covariates <- function(){
  # Function to create covariates based on the grid and buffer.
  neighborhood_radius <- 5000 * ceiling(max(res(crop_lc_rast))) / 2
  neighborhood_radius
  nuthatch_unique <- nuthatch %>% distinct(locality_id, latitude, longitude) 
  # Select only the unique observations within the state
  nuthatch_unique_state <- nuthatch_unique %>% st_as_sf(coords = c("longitude", "latitude"), crs=4326) %>% 
    # transform to state crs projection
    st_transform(crs = crs(state_bound)) %>% 
    st_intersection(state_bound)
  nuthatch_unique_state
  
  nuthatch_buffer <- nuthatch_unique_state %>% 
    # buffer to create neighborhood around each point
    st_buffer(dist = neighborhood_radius)
  nuthatch_buffer
  
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
  #fig_name="data/ebird/site_buffer.png"
  #ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")`
  
  # Create the fractional covariates
  lc_ext_frac <- exact_extract(crop_lc_rast, buff_prj, "frac", progress=F)
  colnames(lc_ext_frac)
  dim(lc_ext_frac)
  max_cols <- max.col(lc_ext_frac)
  max_col_names <- colnames(lc_ext_frac)[max_cols]
  mode_sites <- as.numeric(str_sub(max_col_names, 6, -1))
  # Take mode
  lc_ext_mode <- exact_extract(crop_lc_rast, buff_prj, "mode", progress=F)
  # Manual mode from frac should equal mode using function .
  stopifnot(lc_ext_mode == mode_sites)
  
  # Let's do the same for the EVI data
  evi_buff <- exact_extract(crop_evi_rast, buff_prj, "mean", progress=F)
  evi_df <- data.frame(evi=evi_buff, locality_id=buff_prj$locality_id)
  
  # Then extract from the elevation
  elev_buff <- exact_extract(crop_elev_rast, buff_prj, "mean", progress=F)
  head(elev_buff)
  class(evi_buff)
  elev_df <- data.frame(elev=elev_buff, locality_id=buff_prj$locality_id)
  
  # Mapping unique sites to covariates using the locality_id
  dim(nuthatch_unique_state)
  lc_ext_frac$locality_id <- buff_prj$locality_id
  names(lc_ext_frac)
  #nuthatch_unique_state$landcover <- lc_ext_mode 
  # In this way we will keep the buffer polygons
  nuthatch_unique_covars <- inner_join(buff_prj, lc_ext_frac, by=c("locality_id"))
  nuthatch_unique_covars <- inner_join(nuthatch_unique_covars, evi_df, by=c("locality_id"))
  nuthatch_unique_covars <- inner_join(nuthatch_unique_covars, elev_df, by=c("locality_id"))
  dim(nuthatch_unique_covars)
  nuthatch_unique_covars
  names(nuthatch_unique_covars)
  class(nuthatch_unique_covars)
}
