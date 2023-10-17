#############################################################
# Module for processing PO ebird data and covariates for 
# optimal design usage
#############################################################


inital_auk_processing <- function(base_data_dir, ebd_download_dir){
  # Function to generate and save csv as "nuthatch_filtered_for_occ.csv"
  # and "nuthatch_filtered_po_2019.txt"
  
  # See code in ebird_auk_nuthatch_occupancy
  ebd_nh <- auk_ebd(file.path(base_data_dir, ebd_download_dir, "ebd_US_bnhnut_201901_201912_smp_relJul-2023.txt"),
                    file_sampling = file.path(base_data_dir, ebd_download_dir, "ebd_US_bnhnut_201901_201912_smp_relJul-2023_sampling.txt"))
  ebd_nh %>% auk_date(date = c("2019-01-01", "2019-12-31")) %>%  auk_protocol(protocol = c("Stationary", "Traveling")) %>% auk_complete() -> ebd_nh_filtered
  auk_filter(ebd_nh_filtered, file = file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_po_2019.txt"), 
             file_sampling=file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_2019_sampling.txt"), overwrite=T) 
  nuthatch_obs <- read_ebd(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_po_2019.txt"))
  # TODO: Make plots only using the PO data and not the merged sampling + PO data
  nuthatch_sampling <- read_sampling(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_2019_sampling.txt"))
  nuthatch_zf <- auk_zerofill(nuthatch_obs, nuthatch_sampling, collapse = TRUE)
  # Some 4 million rows
  head(nuthatch_zf)
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
}

main_data_handling_oe <- function(map_prj, map_path, base_data_dir, ebd_download_dir, exp_dir, lc_path, modis_path, elev_path, rast_surface_path, cell_size){
  # Function to do the initial covariate/count data processing and create dataframe and PO data matrix
  # This can/should be run ahead of the OE experiments separately.
  
  state_bound <- load_map(map_prj, map_path, base_data_dir)
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
  fig_name=file.path(exp_dir, "nuthatch_checklists.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  crop_lc_rast <- load_landcover(base_data_dir, lc_path, state_bound)
  # Convert crs to the landcover which we will work in for the rest of the script, though the size of the 
  # grid/raster conversion is probably ok for either
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
  fig_name=file.path(exp_dir, "landcover.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
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
  fig_name=file.path(exp_dir, "evi.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
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
  fig_name=file.path(exp_dir, "elevation.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  subgrid <- create_pp_grid(state_bound, crs(crop_lc_rast), cell_size)
  stopifnot(st_crs(subgrid)==st_crs(crop_lc_rast))
  p <- ggplot() + 
    geom_spatvector(data=subgrid, color=alpha("#B6B6B6", 0.5)) +
    #geom_spatvector(data=subgrid, aes(fill=counts), color=alpha("#B6B6B6", 0.5)) +
    geom_sf(data=state_po, color=alpha("#FFC81C",0.5), size=0.5)+
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7) +
    theme_minimal()+
    ggtitle("2.5km by 2.5km grid over Tennessee")
  fig_name=file.path(exp_dir, "tennessee_grid.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  pp_counts <- aggregate_point_counts(subgrid, state_po_prj)
  pp_counts[which(pp_counts>1)]
  hist(pp_counts[which(pp_counts>1)])
  subgrid$counts <- pp_counts
  
  # Next calculate covariate aggregation based on the grid.
  # TODO: Add ebird covariates from the buuffer: number_observers and duration_minutes
  agg_covars <- generate_gridded_covariates(subgrid, crop_lc_rast, crop_evi_rast, crop_elev_rast, state_po_prj)
  subgrid$evi <- agg_covars$evi_agg
  subgrid$elev <- agg_covars$elev_agg
  subgrid_covars <- bind_cols(subgrid, agg_covars$lc_frac)
  
  # The grill will look a bit lopsided because we make the grid in state crs but transform it to raster crs
  p <- ggplot() + 
    geom_spatvector(data=subgrid_covars, aes(fill=counts), color=alpha("#B6B6B6", 0.5)) +
    scale_fill_viridis_c(begin=0.2, end=1, option="magma",alpha=0.7) +
    theme_minimal()+
    ggtitle("Brown-headed Nuthatch aggregated intensity in Tennessee")
  fig_name=file.path(exp_dir, "aggregated_intensity.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  # rasterized version of site evi 
  evi_rast <- terra::rasterize(subgrid_covars, crop_evi_rast, field="evi")
  # rasterized version of site elev 
  elev_site_rast <- terra::rasterize(subgrid_covars, crop_evi_rast, field="elev")
  elev_site_rast[is.na(elev_site_rast)] <- 0
  p <- ggplot() + 
    #  geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
    geom_spatraster(data=evi_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, na.value="white") +
    theme_minimal()+
    ggtitle("Nuthatch EVI in Tennessee at Nuthatch sites")
  fig_name=file.path(exp_dir, "site_evi.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  p <- ggplot() + 
    geom_sf(data = prj_state, color=alpha("black",0.5), fill='transparent', linewidth=0.7) + 
    geom_spatraster(data=elev_site_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, na.value="white") +
    theme_minimal()+
    ggtitle("Elevation in Tennessee at Nuthatch sites")
  fig_name=file.path(exp_dir, "site_elevation.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  covar_names <- c(colnames(subgrid_covars)[5:(length(colnames(subgrid_covars))-1)])
  # List of all covar labels
  covar_labels <- c("shrubs", "herbaceous_vegetation", "agriculture", "urban", "bare/sparse", "water", "wetland", "closed_forest_evergreen_needle", 
                    "close_forest_decid_broad", "closed_forest_mixed", "closed_forest_unknown", "open_forest_evergreen_needle", 
                    "open_forest_decid_broad leaf", "open_forest_unknown")
  names(covar_names) <- covar_labels
  
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
    fig_name <- paste(covar, "_covar_tenn.png", sep="")
    fig_name <- file.path(exp_dir, fig_name)
    ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  }
  
  covar_df <- st_drop_geometry(subgrid_covars)
  # Let's save the centroids for plotting
  site_centroids <- st_centroid(st_geometry(subgrid_covars))
  #intensity_covars <- covar_df[c("evi", "elev", "frac_111", )]
  intensity_covars <- covar_df[c("evi", "elev", "frac_111", "frac_115", "frac_124")]
  bias_covars <- covar_df[c("frac_50")]#"duration_minutes", "number_observers",  "effort_distance_km", "time_observations_started", 
  k_param_intn <- length(colnames(intensity_covars)) 
  k_param_bias <- length(colnames(bias_covars)) 
  # Need to save a file for use as the rast surface for later
  writeRaster(evi_rast, rast_surface_path, overwrite=T) 
  # Save evi as the rast surface
  data_pack <- list(counts=covar_df$counts, intensity_covars=intensity_covars, bias_covars=bias_covars, k_param_intn=k_param_intn, k_param_bias=k_param_bias,
              site_centroids=site_centroids, subgrid=subgrid, state_po_prj=state_po_prj, rast_surface=evi_rast)
  saveRDS(data_pack, file.path(exp_dir, "data_pack.RDS"))
}

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

load_map <- function(map_proj, map_path, base_data_dir){
  # Load region for optimal design to work on.
  state_bound <- read_sf(file.path(base_data_dir, map_path)) %>% 
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

create_pp_grid <- function(state_bound, rast_crs, cell_size){
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

generate_ebird_pa <- function(data_reps, area_a, surface_data, p_0, pp_posterior, n, sites, link){
  # Add a column of 1s to the auxiliary data
  a_intercept <- rep(1, nrow(surface_data))
  surface_data <- cbind(a_intercept, surface_data)
  occupancy_maps <- matrix(nrow=data_reps, ncol=sites)
  Y_detection <- matrix(nrow=data_reps, ncol=sites)
  theta_reps <- matrix(nrow=data_reps, ncol=sites)  
  for (r in 1:data_reps){
    # The posterior values can be either draws from the PP posterior or the expectation
    param_vec <- pp_posterior[r,]
    if (link=="cloglog"){
      #g_theta <- 1 - exp(-exp(as.matrix(surface_data)%*%t(param_vec)))
      #g_theta <- 1 - exp(-exp(as.matrix(surface_data)%*%param_vec))
      g_theta <- 1 - exp(-exp(as.matrix(surface_data)%*%param_vec + log(area_a)))
    }else{
      stop("No other link function implemented")
    }
    stopifnot(sites==length(g_theta))
    site_presence <- rbinom(sites, 1, g_theta)
    theta_reps[r, ] <- g_theta
    occupancy_maps[r, ] <- site_presence
    Y_detection[r, ] <- rbinom(sites, n, p_0*site_presence)
  } 
  return(list(occupancy=occupancy_maps, Y=Y_detection, theta=theta_reps))
}

fit_point_process_ebird <- function(stan_path, exp_dir, sites, area_a, data_reps, intensity_covars, bias_covars, r_po_data, k_i, k_b, pp_diagnostic,
                                    subgrid, rast_surface){
  # Fit the point process model to the PO data
  stan_models_path <- file.path(stan_path, "stan_nuthatch_models.R")
  source(stan_models_path)
  model_path <- "nuthatch_poisson_process.stan"
  model_string <- nuthatch_poisson_process
  write(model_string, model_path)
  data_site_occ = list(N=sites, X=intensity_covars, y=r_po_data$Y[1,], Z=bias_covars, k_i=k_i, k_b=k_b, area_a=area_a)
  model <- cmdstan_model(model_path) 
  print("Fitting Point Process model to PO data")
  fit <- model$sample(data=data_site_occ, seed=13, chains=3, iter_sampling=2000, iter_warmup=500) 
  print(fit$summary())
  all_params <- c("alpha", "gamma", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "delta[1]")
  #all_params <- c("alpha", "gamma", "beta[1]", "beta[2]", "beta[3]", "delta[1]")
  #intensity_params <- c("alpha", "beta[1]", "beta[2]", "beta[3]")
  intensity_params <- c("alpha", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]")
  params_intercept <- c("alpha", "gamma")
  
  if (exp_args$pp_diagnostic == T){
    posterior <- fit$draws(all_params)
    color_scheme_set("mix-blue-pink")
    p_trace <- mcmc_trace(posterior,
                          facet_args = list(nrow = 2, labeller = label_parsed))
    print(p_trace + facet_text(size = 15))
    
    plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
    print(mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[1]")) + plot_title)
    print(mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[2]")) + plot_title)
    print(mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[3]")) + plot_title)
    print(mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[4]")) + plot_title)
    print(mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c("beta[5]")) + plot_title)
    
    plot_title <- ggtitle(paste("Posterior distributions of detection probability, mean and 90% interval"))
    p_post <- mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars = params_intercept) + plot_title
    print(p_post)
    
    plot_title <- ggtitle(paste("Posterior distributions of detection probability, mean and 90% interval"))
    p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", pars = c("delta[1]")) + plot_title
    print(p_post)
    
    print(mcmc_intervals(fit$draws(), pars=all_params))
    print(mcmc_hist(fit$draws(), pars = all_params))
    print(mcmc_pairs(fit$draws(), pars=all_params))
    print(mcmc_scatter(fit$draws(), pars=c('alpha', 'gamma')))
    print(mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[2]')))
    print(mcmc_scatter(fit$draws(), pars=c('beta[2]', 'beta[3]')))
    print(mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[3]')))
    print(mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[4]')))
    print(mcmc_scatter(fit$draws(), pars=c('beta[1]', 'beta[5]')))
  }
    
  # Examine the ppd for y
  generated_yrep <- fit$draws("y_rep", format="matrix")
  yrep_means <- colMeans(generated_yrep)
  # Create new sf for the preds
  subgrid_geo <- st_geometry(subgrid) %>% st_sf()
  subgrid_geo$preds <- yrep_means
  # Then plot the predictions
  pp_pred_rast <- terra::rasterize(vect(subgrid_geo), rast_surface, field="preds")
  p <- ggplot() + 
    geom_spatraster(data=pp_pred_rast) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.8, na.value="white") +
    #geom_sf(data = prj_state, color=alpha("black",0.4), fill='transparent', linewidth=0.7) + 
    theme_minimal()+
    ggtitle("Predicted counts at nuthatch sites")
  fig_name=file.path(exp_dir, "pp_predictions.png")
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  # Take expectation from the PP params
  pp_posterior <- fit$draws(intensity_params, format="matrix")
  pp_posterior <- matrix(rep(colMeans(pp_posterior), data_reps), nrow=data_reps, byrow = T)
  print("Estimates for intensity (not bias) from pp model")
  print(pp_posterior)
  saveRDS(pp_posterior, file.path(exp_dir, "pp_posterior_ebird.RDS"))
}
  