# Functions for experimental design

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

# Brier based on previous data structures
brier_score <- function(z, z_obs){
  b <- 0
  sites <- ncol(z)
  iter <- nrow(z)
  for (i in 1:sites){
    z_bar_i <- sum(z[, i]) / iter
    if(z_obs[i] == 0){
      b <- b + (z_bar_i)^2 
    }else{
      b <- b + (1-z_bar_i)^2
    }
  }
  return(b/sites)
}

brier_score_stan <- function(z, z_obs){
  # Brier score compatible with stan output
  sites <- length(z)
  brier <- sum((z_obs - z)^2) / sites
  return(brier)
}

design_criteria <- function(criteria, sim_occ, obs_occ){
  if (criteria=="brier"){
    score <- brier_score_stan(sim_occ, obs_occ)
  }
  else{
    print("No other methods implemented currently")
  }
  return(score)
}

estimate_v <- function(model, n_surveys, data_reps, m, sites, sampling_surface, 
                       select_idx, select_sites, r_survey_data, r_po_data,
                       p_logging, params, generated_vars, model_selection, mcmc_iter, debug_stan){
  
  estimate_mat <- matrix(nrow=data_reps, ncol=1)
  for (r in 1:data_reps){
    # For each rth dataset get the m randomly chosen sites
    selected_data <- r_survey_data$Y[r, select_idx]
    #selected_occ <- r_survey_data$occupancy[r, select_idx]
    occ <- r_survey_data$occupancy[r, ]
    # We use all po data for given r
    selected_po <- r_po_data$Y[r,]
    # Here I use the same X covariate for the presence-only data as used for the survey data. The difference is that
    # the survey data here is a subset (select_sites) of the sites, whereas the presence only data 
    # needs covariates for the whole grid to approximate the expected count in the entire region.
    if (model_selection==3){
      data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, X=select_sites$aux_x, Y=selected_data, PO=selected_po, 
                           X_po=sampling_surface$aux_x, Z_po=sampling_surface$aux_z, model_diag=0)
    }
    else if(model_selection==1){
      data_site_occ = list(n_surveys=n_surveys, n_sites=m, total_sites=sites, X=select_sites$aux_x, X_all=sampling_surface$aux_x, Y=selected_data)
    }
    else if(model_selection==4){
      data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, total_sites=sites, Y=selected_data)
    }
    else if(model_selection==5){
      data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, Y=selected_data, PO=selected_po)
    }
    else{
      stop("Other models implementation needs to be checked")
    }
    # For debugging quiet_model can be turned off
    # refresh=0 turns off messages except errors from stan
    # quiet function silences stan output 
    if (debug_stan==F){
      fit <- quiet(model$sample(data=data_site_occ, seed=13, chains=1, 
                                iter_sampling=mcmc_iter, iter_warmup=100, refresh=0, show_messages=F))
    }
    else{
      fit <- model$sample(data=data_site_occ, seed=13, chains=1, 
                                iter_sampling=mcmc_iter, iter_warmup=100)
    }
    if (p_logging==T){
      print("Logging posterior")
      posterior <- fit$draws()
      # TODO: Why is estimate for alpha now rather large??
      print(fit$summary(variables=params))
      
      color_scheme_set("mix-blue-pink")
      p_trace <- mcmc_trace(posterior,  pars = params,
                            facet_args = list(nrow = 2, labeller = label_parsed))
      print(p_trace + facet_text(size = 15))
      
      plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
      p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = params) + plot_title
      print(p_post)
      
      #if(model_selection==3){
      #  plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
      #  p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = c("gamma", "delta")) + plot_title
      #  print(p_post)
      #}
      
      #plot_title <- ggtitle(paste("Posterior distribution of detection probability, mean and 90% interval"))
      #p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", pars = c("p")) + plot_title
      #print(p_post)
      
      # Check issues with divergences
      color_scheme_set("darkgray")
      nuts_fit <- nuts_params(fit)
      diverge_p <- mcmc_parcoord(posterior, pars = params, np = nuts_fit, alpha=.1, 
                                 np_style=parcoord_style_np(div_alpha=1, div_size=.5))
      print(diverge_p)
      
      ppd_count_df <- data.frame(x=sampling_surface$x, y=sampling_surface$y, 
                                 occ_prob=fit$summary(variables=generated_vars[1])$mean)
      
      p <- ggplot(ppd_count_df, aes(x, y, fill=occ_prob)) + 
        geom_tile() +
        scale_fill_viridis(discrete=FALSE) +
        ggtitle("Generated occupancy probability per site")
      print(p)
      
    }
    # Average estimate after R iterations through the datasets.
    # Need to select_idx the sites since it generates for all of them.
    # The 2 index is the ppd for Z
    gen_occupancy <- fit$summary(variables=generated_vars[2])$mean
    print("brier comp")
    print(gen_occupancy)
    print(occ)
    brier <- sum((occ - gen_occupancy)^2) / sites
    # Take the mean of the generated quantity for the posterior estimate
    estimate_mat[r,1] <- design_criteria(criteria="brier", sim_occ=gen_occupancy, obs_occ=occ)
  }
  return(estimate_mat)
}
  

estimate_v_parallel <- function(combined_df, model, n_surveys, m, sites, sampling_surface, 
                       select_idx, select_sites, generated_vars, model_selection, mcmc_iter){
  
  # For each rth dataset get the m randomly chosen sites for the occupancy data, Y surveys, and PO.
  #selected_occ <- combined_df[select_idx]
  occ <- combined_df[1:sites]
  selected_data <- combined_df[sites+select_idx]
  PO_data <- combined_df[(2*sites+1):length(combined_df)] 
  # Here I use the same X covariate for the presence-only data as used for the survey data. The difference is that
  # the survey data here is a subset (select_sites) of the sites, whereas the presence only data 
  # needs covariates for the whole grid to approximate the expected count in the entire region.
  if (model_selection==3){
    data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, X=select_sites$aux_x, Y=selected_data, PO=PO_data, 
                         X_po=sampling_surface$aux_x, Z_po=sampling_surface$aux_z, model_diag=0)
  }
  else if(model_selection==1){
    data_site_occ = list(n_surveys=n_surveys, n_sites=m, total_sites=sites, X=select_sites$aux_x, X_all=sampling_surface$aux_x, Y=selected_data)
  }
  else if(model_selection==4){
    data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, total_sites=sites, Y=selected_data)
  }
  else if(model_selection==5){
    data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, Y=selected_data, PO=PO_data)
  }
  else{
    stop("Other models implementation needs to be checked")
  }
  # refresh=0 turns off messages except errors from stan
  # quiet function silences stan output 
  fit <- model$sample(data=data_site_occ, seed=13, chains=1, 
                            iter_sampling=mcmc_iter, iter_warmup=100, refresh=0, show_messages=F)
  
  gen_occupancy <- fit$summary(variables=generated_vars[2])$mean
  v_est <- design_criteria(criteria="brier", sim_occ=gen_occupancy, obs_occ=occ)
  return(v_est)
}

estimate_v_nuthatch <- function(model, n_surveys, data_reps, m, sites, area_a, intensity_covars, bias_covars, 
                                select_idx, select_sites, r_survey_data, r_po_data,
                                p_logging, params, generated_vars, k_i, k_b, model_selection, mcmc_iter){
  
  estimate_mat <- matrix(nrow=data_reps, ncol=1)
  for (r in 1:data_reps){
    print("rth dataset:")
    print(r)
    # For each rth dataset get the m randomly chosen sites
    selected_data <- r_survey_data$Y[r, select_idx]
    #selected_occ <- r_survey_data$occupancy[r, select_idx]
    occ <- r_survey_data$occupancy[r, ]
    print(selected_data)
    print(selected_occ)
    # We use all po data for given r
    selected_po <- r_po_data$Y[r,]
    # Here I use the same X covariate for the presence-only data as used for the survey data. The difference is that
    # the survey data here is a subset (select_sites) of the sites, whereas the presence only data 
    # needs covariates for the whole grid to approximate the expected count in the entire region.
    if (model_selection==1){
      data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, area_a=area_a, X=select_sites, Y=selected_data, PO=selected_po,
                           X_po=intensity_covars, Z_po=bias_covars, k_i=k_i, k_b=k_b, model_diag=0)
    }else if (model_selection==2){
      data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, area_a=area_a, X=select_sites, Y=selected_data,
                           X_po=intensity_covars, k_i=k_i)
    }else{
      stop("No other models implemented")
    }
    # refresh=0 turns off messages except errors from stan
    # quiet function silences stan output 
    fit <- model$sample(data=data_site_occ, seed=13, chains=1, 
                        iter_sampling=mcmc_iter, iter_warmup=100, refresh=0, show_messages=F)
    if (p_logging==T){
      print("Logging posterior")
      posterior <- fit$draws()
      print(fit$summary(variables=params))
      
      color_scheme_set("mix-blue-pink")
      p_trace <- mcmc_trace(posterior,  pars = params,
                            facet_args = list(nrow = 2, labeller = label_parsed))
      print(p_trace + facet_text(size = 15))
      
      plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
      p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = params) + plot_title
      print(p_post)
      
      # Check issues with divergences
      color_scheme_set("darkgray")
      nuts_fit <- nuts_params(fit)
      diverge_p <- mcmc_parcoord(posterior, pars = params, np = nuts_fit, alpha=.1, 
                                 np_style=parcoord_style_np(div_alpha=1, div_size=.5))
      print(diverge_p)
    }
    gen_occupancy <- fit$summary(variables=generated_vars[1])$mean
    estimate_mat[r,1] <- design_criteria(criteria="brier", sim_occ=gen_occupancy, obs_occ=occ)
  }
}

estimate_v_parallel_nuthatch <- function(combined_df, model, n_surveys, m, sites, area_a, intensity_covars, bias_covars, 
                                         select_idx, select_sites, generated_vars, k_i, k_b, model_selection, mcmc_iter){
  tryCatch({
    # For each rth dataset get the m randomly chosen sites for the occupancy data, Y surveys, and PO.
    occ <- combined_df[1:sites]
    selected_data <- combined_df[sites+select_idx]
    PO_data <- combined_df[(2*sites+1):length(combined_df)] 
    # Here I use the same X covariate for the presence-only data as used for the survey data. The difference is that
    # the survey data here is a subset (select_sites) of the sites, whereas the presence only data 
    # needs covariates for the whole grid to approximate the expected count in the entire region.
    if (model_selection==1){
      data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, area_a=area_a, X=select_sites, Y=selected_data, PO=PO_data,
                           X_po=intensity_covars, Z_po=bias_covars, k_i=k_i, k_b=k_b, model_diag=0)
    }else if (model_selection==2){
      data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, area_a=area_a, X=select_sites, Y=selected_data,
                           X_po=intensity_covars, k_i=k_i)
    }else{
      stop("No other models implemented")
    }
    # refresh=0 turns off messages except errors from stan
    # quiet function silences stan output 
    fit <- "temp"
    v_est <- "temp"
    gen_occupancy <- "temp"
    fit <- model$sample(data=data_site_occ, seed=13, chains=1, 
                        iter_sampling=mcmc_iter, iter_warmup=100, refresh=0, show_messages=F)
    
    gen_occupancy <- fit$summary(variables=generated_vars[1])$mean
    v_est <- design_criteria(criteria="brier", sim_occ=gen_occupancy, obs_occ=occ)
    return(v_est)
  }, error=function(e){
    return(list("Issue computing V, returning_current_vars", n_surveys, selected_occ, selected_data, PO_data, 
                 model_selection, select_sites, fit, v_est, gen_occupancy))
  })
}

distance_func <- function(x, center_coord){
  dist <- sqrt((x[1] - center_coord[1])^2 + (x[2] - center_coord[2])^2)
  return(dist)
}

get_neighbors <- function(surface, l){
  # The rows from dist() are names as the indices for each row from surface, 
  # but these don't correspond to coordinates themselves
  distances <- as.matrix(dist(surface[,1:2]))
  # Data structure for nn; could be df too
  nearest_neighbors <- matrix(nrow=nrow(surface), ncol=l)
  for (row in 1:nrow(surface)){
    # Don't actually need the coords here; we're working just with row indices
    #print(paste(x_coord, y_coord))
    # This gives us the closest neighboring rows (in surface) for a set of coordinates
    min_points <- sort(distances[row, ])[2:(l+1)]
    stopifnot(length(min_points)==l)
    # Take column names, which gives us a character vector annoyingly
    nearest_neighbor_rows <- as.numeric(names(min_points))
    nearest_neighbors[row,] <- nearest_neighbor_rows
  }
  return(nearest_neighbors)
}

get_neighbors_sf_grid <- function(centroids, l){
  # st_distances gives us the nearest distance from one feature to another
  # So if two polygons are touching, distance will be 0 
  # So don't use the polygons, use the centroids of the polygons.
  distances <- st_distance(centroids, centroids)
  nearest_neighbors <- matrix(nrow=length(centroids), ncol=l)
  for (row in 1:length(centroids)){
    # This gives us the closest neighboring rows
    min_points <- order(distances[row, ])[2:(l+1)]
    stopifnot(!(row%in%min_points))
    stopifnot(length(min_points)==l)
    # Take column names, which gives us a character vector annoyingly
    nearest_neighbors[row,] <- min_points 
  }
  return(nearest_neighbors)
}

exchange_coordinates_deterministic <- function(select_id, cur_site_index, new_neighbor){
  # In the deterministic exchange, we pick a site, iterate through its neighbors
  # and then switch that site with n>0 to n=0 and the other site to n>0, for 
  # each neighbors
  select_id[cur_site_index] <- new_neighbor
  return(select_id)
}

exchange_coordinates_stochastic <- function(select_id, current_site, cur_site_index, n_neighbors, l){
  # In the exchange algorithm, we pick a site, find a neighbor close to that site,
  # and then switch that site with n>0 to n=0 and the other site to n>0, for 
  # whatever value of n we are using.
  # Make sure new neighbors are not in current list.
  new_site <- F
  while (new_site==F){
    # Get the neighbors for the selected site
    site_neighbors <- n_neighbors[current_site,]
    # Randomly sample 1 new neighbor
    new_neighbor_id <- sample(1:l, 1)
    new_neighbor <- site_neighbors[new_neighbor_id]
    # Make sure the new neighbor is not already in the list of current ids, to avoid selecting the same site twice
    if (all(new_neighbor!=select_id)){
      new_site <- T
    }
  }
  # Assertion to check that new neighbor is not in any of previous ids
  stopifnot(all(new_neighbor!=select_id))
  # Replace old site id with new site id
  select_id[cur_site_index] <- new_neighbor
  return(select_id)
}

specify_corr <- function(surface_data){
  distances <- as.matrix(dist(surface_data))
  corr <- exp(-distances/3)
  return(corr)
}

get_sampling_surface_simple <- function(k){
  # Define the window of interest
  dim <- c(k, k)
  win <- owin(c(0,dim[1]), c(0,dim[2]))
  
  # set number of pixels to simulate an environmental covariate
  spatstat.options(npixel=c(dim[1],dim[2]))
  
  y0 <- seq(win$yrange[1], win$yrange[2],
            length=spatstat.options()$npixel[2])
  x0 <- seq(win$xrange[1], win$xrange[2],
            length=spatstat.options()$npixel[1])
  multiplier <- 1/dim[2]
  # Get coordinate combinations to make a dataframe,
  # in the same format as other data gen functions
  surface_data <- expand.grid(x0, y0)
  names(surface_data) <- c("x", "y")
  
  aux_x <- outer(x0,y0, function (x,y) multiplier*y + 0*x)
  aux_z <- outer(x0,y0, function (x,y) 0*y + multiplier*x)
  surface_data$aux_x <- c(t(aux_x))
  surface_data$aux_z <- c(t(aux_z))
  return(surface_data)
}

get_sampling_surface_donut <- function(k, deviation){
  # good values for deviation are 2 or 5
  x <- seq(1:k)
  y <- seq(1:k)
  center_coord <- c(k/2, k/2)
  # Generate all possible coordinate points
  sampling_grid <- expand.grid(x, y)
  sampling_grid
  # Then compute distance from each location to center of the grid.
  # This gives a gaussain function that peaks a certain distance (7) away from the centroid
  r <- apply(sampling_grid, 1, distance_func, center_coord=center_coord)
  # Then create x covariate
  x <- exp(-10*((r-7)/deviation)^2)
  # This gives us the data for the auxiliary surveys x_i in the paper, which will be the 
  # covariates used in the paper X_i^T = [1, X_i1]
  x_1 <- qnorm(0.98*x + .01)
  sampling_grid$aux_x <- x_1
  colnames(sampling_grid) <- c("x", "y", "aux_x")
  return(sampling_grid)
}

get_bias_surface_correlated <- function(sampling_grid, aux_cor){
  # Various ways to generate correlated vector from one already existing:
  # https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
  V <- matrix(c(1, aux_cor, aux_cor, 1), nrow=2, ncol=2)
  R <- chol(V)
  aux_z <- rnorm(k*k, 0, 1)
  X <- cbind(sampling_grid$aux_x, aux_z)
  cor_X <- X %*% R 
  print("Correlation between X and Z")
  print(cor(cor_X))
  #cor_X[1:5, 1]
  sampling_grid$aux_z <- cor_X[,2]
  return(sampling_grid)
}

get_bias_surface_exponential <- function(sampling_grid, centroid){
  x <- seq(1:k)
  y <- seq(1:k)
  # Generate all possible coordinate points
  grid_points <- expand.grid(x, y)
  # Then compute distance from centroid to every other location
  r <- apply(grid_points, 1, function(x, center_coord){sqrt((center_coord[1] - x[1])^2 + (center_coord[2] - x[2])^2)}, center_coord=centroid)
  # Then create z covariate
  x <- exp(-2*((r)/5))
  x_1 <- qnorm(0.98*x + .01)
  sampling_grid$aux_z <- x_1
  return(sampling_grid)
}

generate_data_so <- function(data_reps, surface_data, corr_matrix, p_0, alpha, beta, sigma, n, sites, link){
  # Generate some random covariate data that will be used to model theta
  occupancy_maps <- matrix(nrow=data_reps, ncol=sites)
  Y_detection <- matrix(nrow=data_reps, ncol=sites)
  theta_reps <- matrix(nrow=data_reps, ncol=sites)
  for (r in 1:data_reps){
    # Equivalent way to induce correlation  
    #R <- t(chol(corr_matrix)) 
    #theta <- b_0 + sampling_surface$aux*b_1 + R %*% rnorm(sites)
    if (link=="cloglog"){
      g_theta <- 1 - exp(-exp(alpha + surface_data$aux_x*beta))
    }else{
      theta <- rnorm(sites, alpha + surface_data$aux_x*beta, sigma*corr_matrix)
      g_theta <- pnorm(theta)
    }
    site_presence <- rbinom(sites, 1, g_theta)
    theta_reps[r, ] <- g_theta
    occupancy_maps[r, ] <- site_presence
    Y_detection[r, ] <- rbinom(sites, n, p_0*site_presence)
  } 
  return(list(occupancy=occupancy_maps, Y=Y_detection, theta=theta_reps))
}

save_basic_plots <- function(fig_name, p){
  ggsave(fig_name, plot=p, dpi=300, width=7, height=6, units="cm")
}

plot_sites <- function(sampling_surface, r_po_data, select_idx, title, fig_text_size){
  p <- ggplot() + 
    geom_tile(sampling_surface, mapping=aes(x, y, fill=aux_x)) + 
    geom_point(data=r_po_data$Y_coords, mapping=aes(x=x, y=y), size=2, colour="orange") +
    geom_point(data=sampling_surface[select_idx,], aes(x=x, y=y), colour = "white", size = 1.5) +
    scale_fill_viridis(discrete=FALSE, name="x1") +
    labs(title=title, x="", y="") + 
    theme(text=element_text(size=fig_text_size), axis.title = element_text(size = 5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed() 
  return(p)
}

plot_sites_vs_best <- function(sampling_surface, current_site, select_idx, best_select_idx, prev_visits, optimal_visits, title, fig_text_size){
  # Plot current set of sites compared to the best sites. 
  # Best will always be pink
  prev_size <- prev_visits/(length(prev_visits)) + 1.5
  best_size <- optimal_visits/(length(optimal_visits)) + 1
  p <- ggplot(sampling_surface, aes(x, y, fill=aux_x)) + 
    geom_tile() +
    geom_point(data=sampling_surface[select_idx,], aes(x=x, y=y), colour = "white", size = prev_size) +
    geom_point(data=sampling_surface[best_select_idx,], aes(x=x, y=y), colour = "hotpink1", size = best_size, alpha=1) +
    geom_point(data=sampling_surface[current_site,], aes(x=x, y=y), colour = "black", size = 0.5) +
    scale_fill_viridis(discrete=FALSE, name="x1") +
    labs(title=title, x="", y="") + 
    theme(text=element_text(size=fig_text_size), axis.title = element_text(size = 5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed() 
  return(p)
}

plot_po_optimal_sites <- function(sampling_surface, r_po_data, best_select_idx, optimal_visits, title, fig_text_size){
  best_size <- optimal_visits/(length(optimal_visits)) + 1
  p <- ggplot() +
    geom_tile(sampling_surface, mapping=aes(x, y, fill=aux_x)) + 
    geom_point(data=r_po_data$Y_coords, mapping=aes(x=x, y=y), size=2, colour="orange") +
    geom_point(data=sampling_surface[best_select_idx,], aes(x=x, y=y), colour = "white", size = best_size) +
    scale_fill_viridis(discrete=FALSE, name="x1") +
    labs(title=title, x="", y="") + 
    theme(text=element_text(size=fig_text_size), axis.title = element_text(size = 5), legend.key.size = unit(0.25, 'cm')) +
    coord_fixed() 
  return(p)
}

# Another set of functions for plotting design over a map
plot_sites_ebird <- function(site_centroids, rast_surface, site_idx, title){
  p <- ggplot() + 
    geom_spatraster(data=rast_surface) +
    geom_sf(data = site_centroids[site_idx], color=alpha("white",1)) + 
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, na.value="white") +
    theme_minimal()+
    ggtitle(title) +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm'))
  return(p)
}

plot_sites_vs_best_ebird <- function(site_centroids, rast_surface, current_site, select_idx, best_select_idx, prev_visits, optimal_visits, title){
  # Plot current set of sites compared to the best sites. 
  # Best will always be pink
  prev_size <- prev_visits/(length(prev_visits)) + 1.5
  best_size <- optimal_visits/(length(optimal_visits)) + 1
  # Site centroids are 1 dimensional
  p <- ggplot() + 
    geom_spatraster(data=rast_surface) +
    geom_sf(data = site_centroids[select_idx], color=alpha("white",1), size=prev_size) + 
    geom_sf(data = site_centroids[best_select_idx], color=alpha("hotpink",1), size=best_size) + 
    geom_sf(data = site_centroids[current_site], color=alpha("black",1), size=0.5) +
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, na.value="grey") +
    theme_minimal()+
    ggtitle(title) +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm'))
  return(p)
}

plot_po_optimal_sites_ebird <- function(site_centroids, rast_surface, po_data, best_select_idx, optimal_visits, title){
  best_size <- optimal_visits/(length(optimal_visits)) + 1
  p <- ggplot() + 
    geom_spatraster(data=rast_surface) +
    geom_sf(data=po_data, color=alpha("orange",0.5), size=0.5)+
    geom_sf(data = site_centroids[best_select_idx], color=alpha("hotpink",1), size=best_size) + 
    scale_fill_viridis_c(begin=0.2, end=1, option="viridis",alpha=0.7, na.value="grey") +
    theme_minimal()+
    ggtitle(title) +
    theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm'))
  return(p)
}

write_results <- function(random_starts, best_v, best_site_mat, visits, v_df, exp_dir, exp_name){
  # This will be useful if I have multiple criteria
  v_stat_df <- data.frame(v=sum(best_v) / random_starts, v_var=var(best_v))
  best_v <- data.frame(best_v=best_v)
  site_df <- data.frame(t(best_site_mat), t(visits))
  names(site_df) <- rep(paste("rand-start", c(1:random_starts), sep=""), 2)
  names(v_df) <- c("iter", "v", "rand-start")
  xlsx_list <- list("v_stat"=v_stat_df, "best_v"=best_v, "v_iterations"=v_df, "best_sites"=site_df)
  #write.xlsx(xlsx_list, file=file.path(exp_dir, paste("results_", exp_name, ".xlsx", sep="")), rowNames=F)
  write.xlsx(xlsx_list, file=file.path(exp_dir, "results.xlsx"), rowNames=F)
}
