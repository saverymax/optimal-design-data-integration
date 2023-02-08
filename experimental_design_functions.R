# Functions for experimental design

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

brier_score <- function(z, z_obs){
  b <- 0
  sites <- ncol(z)
  iter <- nrow(z)
  for (i in 1:sites)
    z_bar_i <- sum(z[, i]) / iter
    if(z_obs[i] == 0){
      b <- b + (z_bar_i)^2 
    }else{
      b <- b + (1-z_bar_i)^2
    }
  return(b/sites)
}

brier_score_stan <- function(z, z_obs){
  # Brier score more compatible with stan output
  b <- 0
  sites <- length(z)
  for (i in 1:sites)
    z_bar_i <- z[i]
    if(z_obs[i] == 0){
      b <- b + (z_bar_i)^2 
    }else{
      b <- b + (1-z_bar_i)^2
    }
  return(b/sites)
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
                       select_idx, select_sites, selected_data, r_survey_data, r_po_data, 
                       p_logging, params){
  for (r in 1:data_reps){
    # For each rth dataset get the m randomly chosen sites
    selected_data <- r_survey_data$Y[r, select_idx]
    selected_occ <- r_survey_data$occupancy[r, select_idx]
    # Here I use the same X covariate for the presence-only data as used for the survey data. The difference is that
    # the survey data here is a subset (select_sites) of the sites, whereas the presence only data 
    # needs covariates for the whole grid to approximate the expected count in the entire region.
    data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, X=select_sites$aux_x, Y=selected_data, PO=r_po_data$Y[r,], 
                         X_po=sampling_surface$aux_x, Z_po=sampling_surface$aux_z)
    #print(data_site_occ)
    # refresh=0 turns off messages except errors from stan
    # quiet function silences stan output 
    # TODO: Check out issues with divergences...
    fit <- quiet(model$sample(data=data_site_occ, seed=13, chains=1, 
                              iter_sampling=1000, iter_warmup=100, refresh=0, show_messages=F))
    if (p_logging==T){
      print("Logging posterior")
      posterior <- fit$draws()
      print(fit$summary(variables=params))
      
      color_scheme_set("mix-blue-pink")
      p_trace <- mcmc_trace(posterior,  pars = params,
                            facet_args = list(nrow = 2, labeller = label_parsed))
      print(p_trace + facet_text(size = 15))
      
      plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
      p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = c("alpha", "beta")) + plot_title
      print(p_post)
      
      plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
      p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", regex_pars = c("gamma", "delta")) + plot_title
      print(p_post)
      
      plot_title <- ggtitle(paste("Posterior distribution of detection probability, mean and 90% interval"))
      p_post <- mcmc_areas(posterior,  prob = 0.9, point_est="mean", pars = c("p")) + plot_title
      print(p_post)
      
      # TODO: Plot site specific probs
      #ppd_count_df <- data.frame(x=sampling_surface$x, y=sampling_surface$y, 
      #                           occ=fit$summary(variables=generated_vars[1])$mean)
      
      #p <- ggplot(ppd_count_df, aes(x, y, fill=occ)) + 
      #  geom_tile() +
      #  scale_fill_viridis(discrete=FALSE) +
      #  ggtitle("Occupancy probability per site")
      #print(p)
      
    }
    # Average estimate after R iterations through the datasets.
    gen_occupancy <- fit$summary(variables=c(generated_vars))
    # Take the mean of the generated quantity for the posterior estimate
    # TODO: Check if the magnitude of V makes sense
    estimate_mat[r,1] <- design_criteria(criteria="brier", sim_occ=gen_occupancy$mean, obs_occ=selected_occ)
  }
}


distance_func <- function(x, center_coord){
  dist <- sqrt((x[1] - center_coord)^2 + (x[2] - center_coord)^2)
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

exchange_coordinates_deterministic <- function(select_id, current_site, cur_site_index, n_neighbors, l){
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

exchange_coordinates <- function(select_id, current_site, cur_site_index, n_neighbors, l){
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

get_sampling_surface <- function(k){
  # Function generates X and correlated Z auxiliary data
  x <- seq(1:k)
  y <- seq(1:k)
  center_coord <- k/2
  # Generate all possible coordinate points
  sampling_grid <- expand.grid(x, y)
  sampling_grid
  # Then compute distance from each location to center of the grid.
  r <- apply(sampling_grid, 1, distance_func, center_coord=center_coord)
  r
  # Then create x covariate
  x <- exp(-10*((r-7)/5)^2)
  # This gives us the data for the auxiliary surveys x_i in the paper, which will be the 
  # covariates used in the paper X_i^T = [1, X_i1]
  x_1 <- qnorm(0.98*x + .01)
  sampling_grid$aux_x <- x_1
  
  # This gives us auxiliary data for the intensity parameter lambda.
  # Next generate some correlated data for the bias parameter b.
  # Various ways to generate correlated vector from one already existing:
  # https://stats.stackexchange.com/questions/15011/generate-a-random-variable-with-a-defined-correlation-to-an-existing-variables
  V <- matrix(c(1, .8, .8, 1), nrow=2, ncol=2)
  R <- chol(V)
  aux_z <- rnorm(k*k, 0, 1)
  X <- cbind(sampling_grid$aux_x, aux_z)
  cor_X <- X %*% R 
  #sampling_grid$aux_x[1:5]
  print("Correlation between X and Z")
  print(cor(cor_X))
  #cor_X[1:5, 1]
  sampling_grid$aux_z <- cor_X[,2]
  
  # Now generate the plot of auxiliary data for the grid
  colnames(sampling_grid) <- c("x", "y", "aux_x", "aux_z")
  # Heatmap 
  p <- ggplot(sampling_grid, aes(x, y, fill=aux_x)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Initial sampling surface, X covariate")
  print(p)
  
  p <- ggplot(sampling_grid, aes(x, y, fill=aux_z)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Initial sampling surface, Z covariate")
  print(p)
  
  
  return(sampling_grid)
}

generate_data <- function(data_reps, surface_data, corr_matrix, p_0, b_0, b_1, sigma, n, sites){
  # Generate some random covariate data that will be used to model theta
  occupancy_maps <- matrix(nrow=data_reps, ncol=sites)
  Y_detection <- matrix(nrow=data_reps, ncol=sites)
  theta_reps <- matrix(nrow=data_reps, ncol=sites)
  for (r in 1:data_reps){
    # There was an error in the earlier code
    # This models spatial dependence in random effects for each site, which is why we use rnorm
    # In the strict probit model, we can just use X'B directly, ie pnorm(XB)
    
    # Equivalent way to induce correlation
    #R <- t(chol(corr_matrix)) 
    #theta <- b_0 + sampling_surface$aux*b_1 + R %*% rnorm(sites)
    theta <- rnorm(sites, b_0 + surface_data$aux_x*b_1, sigma*corr_matrix)
    g_theta <- pnorm(theta)
    site_presence <- rbinom(sites, 1, g_theta)
    theta_reps[r, ] <- theta
    occupancy_maps[r, ] <- site_presence
    Y_detection[r, ] <- rbinom(sites, n, p_0*site_presence)
  } 
  return(list(occupancy=occupancy_maps, Y=Y_detection, theta=theta_reps))
}

plot_sites <- function(sampling_surface, select_idx, title){
  p <- ggplot(sampling_surface, aes(x, y, fill=aux_x)) + 
    geom_tile() +
    geom_point(data=sampling_surface[select_idx,], aes(x=x, y=y), colour = "white", size = 3) +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle(title)
  print(p)
  return(p)
}

plot_sites_vs_best <- function(sampling_surface, select_idx, best_select_idx, title){
  # Plot current set of sites compared to the best sites. 
  p <- ggplot(sampling_surface, aes(x, y, fill=aux_x)) + 
    geom_tile() +
    geom_point(data=sampling_surface[select_idx,], aes(x=x, y=y), colour = "white", size = 3) +
    geom_point(data=sampling_surface[best_select_idx,], aes(x=x, y=y), colour = "hotpink1", size = 2, alpha=1) +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle(title)
  print(p)
  return(p)
}

