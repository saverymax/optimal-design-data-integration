# Functions for generating and analysing presence-only data.
# Relies on sampling surface data generated in "experimental_design_functions.R" script.

generate_coords <- function(grid, counts){
  # Function to generate coordinates for each observation based on the 
  # count generated from the poisson variable.
  pos_count <- sum(counts)
  cont_coords <- matrix(nrow=pos_count, ncol=2)
  x_coords <- c()
  y_coords <- c()
  k <- .5
  for (site in 1:nrow(counts)){
    if (counts[site,1] > 0){
      # Could make the coordinates dependent on the covariates, which would be interesting.
      # Currently assumes quadrates are centered on integers.
      x_coord <- runif(counts[site], grid$x[site]-k, grid$x[site]+k)
      y_coord <- runif(counts[site], grid$y[site]-k, grid$y[site]+k)
      x_coords <- c(x_coords, x_coord)
      y_coords <- c(y_coords, y_coord)
    } 
  }
  cont_coords[,1] <- x_coords
  cont_coords[,2] <- y_coords
  cont_coords <- as.data.frame(cont_coords)
  names(cont_coords) <- c("x", "y")
  return(cont_coords)
}


generate_ppp_data <- function(surface_data, params, n_sites){
  # Generate some random covariate data per site on a grid that will be used to model lambda and b
  Y_counts <- matrix(nrow=n_sites, ncol=1)
  lambda <- exp(params$alpha + params$beta*surface_data$aux_x)
  eta <- exp(params$gamma + params$delta*surface_data$aux_z)
  b <- eta / (1 + eta)
  # Account for true occupancy maps
  Y_counts[, 1] <- rpois(n_sites, lambda * b)
  # At the same time, generate coordinate points within each site for each observation.
  Y_coords <- generate_coords(surface_data, Y_counts)
  
  p <- ggplot(surface_data, aes(x, y, fill=lambda)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Generated intensity per site")
  print(p)
  
  p <- ggplot(surface_data, aes(x, y, fill=b)) + 
    geom_tile() +
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Generated bias per site")
  print(p)
  
  p <- ggplot() +
    geom_tile(surface_data, mapping=aes(x, y, fill=lambda*b, width=1, height=1), alpha=.6) + 
    scale_fill_viridis(discrete=FALSE) +
    ggtitle("Generated thinning per site, including generated (or observed) individuals") +
    geom_point(data=Y_coords, mapping=aes(x=x, y=y), size=3, col="white") +
    theme(panel.grid.minor = element_line(colour="white")) +
    scale_y_continuous(breaks = seq(0, 20, 1)) +
    scale_x_continuous(breaks = seq(0, 20, 1))
  print(p)
  
  return(list(Y=Y_counts, coordinates=Y_coords))
}

generate_ppp_data_r <- function(surface_data, params, n_sites, data_reps){
  # Generate some random covariate data per site on a grid that will be used to model lambda and b
  # These data will be generated based on "true" or generated occupancy maps, 
  # so that in cases where it is known that there are no species occupying the site,
  # there will be no counts generated there.
  # There will be r (data_reps) datasets generated, for use with the exchange algorithm/monte carlo integration
  # This also takes into account "true" occupancy maps
  Y_po <- matrix(nrow=data_reps, ncol=sites)
  lambdas <- matrix(nrow=data_reps, ncol=sites)
  biases <- matrix(nrow=data_reps, ncol=sites)
  # In the occupancy maps, the reps are stored per row, with each site in each column.
  for (i in 1:data_reps){
    Y_counts <- matrix(nrow=n_sites, ncol=1)
    lambda <- exp(params$alpha + params$beta*surface_data$aux_x)
    eta <- exp(params$gamma + params$delta*surface_data$aux_z)
    b <- eta / (1 + eta)
    Y_counts[, 1] <- rpois(n_sites, lambda * b)
    # Not including occupancy maps anymore
    #Y_counts[, 1] <- rpois(n_sites, lambda * b * occ_maps[i, ])
    # Generate coordinate points within each site for each observation.
    # These aren't saved in this function, except the last iteration for viz purposes
    Y_coords <- generate_coords(surface_data, Y_counts)
    Y_po[i, ] <- Y_counts
    lambdas[i, ] <- lambda
    biases[i, ] <- b
  }
  return(list(Y=Y_po, lambda=lambdas, bias=biases, Y_coords=Y_coords))
}
