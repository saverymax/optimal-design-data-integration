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
  for (site in 1:length(counts)){
    if (counts[site] > 0){
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


generate_ppp_data_r <- function(surface_data, params, n_sites, data_reps, corr_matrix, gp_bool, area_D){
  # Generate some random covariate data per site on a grid that will be used to model lambda and b
  # There will be r (data_reps) datasets replicated, for use with the exchange 
  # algorithm/monte carlo integration and to maintain efficiency in the parallel 
  # function computation of V
  lambda <- (area_D / n_sites) * exp(params$alpha + params$beta*surface_data$aux_x)
  eta <- exp(params$gamma + params$delta*surface_data$aux_z)
  b <- eta / (1 + eta)
  if (gp_bool==T){
    gp <- exp(rnorm(n_sites, 0, corr_matrix))
    Y_counts <- rpois(n_sites, lambda * b * gp)
  }
  else{
    Y_counts <- rpois(n_sites, lambda * b)
  }
  # Generate coordinate points within each site for each observation.
  Y_coords <- generate_coords(surface_data, Y_counts)
  # Create matrix of data reps, though for po data we just copy everything,
  # as opposed to random generation for each rep.
  Y_po <- matrix(rep(Y_counts, data_reps), nrow=data_reps, ncol=n_sites, byrow=T)
  lambdas <- matrix(rep(lambda, data_reps), nrow=data_reps, ncol=n_sites, byrow=T)
  biases <- matrix(rep(b, data_reps), nrow=data_reps, ncol=n_sites, byrow=T)
  return(list(Y=Y_po, lambda=lambdas, bias=biases, Y_coords=Y_coords))
}

