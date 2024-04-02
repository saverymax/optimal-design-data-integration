# Script to pre-generate PO data.
# Need to set intensity params here and use the same ones in the optimal design script

library(optparse)
library(dplyr)
library(ggplot2)
library(viridis)
library(hrbrthemes)

parser <- OptionParser()
parser <- add_option(parser, "--working_dir", type="character", default=".", help="Path to the directory containing script")
parser <- add_option(parser, "--exp_name", type="character", default="po_gen", help="Base name to save data")
parser <- add_option(parser, "--alpha", type="double", default=-2, help="Intercept for intensity")
parser <- add_option(parser, "--beta", type="double", default=0.5, help="Slope for intensity")
parser <- add_option(parser, "--gamma", type="double", default=1, help="Intercept for bias")
parser <- add_option(parser, "--delta", type="double", default=.5, help="Slope for bias")
parser <- add_option(parser, "--intensity_func", type="character", default="donut", help="Intensity function for sampling surface")
parser <- add_option(parser, "--sd", type="double", default=5, help="Standard deviation for donut intensity surface")
parser <- add_option(parser, "--bias_func", type="character", default="exponential", help="Bias function for sampling surface")
parser <- add_option(parser, "--area", type="integer", default=100, help="Area of region D")
parser <- add_option(parser, "--k", type="integer", default=20, help="Number of sites along one side of grid")
parser <- add_option(parser, "--data_reps", type="integer", default=4, help="Number of dataset reps for criterion estimation")
exp_args <- parse_args(parser)
print(exp_args)

source(file.path(exp_args$working_dir, "experimental_design_functions.R"))
source(file.path(exp_args$working_dir, "presence_only_functions.R"))

data_dir <- file.path(exp_args$working_dir, "data", "sim_data")
dir.create(data_dir)

# Area for whole space, which allows us to set area for sites based on number of sites.
area_D <- exp_args$area
# For generating data according to GP
gp_bool <- F
data_reps <- exp_args$data_reps
# k is one side of grid
k <- exp_args$k
sites <- k^2
# alpha=2 indicates "good" quality of auxiliary information
alpha <- exp_args$alpha
beta <- exp_args$beta
gamma <- exp_args$gamma
delta <- exp_args$delta
# standard deviation of surface
deviation <- exp_args$sd
# Select the intensity surface
if (exp_args$intensity_func == "simple"){
  sampling_surface <- get_sampling_surface_simple(k)
}else{
  sampling_surface <- get_sampling_surface_donut(k, deviation)
  if (exp_args$bias_func == "exponential"){
    centroid <- c(10,4)
    sampling_surface <- get_bias_surface_exponential(sampling_surface, centroid)
  }else if (exp_args$bias_func == "correlation"){
    sampling_surface <- get_bias_surface_correlated(sampling_surface, aux_cor)
  }
  else{
    stop("Bias function not implemented")
  }
  # Filter for only a quarter of the grid.
  # If we filter, we need to change the total sites as well
  sites <- sites/4
  sampling_surface <- sampling_surface %>% dplyr::filter(x<11, y<11)
}
print(sampling_surface)
corr_matrix <- specify_corr(sampling_surface[,1:2])

params <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta)
r_po_data <- generate_ppp_data_r(sampling_surface, params, sites, data_reps, corr_matrix, gp_bool, area_D)

r_po_data$params <- list(alpha=alpha, beta=beta, gamma=gamma, delta=delta)

# Save the PO data sets
param_setting <- paste(exp_args$exp_name, "_", "ints-", exp_args$intensity_func, "_peak=", deviation ,"_a=", alpha, "_b=", beta, "_g=", gamma, "_d=", delta, sep="")
p <- ggplot(sampling_surface, aes(x, y, fill=r_po_data$lambda[data_reps,])) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, name="Lambda") +
  ggtitle("Generated intensity per site") +
  theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
  coord_fixed()
fig_name <- file.path(data_dir, paste("intensity-site_", param_setting, ".png", sep=""))
save_basic_plots(fig_name, p)

p <- ggplot(sampling_surface, aes(x, y, fill=r_po_data$bias[data_reps,])) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE, "Bias") +
  ggtitle("Generated bias per site") +
  theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
  coord_fixed()
fig_name <- file.path(data_dir, paste("bias-site_", param_setting,".png", sep=""))
save_basic_plots(fig_name, p)

thinned_intensity <- r_po_data$lambda[data_reps,]*r_po_data$bias[data_reps,]
p <- ggplot() +
  geom_tile(sampling_surface, mapping=aes(x, y, fill=thinned_intensity, width=1, height=1), alpha=.6) + 
  scale_fill_viridis(discrete=FALSE, name="L*b") +
  ggtitle("Generated PP thinning per site, \nincluding observed PO individuals") +
  geom_point(data=r_po_data$Y_coords, mapping=aes(x=x, y=y), size=2, col="orange") +
  theme(panel.grid.minor = element_line(colour="white")) +
  scale_y_continuous(breaks = seq(0, 20, 1)) +
  scale_x_continuous(breaks = seq(0, 20, 1)) +
  theme(text=element_text(size=5), legend.key.size = unit(0.25, 'cm')) +
  coord_fixed()
fig_name <- file.path(data_dir, paste("po-thinning-per-site_", param_setting, ".png", sep=""))
save_basic_plots(fig_name, p)

saveRDS(r_po_data, file=file.path(data_dir, paste(param_setting, ".Rds", sep="")))
