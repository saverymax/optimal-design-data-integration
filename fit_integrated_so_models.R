# Short script to fit the integrated SO model to the final designs generated
# from the design optimization algorithm.

library(ggplot2)
library(dplyr)
library(viridis)
library(cmdstanr) 
library(bayesplot)
library(optparse)
library(openxlsx)


parser <- add_option(parser, "--design_selection", type="character", default="", help="Name of experiment name that corresponds to file from experiment runs to load and use to fit the SO model")
parser <- add_option(parser, "--po_data_file", type="character", default="po_gen_ints-donut_a=-2_b=0.5_g=1_d=0.5.Rds", help="File name of PO data to use in optimization")
parser <- add_option(parser, "--model_selection", type="integer", default=3, help="Occupancy model to use.")

if (exp_args$model_selction!=3){
  stop("Only integrated model (3) is supported in this script")
}

mcmc_iter <- 2000
n_chains <- 3
data_reps <- 20

alpha <- -2

design_results <- read.xlsx(file_path, sheet="best_sites")
# There will be r designs for the r random starts.

pa_eval_data <- generate_data_so(data_reps, sampling_surface, corr_matrix, p_0, alpha, beta, sigma, visits[2], sites, link=link_func)

model_strings <- list(
	"cloglog_site_occupancy"=cloglog_site_occupancy, 
	"poisson_process_constant"=pp_site_occ_constant_no_po,
	"poisson_process_prior"=poisson_process_site_occupancy, 
	"poisson_process_constant_po_prior"=pp_site_occ_constant_po
	)
model_selection <- exp_args$model_selection
model_name <- names(model_strings)[model_selection]
print(paste("Using model", model_name, model_selection))
model_path <- file.path(stan_dir, paste(model_name, ".stan", sep=""))
print(paste("Cmdstan model", model_strings[[model_selection]]))
write(model_strings[[model_selection]], model_path)
model <- cmdstan_model(model_path) 

data_site_occ = list(n_surveys=n_surveys, n_pa_sites=m, n_po_sites=sites, X=select_sites$aux_x, Y=selected_data, PO=selected_po, 
                           X_po=sampling_surface$aux_x, Z_po=sampling_surface$aux_z, model_diag=0)
fit_so_1 <- model_so$sample(data=data_site_occ, seed=13, chains=n_chains, iter_sampling=mcmc_iter, iter_warmup=500)

all_params <- c("alpha", "beta[1]", "gamma", "delta[1]")
params_intercept <- c("alpha", "gamma")
fit_so_1$summary(all_params)
posterior <- fit_so_1$draws(all_params)

color_scheme_set("mix-blue-pink")
p_trace <- mcmc_trace(posterior,
                      facet_args = list(nrow = 2, labeller = label_parsed))
print(p_trace + facet_text(size = 15))

plot_title <- ggtitle(paste("Posterior distributions, with means and 90% interval"))
for (aux_param in all_params){
  mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars=c(aux_param)) + plot_title
}

plot_title <- ggtitle(paste("Posterior distributions of detection probability, mean and 90% interval"))
p_post <- mcmc_areas(posterior,  prob = 0.8, point_est="mean", pars = params_intercept) + plot_title
print(p_post)

mcmc_intervals(fit_so_1$draws(), pars=all_params)
mcmc_hist(fit_so_1$draws(), pars = all_params)
mcmc_pairs(fit_so_1$draws(), pars=all_params)
mcmc_scatter(fit_so_1$draws(), pars=c('alpha', 'gamma'))

# Examine the ppd for y
generated_occ <- fit_so_1$draws("occ_gen", format="matrix")
generated_prob <- fit_so_1$draws("g_theta_gen", format="matrix")
occ_means <- colMeans(generated_occ)
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
occ_mode <- apply(generated_occ, MARGIN=2, FUN="get_mode")
prob_means <- colMeans(generated_prob)
ppd_df <- data.frame(x=sampling_surface$x, y=sampling_surface$y, occ_mode=occ_mode, occ_rep=occ_means, theta_rep=prob_means)
head(ppd_df)

# Plot generated ys and thetas on the sampling grid
p <- ggplot(ppd_df, aes(x, y, fill=occ_rep)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("PPD occupancy per site")
print(p)

# theta_rep is already in prob space
p <- ggplot(ppd_df, aes(x, y, fill=theta_rep)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("PPD occupancy probability per site")
print(p)

p <- ggplot(ppd_df, aes(x, y, fill=occ_mode)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  ggtitle("PPD occupancy mode per site")
print(p)
