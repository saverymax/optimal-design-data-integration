# Little test script with some of my data and parallel library in R.
library(parallel)
n_cores <- detectCores()
clust <- makeCluster(n_cores)

df <- data.frame(c(1,2,3), c(1,2,3))




#combined_df <- cbind(r_survey_data$occupancy, r_survey_data$Y, r_po_data$Y)

#estimate_mat <- parApply(clust, combined_df, 1, FUN=estimate_v_parallel, model, n_surveys, data_reps, m, sites, sampling_surface, 
