library(parallel)

n_cores <- detectCores()
print(n_cores)
clust <- makeCluster(4)
par_func <- function(df){
  a <- lapply(1:1e6, function(x) x^2)
}
start_time <- Sys.time()
a <- par_func(1)
end_time <- Sys.time()
cur_run_time <- end_time - start_time
print(cur_run_time)
df <- data.frame(c(1:10), c(1:10), c(1:10), c(1:10))
  #df <- data.frame(c(1:10), c(1:10))
start_time <- Sys.time()
estimate_mat <- parApply(clust, df, 2, FUN=par_func)

end_time <- Sys.time()
cur_run_time <- end_time - start_time
print(cur_run_time)
