Code for work with optimal design of experiments, using data integration in bayesian site occupancy models.

Currently, this code takes a site-occupancy model with a complementary log-log likelihood, and uses this to model presence-absence data, while using as a prior on the parameters a poisson process model for presence-only data.

The experimental results will be written to a folder named experimental_runs, which you have to create yourself. However, this will include only plots from the exchange and figures of the intial data generation. The rest of the output printed during the experiment is written to std-out, which is also up to you how to control. Typically I write these to properly named .o and .e files when running on the HPC cluster.
