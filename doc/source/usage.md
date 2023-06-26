Basic usage of this code is described here.

The parameters used in the code are as follows:
```
"--working_dir", type="character", help="Path to the directory containing code to source for the main script"
"--exp_name", type="character", help="Name of current experiment, which is used for dir to save output"
"--data_reps", type="integer", default=10, help="Number of dataset reps for criterion estimation"
"--m", type="integer", default=5, help="Number of sites to survey"
"--n", type="integer", default=5, help="Number of time to visit each site"
"--model_selection", type="integer", default=3, help="Occupancy model to use"
"--random_starts", type="integer", default=3, help="Number of random starts to run the exchange"
"--exch_iter", type="integer", default=20, help="Number of iterations of exchange before ending optimization"
"--mcmc_iter", type="integer", default=1000, help="Number of MCMC iterations in Stan"
"--intensity_func", type="character", default="simple", help="Intensity function for sampling surface"
"--bias_func", type="character", default="exponential", help="Bias function for sampling surface"
"--p_logging", action="store_true", default=F, help="Boolean for logging information about posterior estimates"
"--v_parallel", action="store_true", default=F, help="Boolean for parallel computation of V criterion"
"--cores", type="integer", default=4, help="Number of cores to use for parallel processing"
"--alpha", type="double", default=-2, help="Intercept for intensity"
"--beta", type="double", default=2, help="Slope for intensity"
"--gamma", type="double", default=-1, help="Intercept for bias"
"--delta", type="double", default=.5, help="Slope for bias"
"--p", type="double", default=0.7, help="Probability of detection"
"--area", type="integer", default=100, help="Area of region D"
"--k", type="integer", default=20, help="Number of sites along one side of grid"
"--aux_cor", type="double", default=0.8, help="Correlation between auxiliary covariates"
```

