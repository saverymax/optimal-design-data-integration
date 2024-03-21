# Modelling

## Simulated setting

There are 3 models that can be used within the exchange algorithm. All models are based on the site-occupancy likelihood, with varying priors. 

The most basic model treats the intensity as constant across the space, i.e., uses no covariates to model the intensity. This is likely an incorrect specification and is included solely for comparison purposes. 

The next model allows the intensity to vary over the space and is specified by the parameters alpha and beta. Normal or vague priors can be placed over the alpha and beta parameters.

The third model uses the posterior from presence-only data modelled as a non-homogenous poisson process as the prior for the parameters in the site-occupancy model. 

The code also supports site-occupancy models using the probit likelihood, both with and without site-specific effects. However, these have not been tested in the exchange and results are not reported.

## Application
