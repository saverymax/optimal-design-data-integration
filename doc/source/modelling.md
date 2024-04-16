# Modelling

## Application

The two supported models in the application module of the code are shown below. They are written in the Bayesian programming language Stan. The models used in this work are site occupancy models with presence-absence (PA) data only and with PA + presence-only (PO) data as a prior on the intensity parameters.  

The first model shown below does not use the PO data. We can see that there is a simple wide normal prior placed on the parameters (alpha and betas). We do not model bias (delta or gamma) here because it is not necessary when we are not using the PO data.

```

nuthatch_site_occ_no_po <- '
    data{
      real<lower = 0> area_a;
      int<lower = 1> n_pa_sites;
      int<lower = 1> n_po_sites;
      int<lower=0> k_i;
      array[n_pa_sites] int Y;
      array[n_pa_sites] int n_surveys;
      matrix[n_pa_sites, k_i] X;
      matrix[n_po_sites, k_i] X_po; // predictor matrix for intensity
    }
    parameters{
      real alpha;           // intercepts
      vector[k_i] beta;       // vector of params for intensity
      real<lower = 0, upper = 1> p;
    }
    model{
      // priors
      target += normal_lpdf(alpha |   0,10);
      target += normal_lpdf(beta | 0,10);
      vector[n_pa_sites] g_theta;
      g_theta = 1 - exp(-exp(alpha + X*beta + log(area_a)));
      for (i in 1:n_pa_sites) {
        if (Y[i] > 0){
          target += log(g_theta[i]*choose(n_surveys[i], Y[i])*(p^Y[i])*(1-p)^(n_surveys[i]-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys[i]) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      // Required generated quantities
      // We generate over all sites, instead of just those being used for PA
      vector[n_po_sites] g_theta_gen;
      array[n_po_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha + X_po*beta + log(area_a)));
      occ_gen = bernoulli_rng(g_theta_gen);
    }
'
```
The mext model uses presence-only data as a prior. This approach uses the Poisson approximation to the non-homogeneous Poisson process (NHPP) -- as discussed in the paper -- to model the number of PO observations at each site in the region. This Poisson distributed variable is parameterized via the intensity and bias, but the occupancy probability is modelled only as a function of intensity. Wide normal priors are placed on the parameters for the NHPP. 
```

nuthatch_poisson_process_site_occupancy <- '
    data{
      real<lower = 0> area_a;
      int<lower = 1> n_pa_sites;
      int<lower = 1> n_po_sites;
      int<lower=0> k_i;
      int<lower=0> k_b;
      int<lower=0, upper=1> model_diag;
      array[n_pa_sites] int Y;
      array[n_pa_sites] int n_surveys;
      matrix[n_pa_sites, k_i] X;
      array[n_po_sites] int PO;
      matrix[n_po_sites, k_i] X_po; // predictor matrix for intensity
      matrix[n_po_sites, k_b] Z_po; // predictor matrix for bias
    }
    parameters{
      real alpha;           // intercepts
      real gamma;
      vector[k_i] beta;       // vector of params for intensity
      vector[k_b] delta;      // vector of params for bias
      real<lower = 0, upper = 1> p;
    }
    model{
      // priors
      target += normal_lpdf(alpha |   0,10);
      target += normal_lpdf(gamma | 0,10);
      target += normal_lpdf(beta | 0,10);
      target += normal_lpdf(delta | 0,10);
      vector[n_pa_sites] g_theta;
      target += poisson_log_lpmf(PO | log(area_a) + alpha + X_po*beta + gamma + Z_po*delta);
      g_theta = 1 - exp(-exp(alpha + X*beta + log(area_a)));
      for (i in 1:n_pa_sites) {
        if (Y[i] > 0){
          target += log(g_theta[i]*choose(n_surveys[i], Y[i])*(p^Y[i])*(1-p)^(n_surveys[i]-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys[i]) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      vector[n_po_sites] g_theta_gen;
      array[n_po_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha + X_po*beta + log(area_a)));
      occ_gen = bernoulli_rng(g_theta_gen);
    }
'

```

Generally, the Stan command used to fit these during during the exchange algorithm is 
```
fit <- model$sample(data=data_site_occ, seed=13, chains=1, iter_sampling=mcmc_iter, iter_warmup=100, refresh=0, show_messages=F)
```
This can be found in the ```estimate_v_nuthatch``` function in the ```experimental_design_functions.R``` module within the source code. Of course, the variables used to fit this model are in memory during the exchange algorithm. The model is fit R times during the exchange algorithm, and we do not provide instructions here to fit a one-off model. This can be provided upon request.

When running the exchange algorithm these models can be selected using the options 1 (with PO + PA data) and 2 (without PO data, only PA data). For example
```
Rscript optimal_design_application.R --working_dir=. --save_dir=experimental_runs/test_application --data_save_dir=data/ebird --exp_name=test_application --data_reps=4 --m=10 --min_visits=1 --max_visits=4 --vary_visits --random_starts=10 --model_selection=2 --p=0.2 --po_sample_prop=0.05 --exch_iter=40 --mcmc_iter=1000 --v_parallel --cores=4 --run_time=2
```
where ```--model_selection=2``` will specify the use of only PA data within the algorithm.


## Simulated setting

The simulated setting uses the same models as in the applied setting. However, there is only one covariate used to model the intensity in the simulated setting, so that the Stan models change slightly. These are shown below.
```
cloglog_site_occupancy <- '
  data{
      int<lower = 1> n_sites;
      int<lower = 1> total_sites;
      // Distinction in the model for total sites and selected ones
      vector[n_sites] X;
      vector[total_sites] X_all;
      array[n_sites] int Y;
      array[n_sites] int n_surveys;
    }
    parameters{
      real alpha;
      real beta;
      real<lower = 0, upper = 1> p;
    }
    model{
      // priors
      target += normal_lpdf(alpha | 0,10);
      target += normal_lpdf(beta | 0,10);
      vector[n_sites] g_theta;
      g_theta = 1 - exp(-exp(alpha + beta * X));
      for (i in 1:n_sites) {
        if (Y[i] > 0){
          target += log(g_theta[i]*choose(n_surveys[i], Y[i])*(p^Y[i])*(1-p)^(n_surveys[i]-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys[i]) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      vector[total_sites] g_theta_gen;
      array[total_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha + beta * X_all));
      // Posterior predictive distribution
      occ_gen = bernoulli_rng(g_theta_gen);
  }
'

poisson_process_site_occupancy <- '
  data{
      int<lower = 1> n_pa_sites;
      int<lower = 1> n_po_sites;
      int<lower=0, upper=1> model_diag;
      array[n_pa_sites] int Y;
      array[n_pa_sites] int n_surveys;
      vector[n_pa_sites] X;
      vector[n_po_sites] X_po;
      vector[n_po_sites] Z_po;
      array[n_po_sites] int PO;
    }
    parameters{
      real alpha;
      real beta;
      real gamma;
      real delta;
      real<lower = 0, upper = 1> p;
    }
    model{
      // priors
      target += normal_lpdf(alpha |   0,10);
      target += normal_lpdf(beta | 0,10);
      target += normal_lpdf(gamma | 0,10);
      target += normal_lpdf(delta | 0,10);
      // Call prior function for presence-only here.
      vector[n_pa_sites] g_theta;
      target += poisson_log_lpmf(PO | alpha + beta * X_po + gamma + delta * Z_po);
      g_theta = 1 - exp(-exp(alpha + beta * X));
      for (i in 1:n_pa_sites) {
        if (Y[i] > 0){
          target += log(g_theta[i]*choose(n_surveys[i], Y[i])*(p^Y[i])*(1-p)^(n_surveys[i]-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys[i]) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      vector[n_po_sites] g_theta_gen;
      array[n_po_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha + beta * X_po));
      occ_gen = bernoulli_rng(g_theta_gen);
  }
'
```
To select these models in the optimal design algorithm use options 1 (no PO data) and 3 (PO data). For example, 
```
Rscript optimal_design_site_occ.R --working_dir=. --save_dir=experimental_runs/test_run --exp_name=design-test-run --data_reps=4 --m=5 --min_visits=1 --max_visits=4 --vary_visits --model_selection=3 --random_starts=3 --exch_iter=10 --mcmc_iter=1000 --use_sim_po --po_data_file=po_gen_peak_ints-donut_peak=2_a=-2_b=0.5_g=1_d=0.25.Rds --intensity_func="donut" --sd=2 --bias_func="exponential" --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --p=0.2
```
where ```--model_selection=3``` will specify the use of PO data in the site occupancy model used in the exchange algorithm. It is also possible to specify models 2 (constant occupancy probability with no PO data) and 4 (constant occupancy probability with PO data), but these parameterize the intensity and bias with only an intercept and are for debugging purposes only.


