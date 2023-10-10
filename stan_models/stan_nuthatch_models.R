nuthatch_poisson_process <- '
    data{
      real<lower = 0> area_a;
      int<lower = 1> N;
      int<lower=0> k_i; // number of predictors for intensity
      int<lower=0> k_b; // number of predictors for bias
      matrix[N, k_i] X; // predictor matrix for intensity
      matrix[N, k_b] Z; // predictor matrix for bias
      array[N] int y;
    }
    parameters{
      real alpha;
      real gamma;
      vector[k_i] beta;       // vector of params for intensity
      vector[k_b] delta;      // vector of params for bias
    }
    model{
      //priors
      target += normal_lpdf(alpha | 0,10);
      target += normal_lpdf(beta | 0,10);
      target += normal_lpdf(gamma | 0,10);
      target += normal_lpdf(delta | 0,10);

      // likelihood
      // log parameterization so we dont have to exponentiate coefs before passing into poisson distribution
      // https://mc-stan.org/docs/functions-reference/poisson-distribution-log-parameterization.html
      target += poisson_log_lpmf(y |  log(area_a) * (alpha + X * beta + gamma + Z * delta));
    }
    generated quantities{
      //vector[N] lambda_rep;
      //vector[N] lambda_rep_log;
      //vector[N] b_rep;
      //vector[N] lambda_bias;
      array[N] int y_rep;
      //lambda_rep = exp(alpha + beta * X);
      //lambda_rep_log = alpha + beta * X;
      ////.*abs(alpha + beta * X);
      //b_rep = exp(gamma + delta * Z);
      ////.*abs(gamma + delta * Z);
      //// Thin the process
      ////lambda_bias = lambda_rep .* b_rep;
      //lambda_bias = exp(alpha + beta * X + gamma + delta * Z);
      //// User guide has some examples of log poisson generation and poisson PPCs
      y_rep = poisson_log_rng(log(area_a) + alpha + X * beta + gamma + Z * delta);
    }
'

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
      // Required generated quantities
      // We can generate over all sites, instead of just those being used for PA
      vector[n_po_sites] g_theta_gen;
      array[n_po_sites] int occ_gen; 
      // Note that we do not include gamma here because the 
      // probability of occupancy is just P(N(C_i)>0) = G(alpha + beta X)
      g_theta_gen = 1 - exp(-exp(alpha + X_po*beta + log(area_a)));
      // We dont need to generate the occupancy if we take gtheta
      // since we take the mean of the occupancy anyway which is the mean of gtheta
      occ_gen = bernoulli_rng(g_theta_gen);
    }
'