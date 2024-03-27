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
      // Here we need to have a likelihood*prior function for the presence-only data, where I can plug in the value
      // of the current parameter (alpha/beta/gamma/delta), and get the value of 
      // the posterior every time the value of alpha/beta/delta/gamma change.
      // These are the priors for the parameters specified in the model for presence-only data
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
      // Here we need to have a likelihood*prior function for the presence-only data, where I can plug in the value
      // of the current parameter (alpha/beta/gamma/delta), and get the value of 
      // the posterior every time the value of alpha/beta/delta/gamma change.
      // These are the priors for the parameters specified in the model for presence-only data
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
      // Required generated quantities
      // We can generate over all sites, instead of just those being used for PA
      vector[n_po_sites] g_theta_gen;
      array[n_po_sites] int occ_gen; 
      // Note that we do not include gamma here because the 
      // probability of occupancy is just P(N(C_i)>0) = G(alpha + beta X)
      g_theta_gen = 1 - exp(-exp(alpha + beta * X_po));
      // Posterior predictive distribution for occupancy
      occ_gen = bernoulli_rng(g_theta_gen);
  }
'

# Using constant intensity parameterized by only an intercept with no PO data
pp_site_occ_constant_no_po <- '
  data{
      int<lower = 1> n_pa_sites;
      int<lower = 1> total_sites;
      array[n_pa_sites] int n_surveys;
      array[n_pa_sites] int Y;
    }
    parameters{
      real alpha;
      real<lower = 0, upper = 1> p;
    }
    model{
      real g_theta;
      // priors
      target += normal_lpdf(alpha |   0,10);
      g_theta = 1 - exp(-exp(alpha));
      for (i in 1:n_pa_sites) {
        if (Y[i] > 0){
          target += log(g_theta*choose(n_surveys[i], Y[i])*(p^Y[i])*(1-p)^(n_surveys[i]-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta*(1-p)^(n_surveys[i]) + (1 - g_theta));
        }
      }
    }
    generated quantities{
      // We can generate over all sites, instead of just those being used for PA
      real g_theta_gen;
      vector[total_sites] g_theta_vec;
      array[total_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha));
      g_theta_vec = rep_vector(g_theta_gen, total_sites);
      // Posterior predictive distribution for occupancy
      occ_gen = bernoulli_rng(g_theta_vec);
  }
'

# Using constant intensity parameterized by only an intercept with prior PO
pp_site_occ_constant_po <- '
  data{
      int<lower = 1> n_pa_sites;
      int<lower = 1> n_po_sites;
      array[n_pa_sites] int n_surveys;
      array[n_pa_sites] int Y;
      array[n_po_sites] int PO;
    }
    parameters{
      real alpha;
      real<lower = 0, upper = 1> p;
    }
    model{
      real g_theta;
      // priors
      target += normal_lpdf(alpha |   0,10);
      target += poisson_log_lpmf(PO | alpha);
      g_theta = 1 - exp(-exp(alpha));
      for (i in 1:n_pa_sites) {
        if (Y[i] > 0){
          target += log(g_theta*choose(n_surveys[i], Y[i])*(p^Y[i])*(1-p)^(n_surveys[i]-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta*(1-p)^(n_surveys[i]) + (1 - g_theta));
        }
      }
    }
    generated quantities{
      // We can generate over all sites, instead of just those being used for PA
      real g_theta_gen;
      vector[n_po_sites] g_theta_vec;
      array[n_po_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha));
      g_theta_vec = rep_vector(g_theta_gen, n_po_sites);
      // Posterior predictive distribution for occupancy
      occ_gen = bernoulli_rng(g_theta_vec);
  }
'