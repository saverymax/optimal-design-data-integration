site_occupany_detection <- '
  data{
      int<lower = 1> n_surveys;
      int<lower = 1> n_sites;
      vector[n_sites] X;
      array[n_sites] int Y;
    }
    parameters{
      real beta_0;
      real beta_1;
      real<lower = 0, upper = 1> p;
    }
    model{
      // priors
      target += beta_lpdf(p | 1,1);
      target += normal_lpdf(beta_0 | 0,10);
      target += normal_lpdf(beta_1 | 0,10);

      // likelihood
      vector[n_sites] g_theta;
      // DO NOT USE LPDF! Because we do log later
      // For when using model without site specific effects
      g_theta = Phi(beta_0 + beta_1*X);
      for (i in 1:n_sites) {
        if (Y[i] > 0){
          target += log(g_theta[i]*choose(n_surveys, Y[i])*(p^Y[i])*(1-p)^(n_surveys-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      vector[n_sites] theta_rep;
      array[n_sites] int y_rep; 
      // Posterior predictive distributions
      theta_rep = Phi(beta_0 + beta_1*X);
      y_rep = binomial_rng(n_surveys, theta_rep);
  }
'

site_occupany_detection_site_specific <- '
  data{
      int<lower = 1> n_surveys;
      int<lower = 1> n_sites;
      vector[n_sites] X;
      array[n_sites] int Y;
      matrix[n_sites, n_sites] cor_mat;
    }
    parameters{
      real beta_0;
      real beta_1;
      real<lower = 0, upper = 1> p;
      vector[n_sites] theta;
    }
    model{
      // priors
      target += beta_lpdf(p | 1,1);
      target += normal_lpdf(beta_0 | 0,10);
      target += normal_lpdf(beta_1 | 0,10);
      theta ~ multi_normal(beta_0 + beta_1*X, cor_mat);
      // likelihood
      vector[n_sites] g_theta;
      g_theta = Phi(theta);
      for (i in 1:n_sites) {
        // take into account site specific effects
        // See https://mc-stan.org/docs/2_19/stan-users-guide/multivariate-hierarchical-priors-section.html
        // Not a bad idea to parameterize in terms of cholesky
        // Also, chunking the site specific effects would be more interesting and efficient.
        if (Y[i] > 0){
          target += log(g_theta[i]*choose(n_surveys, Y[i])*(p^Y[i])*(1-p)^(n_surveys-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      vector[n_sites] theta_rep;
      vector[n_sites] phi_rep;
      array[n_sites] int y_rep; 
      // Posterior predictive distributions
      theta_rep = multi_normal_rng(beta_0 + beta_1*X, cor_mat);
      phi_rep = Phi(beta_0 + beta_1*X);
      y_rep = binomial_rng(n_surveys, phi_rep);
  }
'
cloglog_site_occupancy <- '
  data{
      int<lower = 1> n_surveys;
      int<lower = 1> n_sites;
      vector[n_sites] X;
      array[n_sites] int Y;
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
          target += log(g_theta[i]*choose(n_surveys, Y[i])*(p^Y[i])*(1-p)^(n_surveys-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      vector[n_sites] lambda_rep;
      array[n_sites] int y_rep; 
      lambda_rep = exp(alpha + beta * X);
      // Posterior predictive distribution
      // Is poisson the right choice here?
      y_rep = poisson_rng(lambda_rep);
  }
'

poisson_process_site_occupancy <- '
  data{
      int<lower = 1> n_surveys;
      int<lower = 1> n_pa_sites;
      int<lower = 1> n_po_sites;
      int<lower=0, upper=1> model_diag;
      vector[n_pa_sites] X;
      array[n_pa_sites] int Y;
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
          target += log(g_theta[i]*choose(n_surveys, Y[i])*(p^Y[i])*(1-p)^(n_surveys-Y[i]));
        }
        else{
          // Compute mixture of no detection and no occupancy
          target += log(g_theta[i]*(1-p)^(n_surveys) + (1 - g_theta[i]));
        }
      }
    }
    generated quantities{
      // Dont need these generated quantities at the moment
      //if (model_diag){
      //vector[n_po_sites] lambda_rep;
      //vector[n_po_sites] b_rep;
      //vector[n_po_sites] lambda_bias_rep;
      //array[n_po_sites] int y_rep; 
      //lambda_rep = exp(alpha + beta * X_po);
      //b_rep = exp(gamma + delta * Z_po);
      //lambda_bias_rep = exp(alpha + beta * X_po + gamma + delta * Z_po);
      //y_rep = poisson_log_rng(alpha + beta * X_po + gamma + delta * Z_po);
      //}
      // Required generated quantities
      // We can generate over all sites, instead of just those being used for PA
      vector[n_po_sites] g_theta_gen;
      array[n_po_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha + beta * X_po));
      // Posterior predictive distribution for occupancy
      occ_gen = bernoulli_rng(g_theta_gen);
  }
'
