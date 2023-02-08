
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

