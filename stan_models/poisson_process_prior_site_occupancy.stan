
  data{
      int<lower = 1> n_surveys;
      int<lower = 1> n_sites;
      int<lower = 1> total_sites;
      // Distinction in the model for total sites and selected ones
      vector[n_sites] X;
      vector[total_sites] X_all;
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
      //vector[n_sites] lambda_rep;
      //array[n_sites] int y_rep; 
      //lambda_rep = exp(alpha + beta * X);
      // Is poisson the right choice here?
      //y_rep = poisson_log_rng(alpha + beta * X);
      vector[total_sites] g_theta_gen;
      array[total_sites] int occ_gen; 
      g_theta_gen = 1 - exp(-exp(alpha + beta * X_all));
      // Posterior predictive distribution
      occ_gen = bernoulli_rng(g_theta_gen);
  }

