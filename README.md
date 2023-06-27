This is the repository for work on optimal design of sampling locations for biodiversity monitoring applications, using data integration in Bayesian site occupancy models.

Currently, this code implements an exchange algorithm to search for optimal spatial locations in a region of interest. 
To do so, it takes a site-occupancy model with a complementary log-log likelihood, and uses this to model presence-absence data, while using as a prior on the parameters a poisson process model for presence-only data.
Integrating data in this fashion allows us to incorporate more information than would otherwise be gained by using presence-absence data alone.

Further information, including instructions regarding how to run the code and more in-depth implementation details can be found 
at the Sphinx documentation: https://saverymax.github.io/optimal-design-data-integration/ 

