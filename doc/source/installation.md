# Installation

## Dependencies

The R package versions used in this work are listed below. All results reported in the paper were generated from code using these software versions. The script for the applied work requires a few extra spatial packages.

- cmdstanr: 0.5.2
- openxlsx: 4.2.5
- optparse:1.7.1       
- spatstat:2.3-4
- reshape2: 1.4.4        
- tidyverse: 1.3.1
- tidyr: 1.2.0
- dplyr: 1.0.9           
- bayesplot: 1.9.0       
- ggplot2: 3.3.6
- viridis: 0.6.2         
- kableExtra: 1.3.4

For application: 

- tidyterra: 0.4.0       
- terra: 1.7-55          
- sf: 1.0-14

We plan to create a Rocker (Docker for R) container for this application in the near future.

## CmdStan

This code requires CmdStanR and CmdStan. To install these, a few extra steps beyond the typical R pacakge installation command are required. Please make sure you have these packages available on your system following the instructions here: <https://mc-stan.org/cmdstanr/articles/cmdstanr.html>. Note that it is necessary to install CmdStan separately after installing CmdStanR.