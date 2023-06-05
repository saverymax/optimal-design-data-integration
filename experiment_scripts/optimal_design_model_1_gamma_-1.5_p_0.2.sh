#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N optimal_design_model_1_gamma_-1.5_p_0.2
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=100gb

module load CmdStanR
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/optimal_design_site_occ.R --working_dir=/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration --data_reps=96 --m=5 --model_selection=1 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_fun="donut" --v_parallel --cores=48 --alpha=-2 --beta=2 --gamma=-1.5 --delta=0.5 --p=0.2
