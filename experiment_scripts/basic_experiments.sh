#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N optimal_design_test_run
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb

##module load R
##Rscript $WORKDIR/experiments/test_script.R
module load CmdStanR
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --data_reps=100 --m=5 --model_selection=3 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_fun="simple" --v_parallel
