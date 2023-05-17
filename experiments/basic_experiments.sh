#!/bin/bash
#PBS -N test_job
#PBS -o /data/leuven/459/vsc459/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/leuven/459/vsc459/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=20gb


module load CmdStanR
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
##Rscript $WORKDIR/optimal_design_site_occ.R --data_reps=100 --m=5 --model_selection=3 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_fun="simple" --v_parallel
Rscript $WORKDIR/experiments/test_script.R
##R $WORKDIR/test_script.R > $WORKDIR/job_output/out_${PBS_JOBID}.txt
