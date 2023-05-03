#!/bin/bash
#PBS -N test_job
#PBS -o /data/leuven/459/vsc459/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/leuven/459/vsc459/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=1

module load R
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
Rscript $WORKDIR/experiments/test_script.R
##R $WORKDIR/test_script.R > $WORKDIR/job_output/out_${PBS_JOBID}.txt
