#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/evaluation/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/evaluation/
#PBS -N exp_5_misspec_analysis
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=20gb


module load CmdStanR
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/fit_integrated_so_models.R
