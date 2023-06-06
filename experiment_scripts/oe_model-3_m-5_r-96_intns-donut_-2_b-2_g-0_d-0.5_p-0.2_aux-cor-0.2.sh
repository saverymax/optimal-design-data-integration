#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N oe_model-3_m-5_r-96_intns-donut_-2_b-2_g-0_d-0.5_p-0.2_aux-cor-0.2
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=100gb

module load CmdStanR
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/optimal_design_site_occ.R --working_dir=/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration --exp_name=oe_model-3_m-5_r-96_intns-donut_-2_b-2_g-0_d-0.5_p-0.2_aux-cor-0.2 --data_reps=96 --m=5 --model_selection=3 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_fun="donut" --v_parallel --cores=48 --alpha=-2 --beta=2 --gamma=0 --delta=0.5 --p=0.2 --aux_cor=0.2
