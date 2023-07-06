#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N oe_model-3_m-5_n-1_r-96_intns-donut_a--2_b-0.5_g-1_d-0.25_p-0.2
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=100gb

module load CmdStanR
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/optimal_design_site_occ.R --working_dir=/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration --exp_name=oe_model-3_m-5_n-1_r-96_intns-donut_a--2_b-0.5_g-1_d-0.25_p-0.2 --data_reps=96 --m=5 --min_visits=0 --max_visits=1 --model_selection=3 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_func="donut" --bias_func="exponential" --v_parallel --cores=48 --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --p=0.2 --aux_cor=0.8
