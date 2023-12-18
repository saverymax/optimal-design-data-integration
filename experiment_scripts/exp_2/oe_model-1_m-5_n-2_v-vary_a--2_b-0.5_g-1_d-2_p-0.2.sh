#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N oe_model-1_m-5_n-2_v-vary_a--2_b-0.5_g-1_d-2_p-0.2
#PBS -l walltime=6:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=100gb

module load CmdStanR
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/optimal_design_site_occ.R --working_dir=/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration --exp_name=oe_model-1_m-5_n-2_v-vary_a--2_b-0.5_g-1_d-2_p-0.2 --data_reps=96 --m=5 --min_visits=1 --max_visits=2 --vary_visits --model_selection=1 --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file=po_gen_ints-donut_a=-2_b=0.5_g=1_d=2.Rds --intensity_func="donut" --bias_func="exponential" --v_parallel --cores=48 --alpha=-2 --beta=0.5 --gamma=1 --delta=2 --p=0.2
