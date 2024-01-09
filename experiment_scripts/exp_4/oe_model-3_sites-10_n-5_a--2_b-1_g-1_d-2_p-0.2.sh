#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N oe_model-3_sites-10_n-5_a--2_b-1_g-1_d-2_p-0.2
#PBS -l walltime=13:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=100gb

module load CmdStanR
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/optimal_design_site_occ.R --working_dir=/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration --save_dir=experimental_runs/exp_4 --exp_name=oe_model-3_sites-10_n-5_a--2_b-1_g-1_d-2_p-0.2 --data_reps=96 --m=10 --min_visits=1 --max_visits=5 --vary_visits --model_selection=3 --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file=po_gen_peak_ints-donut_a=-2_b=1_g=1_d=2.Rds --intensity_func="donut" --bias_func="exponential" --v_parallel --cores=48 --alpha=-2 --beta=1 --gamma=1 --delta=2 --p=0.2
