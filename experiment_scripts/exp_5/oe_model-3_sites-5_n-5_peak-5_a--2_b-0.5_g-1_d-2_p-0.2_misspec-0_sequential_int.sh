#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N exp_5_oe_model-3_sites-5_n-5_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0_sequential_int
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=100gb

module load CmdStanR/0.5.2-foss-2022a-R-4.2.1
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/optimal_design_site_occ.R --working_dir=/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration --save_dir=experimental_runs/exp_5 --exp_name=oe_model-3_sites-5_n-5_peak-5_a--2_b-0.5_g-1_d-2_p-0.2_misspec-0_sequential_int --data_reps=96 --m=5 --min_visits=1 --max_visits=5 --vary_visits --model_selection=3 --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file=po_gen_ints-donut_peak=5_a=-2_b=0.5_g=1_d=2_e=0.Rds --intensity_func="donut" --sd=5 --bias_func="exponential" --v_parallel --cores=48 --alpha=-2 --beta=0.5 --gamma=1 --delta=2 --p=0.2 --misspec_run="sequential" --misspec=0 --posterior_file=pp_posterior_po_gen_ints-donut_peak=5_a=-2_b=0.5_g=1_d=2_e=0_int.Rds
