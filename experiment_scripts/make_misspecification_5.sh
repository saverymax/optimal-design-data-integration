exp_dir=exp_5
rm -rf exp_5
mkdir $exp_dir
rm run_$exp_dir.sh
models="1 3"
alpha="-2"
beta="0.5"
gamma="1"
delta="0.25"
misspec="0 0.01 0.05 0.1"
p="0.2"
n_surveys="5"
total_sites="5"
deviation="5"
run_path="sequential"
# Might be nice to write job output to the experimental run dir but it's nice to leave that dir created by the R script
# so as to seperate the HPC and local run capabilities.
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
cores=48
intensity="donut"
bias="exponential"
for m in $models
do
for sites in $total_sites
do
for mc in $misspec
do
for r in $run_path
do
exp_name=oe_model-${m}_sites-${sites}_n-${n_surveys}_peak-${deviation}_a-${alpha}_b-${beta}_g-${gamma}_d-${delta}_p-${p}_misspec-${mc}_${run_path}
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N ${exp_dir}_${exp_name}
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --save_dir=experimental_runs/$exp_dir --exp_name=$exp_name --data_reps=$(($cores*2)) --m=$sites --min_visits=1 --max_visits=$n_surveys --vary_visits --model_selection=$m --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file="po_gen_peak_ints-${intensity}_peak=${deviation}_a=${alpha}_b=${beta}_g=${gamma}_d=${delta}_e=${mc}.Rds" --intensity_func=\"$intensity\" --sd=${deviation} --bias_func=\"$bias\" --v_parallel --cores=$cores --alpha=${alpha} --beta=${beta} --gamma=${gamma} --delta=${delta} --p=$p --misspec_run=\"$run_path\" --misspec=${mc} --gamma_integration --posterior_file=pp_posterior_po_gen_ints-donut_peak=5_a=${alpha}_b=${beta}_g=${gamma}_d=${delta}_e=${mc}.Rds" > $exp_dir/$exp_name.sh
echo "qsub $exp_dir/$exp_name.sh" >> run_$exp_dir.sh
done
done
done
done
