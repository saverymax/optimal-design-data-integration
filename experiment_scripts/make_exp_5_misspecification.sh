exp_dir=exp_5
rm -rf exp_5
mkdir $exp_dir
rm run_$exp_dir.sh
models="3"
alpha="-2"
beta="0.5"
#beta="0.5 1"
gamma="1"
delta="0.25 2"
misspec="0 0.01 0.05 0.1"
p="0.2"
n_surveys="5"
total_sites="5"
deviation="5"
# TODO: Add both run paths and gamma integration
run_path=("oracle" "sequential")
# Might be nice to write job output to the experimental run dir but it's nice to leave that dir created by the R script
# so as to seperate the HPC and local run capabilities.
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
cores=48
intensity="donut"
bias="exponential"
for b in $beta
do 
for d in $delta
do
for m in $models
do
for sites in $total_sites
do
for mc in $misspec
do
for r in "${run_path[@]}"
do
if [[ "$r" == "sequential" ]]; then
# Separate files for integration and no integration
int_status=no_int
exp_name=oe_model-${m}_sites-${sites}_n-${n_surveys}_peak-${deviation}_a-${alpha}_b-${b}_g-${gamma}_d-${d}_p-${p}_misspec-${mc}_"$r"_$int_status
echo $exp_name
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N ${exp_dir}_${exp_name}
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR/0.5.2-foss-2022a-R-4.2.1
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --save_dir=experimental_runs/$exp_dir --exp_name=$exp_name --data_reps=$(($cores*2)) --m=$sites --min_visits=1 --max_visits=$n_surveys --vary_visits --model_selection=$m --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file="po_gen_ints-${intensity}_peak=${deviation}_a=${alpha}_b=${b}_g=${gamma}_d=${d}_e=${mc}.Rds" --intensity_func=\"$intensity\" --sd=${deviation} --bias_func=\"$bias\" --v_parallel --cores=$cores --alpha=${alpha} --beta=${b} --gamma=${gamma} --delta=${d} --p=$p --misspec_run=\"$r\" --misspec=${mc} --posterior_file=pp_posterior_po_gen_ints-donut_peak=5_a=${alpha}_b=${b}_g=${gamma}_d=${d}_e=${mc}_$int_status.Rds" > $exp_dir/$exp_name.sh
echo "qsub $exp_dir/$exp_name.sh" >> run_$exp_dir.sh

int_status=int
exp_name=oe_model-${m}_sites-${sites}_n-${n_surveys}_peak-${deviation}_a-${alpha}_b-${b}_g-${gamma}_d-${d}_p-${p}_misspec-${mc}_"$r"_$int_status
echo $exp_name
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N ${exp_dir}_${exp_name}
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR/0.5.2-foss-2022a-R-4.2.1
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --save_dir=experimental_runs/$exp_dir --exp_name=$exp_name --data_reps=$(($cores*2)) --m=$sites --min_visits=1 --max_visits=$n_surveys --vary_visits --model_selection=$m --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file="po_gen_ints-${intensity}_peak=${deviation}_a=${alpha}_b=${b}_g=${gamma}_d=${d}_e=${mc}.Rds" --intensity_func=\"$intensity\" --sd=${deviation} --bias_func=\"$bias\" --v_parallel --cores=$cores --alpha=${alpha} --beta=${b} --gamma=${gamma} --delta=${d} --p=$p --misspec_run=\"$r\" --misspec=${mc} --posterior_file=pp_posterior_po_gen_ints-donut_peak=5_a=${alpha}_b=${b}_g=${gamma}_d=${d}_e=${mc}_$int_status.Rds" > $exp_dir/$exp_name.sh
echo "qsub $exp_dir/$exp_name.sh" >> run_$exp_dir.sh

else
exp_name=oe_model-${m}_sites-${sites}_n-${n_surveys}_peak-${deviation}_a-${alpha}_b-${b}_g-${gamma}_d-${d}_p-${p}_misspec-${mc}_"$r"
echo $exp_name
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N ${exp_dir}_${exp_name}
#PBS -l walltime=10:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR/0.5.2-foss-2022a-R-4.2.1
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --save_dir=experimental_runs/$exp_dir --exp_name=$exp_name --data_reps=$(($cores*2)) --m=$sites --min_visits=1 --max_visits=$n_surveys --vary_visits --model_selection=$m --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file="po_gen_ints-${intensity}_peak=${deviation}_a=${alpha}_b=${b}_g=${gamma}_d=${d}_e=${mc}.Rds" --intensity_func=\"$intensity\" --sd=${deviation} --bias_func=\"$bias\" --v_parallel --cores=$cores --alpha=${alpha} --beta=${b} --gamma=${gamma} --delta=${d} --p=$p --misspec_run=\"$r\" --misspec=${mc}" > $exp_dir/$exp_name.sh
echo "qsub $exp_dir/$exp_name.sh" >> run_$exp_dir.sh

fi
done
done
done
done
done
done
