exp_dir=exp_4
rm -rf exp_4
mkdir $exp_dir
rm run_$exp_dir.sh
models="1 3"
alpha="-2"
beta="0.5"
# Permute the bias params
gamma="1"
delta="2"
p="0.2"
n_surveys="2 5"
total_sites="2 5 7 10"
# Might be nice to write job output to the experimental run dir but it's nice to leave that dir created by the R script
# so as to seperate the HPC and local run capabilities.
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
cores=48
intensity="donut"
bias="exponential"
for m in $models
do
for g in $gamma
do
for d in $delta
do
for n in $n_surveys
do
for sites in $total_sites
do
exp_name=oe_model-${m}_sites-${sites}_n-${n}_a-${alpha}_b-${beta}_g-${g}_d-${d}_p-${p}
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N $exp_name
#PBS -l walltime=3:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --exp_name=$exp_name --data_reps=$(($cores*2)) --m=$sites --min_visits=1 --max_visits=$n --vary_visits --model_selection=$m --random_starts=10 --exch_iter=40 --mcmc_iter=1000 --use_sim_po --po_data_file="po_gen_ints-${intensity}_a=${alpha}_b=${beta}_g=${g}_d=${d}.Rds" --intensity_func=\"$intensity\" --bias_func=\"$bias\" --v_parallel --cores=$cores --alpha=$alpha --beta=$beta --gamma=$g --delta=$d --p=$p" > $exp_dir/$exp_name.sh
echo "qsub $exp_dir/$exp_name.sh" >> run_$exp_dir.sh
done
done
done
done
done
