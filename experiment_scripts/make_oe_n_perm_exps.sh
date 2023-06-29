rm run_initial_model_comparison.sh
models="1 3"
#alpha="-2 -1.5 -1 -0.5"
alpha="-2"
beta="0.5"
gamma="-1"
#gamma="-2 -1 -0.5 1"
#delta="-0.5"
delta="-0.25 0.25 2"
#detection_p="0.2 0.7"
detection_p="0.2"
n_surveys="5 10"
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
for p in $detection_p
do
for n in $n_surveys
do
exp_name=oe_model-${m}_m-5_n-${n}_r-96_intns-${intensity}_a-${alpha}_b-${beta}_g-${g}_d-${d}_p-${p}
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N $exp_name
#PBS -l walltime=02:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --exp_name=$exp_name --data_reps=$(($cores*2)) --m=5 --n=$n --model_selection=$m --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_func=\"$intensity\" --bias_func=\"$bias\" --v_parallel --cores=$cores --alpha=$alpha --beta=$beta --gamma=$g --delta=$d --p=$p" --aux_cor=0.8 > $exp_name.sh
echo "qsub $exp_name.sh" >> run_initial_model_comparison.sh
done
done
done
done
done
