rm run_initial_model_comparison.sh
models="1 3"
#alpha="-2 -1.5 -1 -0.5"
alpha="-2"
beta="2"
gamma="-2 -1.5 -1 -0.5 0"
delta="0.5"
detection_p="0.2 0.7"
# Might be nice to write job output to the experimental run dir but it's nice to leave that dir created by the R script
# so as to seperate the HPC and local run capabilities.
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
cores=48
intensity="donut"
for m in $models
do
for g in $gamma
do
for p in $detection_p
do
exp_name=optimal_design_model_${m}_gamma_${g}_p_${p}
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N $exp_name
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --data_reps=$(($cores*2)) --m=5 --model_selection=$m --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_fun=\"$intensity\" --v_parallel --cores=$cores --alpha=$alpha --beta=$beta --gamma=$g --delta=$delta --p=$p" > $exp_name.sh
echo "qsub $exp_name.sh" >> run_initial_model_comparison.sh
#echo "qsub optimal_design_model_${m}_gamma_${g}.sh" >> run_initial_model_comparison.sh
#echo "qsub optimal_design_model_${m}_alpha_${a}.sh" >> run_initial_model_comparison.sh

done
done
done
