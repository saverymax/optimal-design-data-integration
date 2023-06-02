rm run_initial_model_comparison.sh
models="1 3"
alpha="-2 -1.5 -1 -0.5"
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
cores=48
intensity="donut"
for m in $models
do
for a in $alpha
do
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N optimal_design_model_$m
#PBS -l walltime=01:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR
Rscript $WORKDIR/optimal_design_site_occ.R --working_dir=$WORKDIR --data_reps=$(($cores*2)) --m=5 --model_selection=$m --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --intensity_fun=\"$intensity\" --v_parallel --cores=$cores --alpha=$a --beta=2 --gamma=-1 --delta=0.5" > optimal_design_model_${m}_alpha_${a}.sh
echo "qsub optimal_design_model_${m}_alpha_${a}.sh" >> run_initial_model_comparison.sh

done
done
