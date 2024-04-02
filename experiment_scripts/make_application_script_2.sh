# No vary visits
exp_dir=application_2
rm -rf $exp_dir
mkdir $exp_dir
rm run_applied_model_comparison_2.sh
run_time=20
models="2"
sites="5 10"
#sites="5 10 15 20"
detection_p="0.2"
n_surveys="2 4"
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
cores=48
po_p=1
save_dir=data/ebird
for model in $models
do
for m in $sites
do
for n in $n_surveys
do
exp_name=oe_application_model-${model}_m-${m}_n-${n}_r-$(($cores*2))_p-${detection_p}
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N ${exp_dir}_${exp_name}
#PBS -l walltime=21:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR/0.5.2-foss-2022a-R-4.2.1
Rscript $WORKDIR/optimal_design_application.R --working_dir=$WORKDIR --save_dir=experimental_runs/$exp_dir --data_save_dir=$save_dir --exp_name=$exp_name --data_reps=$(($cores*2)) --m=$m --min_visits=1 --max_visits=$n --random_starts=10 --model_selection=$model --p=$detection_p --po_sample_prop=$po_p --exch_iter=40 --mcmc_iter=1000 --v_parallel --cores=$cores --run_time=$run_time" > $exp_dir/$exp_name.sh
echo "qsub $exp_dir/$exp_name.sh" >> run_applied_model_comparison_2.sh
done
done
done
