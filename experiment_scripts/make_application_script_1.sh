# No vary visits
exp_dir=application_1
rm -rf $exp_dir
mkdir $exp_dir
rm run_applied_model_comparison_1.sh
run_time=20
models="1 2"
m=10
detection_p="0.2"
po_sample_prop="0.01 0.05 0.1 0.2 0.3 0.5 0.7 1"
n_surveys="3"
WORKDIR=$VSC_DATA/projects/optimal_design_presence_only/optimal-design-data-integration
cores=48
BASE_DATA_DIR=$VSC_DATA/projects/data
map_file=us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg
landcover_file=copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif
modis_file=modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif
elev_file=elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif
save_dir=data/ebird
ebird_dir=ebd_US_bnhnut_201901_201912_smp_relJul-2023
for model in $models
do
for po_p in $po_sample_prop
do
exp_name=oe_application_model-${model}_m-${m}_n-${n_surveys}_r-$(($cores*2))_p-${detection_p}_po_prop-${po_p}
echo "#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N ${exp_dir}_${exp_name}
#PBS -l walltime=21:00:00
#PBS -l nodes=1:ppn=$cores
#PBS -l mem=100gb

module load CmdStanR/0.5.2-foss-2022a-R-4.2.1
Rscript $WORKDIR/optimal_design_application.R --working_dir=$WORKDIR --save_dir=experimental_runs/$exp_dir --base_data_dir=$BASE_DATA_DIR --data_save_dir=$save_dir --ebird_data_dir=$ebird_dir --map_file=$map_file --landcover_file=$landcover_file --modis_file=$modis_file --elevation_file=$elev_file --exp_name=$exp_name --data_reps=$(($cores*2)) --m=$m --min_visits=1 --max_visits=$n_surveys --random_starts=10 --model_selection=$model --p=$detection_p --po_sample_prop=$po_p --exch_iter=40 --mcmc_iter=1000 --v_parallel --cores=$cores --run_time=$run_time" > $exp_dir/$exp_name.sh
echo "qsub $exp_dir/$exp_name.sh" >> run_applied_model_comparison_1.sh
done
done
