#!/bin/bash
#PBS -o /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -e /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/job_output/
#PBS -N oe_application_model-1_m-5_n-5_r-96_p-0.2
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=150gb

module load CmdStanR
Rscript /data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration/optimal_design_application.R --working_dir=/data/gent/459/vsc45956/projects/optimal_design_presence_only/optimal-design-data-integration --base_data_dir=/data/gent/459/vsc45956/projects/data --data_save_dir=data/ebird --ebird_data_dir=ebd_US_bnhnut_201901_201912_smp_relJul-2023 --map_file=us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg --landcover_file=copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif --modis_file=modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif --elevation_file=elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif --exp_name=oe_application_model-1_m-5_n-5_r-96_p-0.2 --data_reps=96 --m=5 --min_visits=1 --max_visits=5 --random_starts=3 --model_selection=1 --p=0.2 --exch_iter=20 --mcmc_iter=1000 --v_parallel --cores=48
