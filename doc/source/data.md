# Data

This describes the data used in both the simulated and applied experiments. The data for the applied setting is described first. For all data downloaded in these instructions, it will be easiest to place them all in the same directory that you later provide to the script to pre-process the data. 

HOWEVER, it is not necessary to follow these steps if you would 
like to use our pre-processed data set. Within this repository, the data can be found in the ```optimal-design-data-integration/data/ebird/data_pack.RDS``` file. This should be provided to the ```optimal_design_application.R``` script as described in the usage section of this documentation. You can skip to the [usage section](usage) in this case.

We have also provided all data necessary to reconstruct the ```data_pack.RDS``` in the zipped ```raw_data``` folder in the figshare repository <https://doi.org/10.6084/m9.figshare.25450108>. The data processing proceeds as: download -> pre-process -> run experiments. We have provided instructions and preliminary files for each step of the process, so the user can choose where to start to replicate the work. To run the experiment code, you can use the data we provide via the file ```optimal-design-data-integation/data/ebird/data_pack.RDS```. Alternatively, the instructions below describe how to download, process, and create this file and data from scratch.  

## Additional packages

If you would like to recreate the data used in this work, a few extra packages are required:
- auk
- exactextracr
- lubridate


## Map data

It will be essential to have a file defining your region of interest before you can proceed with the steps below. We download a .gpkg file from <https://apps.nationalmap.gov/downloader/> for Tennesse in the work here, but you are welcome to define your own region. The code here uses read_sf() from the sf package in R, which is compatible with multiple map file types (.shp, .gpkg). Place this file in a base data directory that you will you use to store all the data for running the design code.

## Land cover data

Land cover data can be acquired from numerous sources. These include products from MODIS, Sentinel-II, Landsat, or Proba-V. In this documentation we
provide limited instructions for obtaining the Proba-V product from Copernicus.

To download the PROBA-V satellite data, go to <https://land.copernicus.eu/global/products/lc> and use the world viewer (<https://lcviewer.vito.be/2019>). Select the appropriate section/s 
of the globe for your use-case. In our case, we simply need one section that appropriately covers Tennessee. Information regarding the product itself can be found in the documentation: <https://land.copernicus.eu/global/sites/cgls.vito.be/files/products/CGLOPS1_PUM_LC100m-V3_I3.4.pdf>. Important for this analysis are the classifications of each landcover type.


## Enhanced Vegetation Index (EVI)

There are various vegetation products available online. We use the MODIS product from NASA. There are numerous ways to download MODIS data. We follow the steps below:
	
	1. First save a .shp file of the area you want to do this in R, we can do it in such a way:
		```
		library(terra)
		library(rnaturalearth)

		us_map <- ne_countries(country = "united states of america", returnclass = "sf") %>% st_transform(crs=map_proj)
		us_vect <- vect(us_map)
		writeVector(us_vect, "us.shp")
		```

	2. Then, navigate to <https://lpdaac.usgs.gov/products/mcd12q2v061/> and select ACCESS DATA.

	3. Select the AppEEARS tool download data button. 

	4. At the top of the menu select Extract. This will take you to a form. Drag your .shp (folder)
		to select the area of interest. Next, select your data and product. For example, in this work we use "Combined MODIS land cover dynamics". Select all layers.

	5. Submit your request and wait for the order to complete.

EVI is described here: <https://www.usgs.gov/landsat-missions/landsat-enhanced-vegetation-index>.

Another option is to use the Copernicus data viewer: See <https://land.copernicus.eu/en/products/vegetation> and use the NDVI product. Regardless, once you have downloaded the .tif file, place it within the data directory that you are storing your data related to this project.

## Elevation

In this work we use the AppEARS app to download elevation data for our region of interest, although many sources of elevation data are available at various resolutions. The interface to download the data is available at <https://appeears.earthdatacloud.nasa.gov/task/area>. Here you can provide a polygon, and download the ASTER GDEM elevation product. ASTER GDEM is described here: <https://lpdaac.usgs.gov/products/astgtmv003/>. For users interested in other sources of elevation data, see the FAQ here: <https://www.usgs.gov/faqs/where-can-i-get-global-elevation-data>.

## eBird 

The previous data processing steps should be completed before the eBird data is downloaded. In this work we use Brown-headed Nuthatch data within Tennessee.

The Brown-headed Nuthatch dataset is first downloaded from the eBird site <https://science.ebird.org/en/use-ebird-data/download-ebird-data-products>. You must initally ask for access to the eBird data. Once access is granted, when we download the data we filter first by species (brown-headed nuthatch) and year (2019) via the downloader user interface.

## Processing data

Once the eBird data is downloaded and the previously mentioned data sources are also downloaded, we can run the script that prepares the data for use within the exchange algorithm. Importantly, in the command shown below, change the ```--base_data_dir``` to the location where you have saved all your data. **For example**, once you have unzipped the ```raw_data``` folder in the figshare repository, you will have a directory ```data``` with a subdirectory ```data/copernicus_landcover```, for example. As ```base_data_dir```, provide the path to (including) ```data``` that contains the rest of the subdirectories. 

To run this processing with the data we use in this project, leave the other CLI arguments as they are and give the path to ```data``` from the unizpped ```raw_data``` file as the base data directory (for arguments ```--ebird_data_dir, --map_file, --landcover_file, --modis_file, --elevation_file```). Navigate to the ```optimal-design-data-integration``` source code directory if you have not done so yet. Then run
```
Rscript data_utils/process_ebird.R --working_dir=. --base_data_dir=your/base/data/directory --ebird_data_dir="ebd_US_bnhnut_201901_201912_smp_relJul-2023" --map_file="us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg" --landcover_file="copernicus_landcover/Discrete-Classification-map_EPSG-4326.tif" --modis_file="modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif" --elevation_file="elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif" --save_dir="data/ebird" --data_pack --pp_fit --pp_diagnostic --auk_process --gamma_integration
```
The process_ebird.R script performs a number of data processing steps using the R package auk. These steps include verifying the correct time period is filtered, selecting complete checklists (those with all species reported), and only those with stationary and travelling protocols (excluding historical and incidental). For pure PO data modelling we can use the incidental observations as well, but for this occupancy analysis we'll stick to the 2 protocols.

Once the auk steps are performed, preliminary files are saved (```"nuthatch_filtered_for_occ.csv"``` being the most important one), and it is not necessary to use the ```auk_process``` option if you run the script again. It does take some time to run this first step though. Once these files are saved, there is a sequence of spatial data processing steps as well, to create the discretized covariates for the site-occupancy model and the non-homogeneous poisson process. After these steps are performed, ```data_pack.RDS1``` is saved, which contains all data necessary to run the applied experiments without relying on the original data files such as the elevation and EVI tifs. If ```--pp_fit``` is included, the posterior from the non-homogeneous poisson process (NHPP) is fit and if ```--gamma_integration``` is included, the posterior of alpha is averaged over gamma to remove the correlation between the two parameters. The NHPP posterior is required in order to have an estimate of the parameters to generate the PA data for the integration for approximation of the expected utility function, so you must include this option the first time you run the script. Gamma integration is optional but recommended as without it the alpha estimate is unreliable and may result in additional noise in the generated PA datasets used in the exchange algorithm.

Once the ```process_ebird.R``` script is sucessefully run and the ```data_pack.RDS``` file generated, you are ready to run the optimal design algorithm applied to the Brown-headed nuthatch example or your own use-case.

For further reference, the CLI options for the ```process_ebird.R``` script are listed below:
```
Options:
        -h, --help
                Show this help message and exit

        --working_dir=WORKING_DIR
                Path to the directory containing code to source for the main script

        --base_data_dir=BASE_DATA_DIR
                Path to the directory containing covariate data

        --ebird_data_dir=EBIRD_DATA_DIR
                Name of directory containing processed ebird data csv's, within basedir

        --save_dir=SAVE_DIR
                Directory to save data within basedir

        --map_file=MAP_FILE
                Nmae of file containing processed US geopackage fil

        --landcover_file=LANDCOVER_FILE
                Name of landcover tif

        --modis_file=MODIS_FILE
                Name of modis EVI tif

        --elevation_file=ELEVATION_FILE
                Name of elevation tif

        --data_reps=DATA_REPS
                Number of dataset reps for criterion estimation

        --cell_size=CELL_SIZE
                Size of one side of cell in point process grid

        --auk_process
                Boolean for running the initial auk filtering steps to generate smallers csv's

        --data_pack
                Boolean for generating the data pack if not already generated

        --pp_fit
                Boolean to fit the Point Process posterior after data saving steps

        --pp_diagnostic
                Boolean for printing diagnostics for point process model

        --gamma_integration
                Boolean for integrating over gamma in the PP
```


## Simulated data

The above instructions describe how to pre-process the data for the applied case-study. To run the simulated experiments, the data processing is simpler but there are still a few pre-processing steps we need to take. If you would like to skip these steps, we have provided the data you need to get started with the optimal design in ```data/sim_data```. See the instructions in the usage section in that case.  

To generate the data for the simulated experiments yourself, you have to first run the script ```generate_po_data.R```. The bash script ```experiment_scripts/generate_po_datasets.sh``` contains the command that generates all PO datasets used in this work. However, the same parameter settings that create the PO data must also be used to create the PA data, so that the two processes match. We pre-generate all PO data for a number of parameter combinations so that any experiment for a particular combination of parameters can just load the correct PO dataset, instead of having to generate a new dataset everytime we call the script. Therefore, to generate one particular PO dataset, we can run
```
Rscript generate_po_data.R --working_dir=. --exp_name=po_gen --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --intensity_func="donut" --bias_fun="exponential" --data_reps=96
```
This will simulate PO data for the particular values of the parameters as seen in the command. The bash script ```experiment_scripts/generate_po_datasets.sh``` will run the R script with all parameter combinations that we use in the experiments in the associated paper.

Then, once the dataset/s are generated, we can run the optimal design using a particular PO simulated dataset. This is admittedly a bit difficult since we have to specify the same parameter sets that are used to generate the PO data in the design script, where these same parameters (alpha and beta) will be used to generate the PA data. For example given that we have ```--alpha=-2 --beta=0.5 --gamma=1 --delta=0.25``` above, in the design we need to run the optimal design script with
```
Rscript optimal_design_site_occ.R --working_dir=. --exp_name=params_a--2_b-0.5_g-1_d-0.25_p-0.2 --data_reps=96 --m=5 --min_visits=1 --max_visits=10 --vary_visits --model_selection=2 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --use_sim_po --po_data_file=po_gen_ints-donut_a=-2_b=0.5_g=1_d=0.25.Rds --intensity_func="donut" --bias_func="exponential" --v_parallel --cores=48 --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --p=0.2 --aux_cor=0.8
```
for example. This is explained more in the usage section. We need to make sure we use the --use_sim_po flag and the correct dataset for the specified parameters given to ```--po_data_file```:
```
--po_data_file=po_gen_ints-donut_a=-2_b=0.5_g=1_d=0.25.Rds
```
This is the dataset generated by the call of ```generate_po_data.R```. Once you have done this for the parameter values you are interested in, you are ready to run the simulated design experiments (or you can just use the pre-processed data already provided as there are many dataset combinations already there).

For further reference, the CLI options for the ```generate_po_data.R``` script are listed below:
```
Options:
        -h, --help
                Show this help message and exit

        --working_dir=WORKING_DIR
                Path to the directory containing script

        --exp_name=EXP_NAME
                Base name to save data

        --alpha=ALPHA
                Intercept for intensity

        --beta=BETA
                Slope for intensity

        --gamma=GAMMA
                Intercept for bias

        --delta=DELTA
                Slope for bias

        --intensity_func=INTENSITY_FUNC
                Intensity function for sampling surface

        --sd=SD
                Standard deviation for donut intensity surface

        --bias_func=BIAS_FUNC
                Bias function for sampling surface

        --area=AREA
                Area of region D

        --k=K
                Number of sites along one side of grid

        --data_reps=DATA_REPS
                Number of dataset reps for criterion estimation
```
### Simulate data for misspecification

The PO data for the misspecification experiments is generated separately from the previous section. To generate this data, you must run ```generate_po_data_misspec.R```, for example 
```
Rscript generate_po_data_misspec.R --working_dir=. --save_dir="misspec" --exp_name=po_gen_misspec --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --epsilon=1.5 --intensity_func="donut" --sd=5 --bias_fun="exponential" --data_reps=96 --gamma_reps=200
```
It is the same as above, only now we provide an ```--epsilon``` argument as the known parameter value for the additional covariate in the data generating model of the bias. The number of iterations for the Monte Carlo integration over gamma is also set by ```--gamma_reps```. These additional options are listed below, but are otherwise the same as in the previous section.
```
Options:
        --epsilon=EPSILON 
                Slope for covariate inducing misspec

        --gamma_reps=GAMMA_REPS
                Number of MC iterations for integration over gamma

        --save_dir=SAVE_DIR
                Name of folder to save data within the data/sim_data directory
```
As before, the bash script ```experiment_scripts/generate_po_datasets_misspec.sh``` will generate data and NHPP posteriors for the parameter combinations that we use in the experiments in the paper.

What ```generate_po_data_misspec.R``` does is a little different than in the previous section. Firstly, it generates and saves PO dataset using an additional covariate that will be not be included in the models in the design script. Multiple PO datasets are generated for one fixed set of parameters, each for a different level of exponential decay (rho in the paper, with ```--epsilon=1.5```) hardcoded in the script. This is to test the sensitivity of the design algorithm to misspecification in the working model. The second thing this script does is estimate alpha from the PO data, to be used for PA dataset generation in the design algorithm. This is done via Monte Carlo integration over gamma, by iteratively drawing from a prior placed on gamma, holding gamma to this constant value in the NHPP, and averaging multiple NHPP estimates for alpha when iteratively fixing gamma. This posterior average estimate is saved to a ```.Rds``` file, such as ```pp_posterior_po_gen_ints-donut_peak=5_a=-2_b=0.5_g=1_d=0.25_e=0.1_int``` in the example used above. ```e=0.1``` refers to the amount of exponential decay. ```int``` refers to the use of MC integration. Both ```int``` and ```no_int``` versions will be saved, where ```no_int``` just takes the estimate of alpha directly from the NHPP even though this estimate is highly correlated with gamma. This posterior file, as well as the PO dataset, will need to be provided to the design script for PA generation. The datasets created and used in the paper are already available in ```data/sim_data/misspec```.