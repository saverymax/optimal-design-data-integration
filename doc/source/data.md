# Data

This describes the data used in both the simulated and applied experiments. The data for the applied setting is described first. For all data downloaded in these instructions, it will be easiest to place them all in the same directory that you later provide to the script to run the design algorithm.

For the purpose of the blind review, we have provided all data necessary to run the application in the directory above the main source-code repository for this project. It is not possible to put the data in the main directory (ie how it would be structured in the github repository) as the datasets are too large. Once this code and data is made public, the data will be stored in a separate public-facing repository.

To run the code, you can use the data we provide. Alternatively, the instructions below describe how to download this data from scratch.

## Map data

It will be essential to have a file defining your region of interest before you can proceed with the steps below. We download a .gpkg file from <https://apps.nationalmap.gov/downloader/> for Tennesse in the work here, but you are welcome to define your own region. The code here uses read_sf() from the sf package in R, which is compatible with multiple map file types (.shp, .gpkg). Place this file in a base data directory that you will you use to store all the data for running the code the instructions here describe.

## Land cover data

Land cover data can be acquired from numerous sources. These include products from MODIS, Sentinel-II, Landsat, or Proba-V. In this documentation we
provide limited instructions for obtaining data the Proba-V product from Copernicus.

To download the PROBA-V satellite data, go to <https://land.copernicus.eu/global/products/lc> and use the world viewer (<https://lcviewer.vito.be/2019>). Select the appropriate section/s 
of the globe for your use-case. On this case, we simply need one section that appropriately covers Tennessee. Information regarding the product itself can be found in the documentation: <https://land.copernicus.eu/global/sites/cgls.vito.be/files/products/CGLOPS1_PUM_LC100m-V3_I3.4.pdf>. Important for this analysis are the classifications of each landcover type.


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

EVI is described here: <https://www.usgs.gov/landsat-missions/landsat-enhanced-vegetation-index#:~:text=EVI%20is%20similar%20to%20Normalized,in%20areas%20with%20dense%20vegetation>.

Another option is to use the Copernicus data viewer: See https://land.copernicus.eu/global/products/lc and use the NDVI product. Regardless, once you have downloaded the .tif file, place it within the data directory that you are storing your data related to this project.

## Elevation

In this work we use the AppEARS app to download elevation data for our region of interest, although many sources of elevation data are available at various resolutions. The interface to download the data is available at <https://appeears.earthdatacloud.nasa.gov/task/area>. Here you can provide a polygon, and download the ASTER GDEM elevation product. ASTER GDEM is described here: <https://lpdaac.usgs.gov/products/astgtmv003/>. For users interested in other sources of elevation data, see the FAQ here: <https://www.usgs.gov/faqs/where-can-i-get-global-elevation-data>.

## eBird 

The previous data processing steps should be completed before the eBird data is downloaded. In this work we use Brown-headed Nuthatch data within Tennessee.


The Brown-headed Nuthatch dataset is first downloaded from the eBird site <https://science.ebird.org/en/use-ebird-data/download-ebird-data-products>. You must initally ask for access to the eBird data. Once access is granted, when we download the data we filter first by species (brown-headed nuthatch) and year (2019) via the downloader user interface.

Once the eBird data is downloaded and the previously mentioned data sources are also downloaded, we can run the script that prepares the data for use within the exchange algorithm. Importantly, change the ```--base_data_dir``` to the location where you have saved all your data. To run this processing with the data we use in this project, leave the other CLI arguments as they are and place the data we provide in your base data directory.
```
Rscript data_utils/process_ebird.R --working_dir=. --base_data_dir=your/base/data/directory --ebird_data_dir="ebd_US_bnhnut_201901_201912_smp_relJul-2023" --map_file="us_states/GOVTUNIT_Tennessee_State_GPKG/GOVTUNIT_Tennessee_State_GPKG.gpkg" --landcover_file="copernicus_landcover/W100N40_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif" --modis_file="modis_landcover_dynamics/MCD12Q2.061_EVI_Area_0_doy2019001_aid0001.tif" --elevation_file="elevation_aster/ASTGTM_NC.003_ASTER_GDEM_DEM_doy2000061_aid0001.tif" --save_dir="data/ebird" --data_pack --pp_fit --pp_diagnostic
```
The process_ebird.R script performs a number of data processing steps using the R package auk. These steps include verifying the correct time period is filtered, selecting complete checklists (those with all species reported), and only those with stationary and travelling protocols (excluding historical and incidental). For pure PO data modelling we can use the incidental observations as well, but for this occupancy analysis we'll stick to the 2 protocols.

Once the auk steps are performed, there is also a sequence of spatial data processing steps as well, to create the discretized covariates for the site-occupancy model and the non-homogeneous poisson process. After these steps are performed, data_pack.RDS is saved, which contains all data necessary to run the applied experiments without relying on the original data files such as the elevation and EVI tifs.

TODO: OOPS, once I run process ebrid, I don't need to call anything more. Need to fix this is the scripts.

## Simulated data

Describe where the generated PO data is.
The simulated data that you need to run the code associated with the simulated experiments can be found in the directory within this repository ```data/sim_data```. We also provide instructions to generate this data yourself.

To run the simulated data experiments, we have to first run the script generate_po_data.R script. The bash script generate_po_datasets.sh contains the command to do so. The reason that we have to run this script is to create the PO datasets that will be used throughout the optimal design. But the same parameters that create the PO data must also be used to create the PA data. We pregenerate the PO data so that any experiment that particular combination of parameters can just load the correct dataset. To generate the PO dataset, we can run

```
Rscript generate_po_data.R --working_dir=. --exp_name=po_gen --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --intensity_func="donut" --bias_fun="exponential" --data_reps=96
```

The bash script in the experiment_scripts directory (experiment_scripts/generate_po_datasets.sh) will run the R script with all parameter combinations that we use in the experiments in the associated paper.

Then, we can run the optimal design using a particular PO simulated dataset. This is admittedly a bit difficult since we have to specify the same parameter sets that are used to generate the PO data in the design script, where these same parameters (alpha and beta) will be used to generate the PA data. For example given that we have ```--alpha=-2 --beta=0.5 --gamma=1 --delta=0.25``` above, in the design we need to run the script with
```
Rscript optimal-design-data-integration/optimal_design_site_occ.R --working_dir=. --exp_name=oe_model-4_m-5_n-10_r-96_intns-donut_a--2_b-0.5_g-1_d-0.25_p-0.2 --data_reps=96 --m=5 --min_visits=1 --max_visits=10 --vary_visits --model_selection=4 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --use_sim_po --po_data_file=po_gen_ints-donut_a=-2_b=0.5_g=1_d=0.25.Rds --intensity_func="donut" --bias_func="exponential" --v_parallel --cores=48 --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --p=0.2 --aux_cor=0.8
```
for example. We also need to make sure we use the --use_sim_po flag with the correct dataset for the specified parameters:
```
--po_data_file=po_gen_ints-donut_a=-2_b=0.5_g=1_d=0.25.Rds
```
