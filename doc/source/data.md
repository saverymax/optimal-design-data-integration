# Data

## eBird 

The optimal design procedure available in this project can be applied to real or simulated data. In the associated paper with 
this documentation, we apply the algorithm to eBird data in Tennessee. 

To download and process the data for the experiments or your own usage, you can follow the following steps. 

	1. 
	2. 
	3. When the data is downloaded and placed in the appropriate directory, run the processing script process_ebird.R
		```Rscript process_ebird.R```

## Land cover data

Land cover data can be acquired from numerous sources. These include products from MODIS, Sentinel-II, Landsat, or Proba-V. In this documentation we
provide limited instructions for obtaining data from MODIS or Proba-V.

There are numerous ways to download MODIS data. 
	
	1. First save a .shp file of the area you want to do this in R, we can do it in such a way:
		```
		library(terra)
		library(rnaturalearth)

		us_map <- ne_countries(country = "united states of america", returnclass = "sf") %>% st_transform(crs=map_proj)
		us_vect <- vect(us_map)
		writeVector(us_vect, "us.shp")
		```

	2. Then, navigate to https://lpdaac.usgs.gov/products/mcd12q2v061/ and select ACCESS DATA.

	3. Select the AppEEARS tool download data button. 

	4. At the top of the menu select Extract. This will take you to a sort of form. The key thing to do here is to drag your .shp (folder)
		to select the area of interest. Next, select your data and product. For example, in this work we use "Combined MODIS land cover dynamics"
		and "Combined MODIS land cover type"
		and selected all layers.

	5. Submit your request and wait for the order to complete.


Another option is to use the Copernicus data viewer: See https://land.copernicus.eu/global/products/lc and https://lcviewer.vito.be/2019

We can also access 
## Simulated data

To run the simulated data experiments, we have to first run the script generate_po_data.R script. The bash script generate_po_datasets.sh contains the command to do so. The reason that we have to run this script is to create the PO datasets that will be used throughout the optimal design. But the same parameters that create the PO data must also be used to create the PA data. We pregenerate the PO data so that any experiment that particular combination of parameters can just load the correct dataset. To generate the PO dataset, we can run

```
Rscript generate_po_data.R --working_dir=. --exp_name=po_gen --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --intensity_func="donut" --bias_fun="exponential" --data_reps=96
```

The bash script will run the R script with all parameter combinations that we use in the experiments.

Then, we can run the optimal design using a particular PO simulated dataset. This is admittedly a bit difficult since we have to specify the same parameter sets that are used to generate the PO data in the design script, where these same parameters (alpha and beta) will be used to generate the PA data. For example given that we have ```--alpha=-2 --beta=0.5 --gamma=1 --delta=0.25``` above, in the design we need to run the script with
```
Rscript optimal-design-data-integration/optimal_design_site_occ.R --working_dir=. --exp_name=oe_model-4_m-5_n-10_r-96_intns-donut_a--2_b-0.5_g-1_d-0.25_p-0.2 --data_reps=96 --m=5 --min_visits=1 --max_visits=10 --vary_visits --model_selection=4 --random_starts=3 --exch_iter=20 --mcmc_iter=1000 --use_sim_po --po_data_file=po_gen_ints-donut_a=-2_b=0.5_g=1_d=0.25.Rds --intensity_func="donut" --bias_func="exponential" --v_parallel --cores=48 --alpha=-2 --beta=0.5 --gamma=1 --delta=0.25 --p=0.2 --aux_cor=0.8
```
for example. We also need to make sure we use the --use_sim_po flag with the correct dataset for the specified parameters:
```
--po_data_file=po_gen_ints-donut_a=-2_b=0.5_g=1_d=0.25.Rds
```
