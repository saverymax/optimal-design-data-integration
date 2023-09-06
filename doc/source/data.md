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