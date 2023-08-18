library(terra)

po_path <- file.path("data", "ebird", "ebdata.csv")
df <- read.csv(po_path)

# https://www.usgs.gov/programs/gap-analysis-project/science/species-data-download

shp_file <- file.path("data", "us_regions", "regions", "GAP_Regions.shp")
p <- terra::vect(shp_file)



