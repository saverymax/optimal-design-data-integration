#############################################################
# Module for processing PO Ebird data and covariates for 
# optimal design usage
#############################################################

library(terra)
library(tidyterra)
library(sf)
library(auk)

load_ebird <- function(base_data_dir, ebd_download_dir){
  # Load pre-processed PO data
  nuthatch <- read_csv(file.path(base_data_dir, ebd_download_dir, "nuthatch_filtered_for_occ.csv")) %>% 
    mutate(year = year(observation_date),
           # occupancy modeling requires an integer response
           species_observed = as.integer(species_observed))
  
}

load_landcover <- function(){
}


load_evi <- function(){
}

create_pp_grid <- function(){
  
  # Then plot over the state
  # Set up a grid first
  # Shoudl be 2500 x 2500 meters
  crs(state_bound)
  state_grid <- state_bound %>% st_make_grid(cellsize=c(2500,2500), what = "polygons", crs = "ESRI:102003")
  state_grid
  # get points in state
  # One way to do it
  # state_pp_within <- st_within(nuthatch_sf, state_bound, prepared = T, sparse=F)
  # state_pp <- nuthatch_sf[state_pp_within]
  state_pp <- st_intersection(nuthatch_sf, state_bound)
  class(state_pp)
  # Using [] to select the cells is the way to go
  subgrid <- state_grid[state_bound]
  plot(state_bound)
  plot(state_pp, pch = 19, cex = 0.5, col = alpha("orange", 0.5), add = TRUE)
  plot(subgrid, col = alpha("black", 0.0001), add=T)
  #plot(state_grid, col = alpha("black", 0.0001), add=T)
  
  # More advanced using tidyterra and ggplot
  p <- ggplot() + 
    #geom_spatvector(data=state_grid, fill = 'transparent', colour="lightblue") +
    geom_sf(data = state_bound, color=alpha("white",0.9)) + 
    geom_sf(data=state_pp, color=alpha("#FFC81C",0.2), size=0.5)+
    ggtitle("Brown-headed Nuthatch checklists in Tennessee") +
    theme_minimal() +
    theme(text=element_text(size=10)) +
    coord_sf()
  print(p)
  fig_name="data/ebird/ebd_checklist_nuthatch.png"
  # To save the degree symbol, we use cairo: 
  # https://www.andrewheiss.com/blog/2017/09/27/working-with-r-cairo-graphics-custom-fonts-and-ggplot/
  # This is also used for embedding custom fonts in pngs/svgs in general
  # To show it in R studio, change graphics backend to cairo
  ggsave(fig_name, plot=p, dpi=300, width=15, height=8, units="cm", bg="white", device="png", type="cairo")
  
  
  # Then create an intensity map of the points per cell
  point_counts <- st_intersects(subgrid, state_pp, sparse=F)
  dim(point_counts)
  # Count the number of points in each cell
  pp_counts <- apply(point_counts, MARGIN=1, FUN=sum)
  length(pp_counts)
  which(pp_counts>0)
}
