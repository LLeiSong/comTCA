# Title     : Script to do experiments for landscape connectivity evaluation
# Objective : Use virtual experiments to test and compare different 
#             methods/models for connectivity conservation planning.
# Created by: Lei Song
# Created on: 06/14/23

# Load packages and set directories
library(here)
library(sf)
library(terra)
library(dplyr)
test_dir <- here("data/test")
ele_dir <- here('data/elephant')

############# Make virtual experiments ##################
## Virtual landscape condition: 
## habitat patches, movement resistance, and disturbances
#########################################################

# Make boundaries
ext <- read_sf(file.path(ele_dir, "protected_area.geojson"))
write_sf(ext, file.path(test_dir, "pas.geojson"))
ext <- ext %>% slice(1) %>% select(name)
write_sf(ext, file.path(test_dir, "aoi.geojson"))

ext <- read_sf(file.path(test_dir, "aoi.geojson"))
pas <- read_sf(file.path(test_dir, "pas.geojson"))
pas <- st_crop(pas, ext)

# Make movement resistance
suit <- rast(file.path(ele_dir, "landscape_utility_1km_integrated.tif"))
suit <- crop(suit["prediction"], ext)
writeRaster(suit, file.path(test_dir, "suitability.tif"))


