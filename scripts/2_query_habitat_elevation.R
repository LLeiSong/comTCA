# Title     : Query suitable habitat and elevation of species from Red list.
# Objective : New script to query suitable habitat and elevation using 
#             RedList API (I finally get one!).
# Created by: Lei Song
# Created on: 07/26/23

# Load libraries
# install.packages("rredlist")
library(sf)
library(dplyr)
library(rredlist)
library(here)
library(dplyr)
library(terra)
library(stringr)
library(parallel)
select <- dplyr::select

# Set up the API key if need
# rl_use_iucn()

# Define paths
data_dir <- here("data/expert_range_maps")
range_map_dir <- file.path(data_dir, "clean_range_maps")
species_info_dir <- file.path(data_dir, "species_info")
if (!dir.exists(species_info_dir)) dir.create(species_info_dir)

# Read polygons
mammals <- read_sf(file.path(range_map_dir, "mammals_tz_relevant.geojson"))
amphibians <- read_sf(file.path(range_map_dir, "amphibians_tz_relevant.geojson"))
birds <- read_sf(file.path(range_map_dir, "birds_tz_relevant.geojson"))
reptiles <- read_sf(file.path(range_map_dir, "reptiles_tz_relevant.geojson"))

# Define the function to query
query_species_info <- function(sci_name){
    # Habitats
    hbts <- rl_habitats(sci_name)
    if (length(hbts$result) == 0){
        warning(sprintf("No habitat found for %s", sci_name))
        Sys.sleep(60)
        hbts <- rl_habitats(sci_name)}
    
    # Elevation
    sps <- rl_search(name = sci_name)
    elevation_upper <- sps$result$elevation_upper
    elevation_lower <- sps$result$elevation_lower
    
    # Return a list
    list(name = sci_name,
         habitats = hbts$result,
         elevation_upper = elevation_upper,
         elevation_lower = elevation_lower)}

# Make the query
mammals_info <- lapply(mammals$sci_name, query_species_info)
amphibians_info <- lapply(amphibians$sci_name, query_species_info)
birds_info <- lapply(birds$sci_name, query_species_info)
reptiles_info <- lapply(reptiles$sci_name, query_species_info)

# Save out
save(mammals_info, file = file.path(species_info_dir, "mammals_info.rda"))
save(amphibians_info, file = file.path(species_info_dir, "amphibians_info.rda"))
save(birds_info, file = file.path(species_info_dir, "birds_info.rda"))
save(reptiles_info, file = file.path(species_info_dir, "reptiles_info.rda"))
