# Title     : Subset expert range maps
# Objective : To get the species relevant to Tanzania from the downloaded
#             expert range maps.
# Created by: Lei Song
# Created on: 07/26/23
# Note      : Expert range maps were queried from Red List as the data
#             last updated : 9th December 2022.

# Load libraries
library(sf)
library(dplyr)
library(here)
library(dplyr)
library(terra)
library(stringr)
library(purrr)
select <- dplyr::select

# Define paths
data_dir <- here("data/expert_range_maps")
range_map_dir <- file.path(data_dir, "clean_range_maps")
if (!dir.exists(range_map_dir)) dir.create(range_map_dir)

# Turn off s2
sf_use_s2(FALSE)

# First, clean the expert range map polygons
bry <- read_sf(here("data/geoms/tanzania_coarse_bry.geojson"))

# Define the function to avoid meaningless duplication
clean_species <- function(species, bry, dst_path){
    message('Clean up expert range maps.')
    ## Subset the interested ones
    ## presence: extant (1), Possibly extant (2 or 3)
    ##           others: Possibly Extinct (4), and Extinct (5)
    ## origin: native (1) and reintroduced (2), and Assisted Colonisation (6)
    ##         others: introduced (3), Vagrant (4), and Origin Uncertain (5)
    ## season: resident (1), breeding season (2) and non-breeding season (3)
    ##         others: Passage(4), and Seasonal Occurrence Uncertain (5)
    ## category: exclude EW and EX
    species <- species %>% 
        filter(marine == "false") %>% 
        filter(presence %in% 1:3) %>% 
        filter(origin %in% c(1:2, 6)) %>% 
        filter(seasonal %in% 1:3) %>% 
        filter(!category %in% c("EW", "EX"))
    
    # No issues, then subset useful information only
    species <- species %>% 
        select(sci_name, category)
    
    # Get species for tanzania
    message("Get species relevant to study area")
    species_tz <- species %>% 
        slice(unique(unlist(st_intersects(bry, .))))
    
    # Mosaic polygons for the same species
    species_tz <- do.call(
        rbind, lapply(unique(species_tz$sci_name), function(sname){
            this_species <- species_tz %>% filter(sci_name == sname)
            if (nrow(this_species) > 1){
                st_union(this_species) %>% st_sf() %>% 
                    mutate(sci_name = this_species$sci_name[1],
                           category = this_species$category[1]) %>% 
                    select(sci_name, category, geometry)
            } else this_species}))
    message(sprintf("No.%s of species left.", nrow(species_tz)))
    
    # Save out
    message(sprintf("Save out to %s", dst_path))
    write_sf(species_tz, dst_path)
}

# Get species exist in Tanzania
## -- Terrestrial mammals
mammals <- read_sf(
    file.path(data_dir, "MAMMALS_TERRESTRIAL_ONLY", 
              "MAMMALS_TERRESTRIAL_ONLY.shp"))

## -- amphibians
amphibians <- read_sf(file.path(data_dir, "AMPHIBIANS", "AMPHIBIANS.shp"))

## -- Birds
birds <- read_sf(file.path(data_dir, "BOTW.gdb"), 
                 layer = "All_Species") %>% 
    left_join(read_sf(
        file.path(data_dir, "BOTW.gdb"), 
        layer = "Taxonomic_checklist"),
        by = c("sisid" = "SISID")) %>% 
    rename(category = RL_Category) %>% 
    mutate(marine = "false") # format for the function

# Fix some MULTISURFACE polygons
birds_good <- birds %>% 
    filter(st_is(. , "MULTIPOLYGON") | st_is(., "POLYGON"))
birds_bad <- birds %>% 
    filter(!(st_is(. , "MULTIPOLYGON") | st_is(., "POLYGON")))
st_write(birds_bad, file.path(data_dir, "birds_bad.geojson"), delete_dsn = TRUE)
birds_bad <- st_read(file.path(data_dir, "birds_bad.geojson"))
birds <- rbind(birds_good %>% rename(), birds_bad)
rm(birds_good, birds_bad)
file.remove(file.path(range_map_dir, "birds_bad.geojson"))

# Subset mammals, amphibians, and birds
clean_species(mammals, bry, here(range_map_dir, "mammals_tz_relevant.geojson"))
rm(mammals)
clean_species(amphibians, bry, here(range_map_dir, "amphibians_tz_relevant.geojson"))
rm(amphibians)
clean_species(birds, bry, here(range_map_dir, "birds_tz_relevant.geojson"))
rm(birds)

## -- Reptiles, do separately
category <- read.csv(
    file.path(data_dir, "GARD", "category.csv"),
    stringsAsFactors = FALSE) %>% 
    select(binomial, predicted_cat)
reptiles <- read_sf(
    file.path(data_dir, "GARD", "Gard_1_7_ranges.shp")) %>% 
    left_join(category, by = "binomial") %>% 
    rename(sci_name = binomial, category = predicted_cat) %>% 
    filter(!category %in% c("EW", "EX"))
reptiles_tz <- reptiles %>% st_make_valid() %>% 
    slice(unique(unlist(st_intersects(bry, .))))

# Check the filters based on Red list polygons
reptiles_iucn <- read_sf(
    file.path(data_dir, "REPTILES", "REPTILES_PART1.shp")) %>% 
    select(-OBJECTID) %>% 
    rbind(read_sf(
        file.path(data_dir, "REPTILES", "REPTILES_PART2.shp"))) %>% 
    select(sci_name, marine, presence, origin, seasonal, category) %>% 
    st_drop_geometry()

reptiles_selected <- reptiles_iucn %>% 
    filter(marine == "false") %>% 
    filter(presence %in% 1:3) %>% 
    filter(origin %in% c(1:2, 6)) %>% 
    filter(seasonal %in% 1:3) %>% 
    select(sci_name) %>% unique()

reptiles_tz <- reptiles_tz %>% left_join(reptiles_selected, "sci_name") %>% 
    select(sci_name, category)

write_sf(reptiles_tz, file.path(range_map_dir, "reptiles_tz_relevant.geojson"))
