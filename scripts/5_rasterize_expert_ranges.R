# Title     : The script to rasterize expert range maps
# Objective : The main function to rasterize the expert range maps. This is the
#           : downstream script of rasterize_expert_ranges_scheduler.R
# Created by: Lei Song
# Created on: 07/26/23

# Load libraries
library(sf)
library(terra)
library(parallel)
library(optparse)
library(dplyr)

# Set directories
data_dir <- "/scratch/lsong36/comTCA/data"
dst_dir <- file.path(data_dir, "expert_range_maps/clean_range_rasters")
if (!dir.exists(dst_dir)) dir.create(dst_dir)

# Parse inline arguments
option_list <- list(
    make_option(c("-s", "--species"), 
                action = "store", type = 'character',
                help = "The species to process."),
    make_option(c("-t", "--taxon"), 
                action = "store", type = 'character',
                help = "The taxon of the species."))
opt <- parse_args(OptionParser(option_list = option_list))
species <- opt$species
taxon <- opt$taxon

message(sprintf("Species %s in taxon %s.", species, taxon))

# Sub-directories
tz_dir <- file.path(dst_dir, "tanzania", taxon)
outside_dir <- file.path(dst_dir, "global", taxon)
if (!dir.exists(tz_dir)) dir.create(tz_dir, recursive = TRUE)
if (!dir.exists(outside_dir)) dir.create(outside_dir, recursive = TRUE)

# Rasterize the range polygons
## Grid template
habitat_tz <- rast(file.path(data_dir, "habitats/habitat_tz_refined_final.tif"))
habitat_outside <- rast(file.path(data_dir, "habitats/habitat_outside_lvl2.tif"))

# Read polygons
plys <- st_read(file.path(data_dir, "expert_range_maps/clean_range_maps",
                          sprintf("%s_tz_relevant.geojson", taxon)))

if (taxon != "birds"){
    # Write out temporary file
    ply <- plys %>% filter(sci_name == species)
    species_nm <- gsub(" ", "_", species)
    fname <- file.path(tempdir(), sprintf("%s.geojson", species_nm))
    st_write(ply, fname)
    
    # For Tanzania
    te <- ext(habitat_tz); tr <- res(habitat_tz)
    command <- sprintf(
        paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
               "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
        te[1], te[3], te[2], te[4], tr[1], tr[2], 
        fname, file.path(tz_dir, sprintf("%s.tif", species_nm)))
    system(command)
    
    Sys.sleep(30) # Wait for writing
    
    # For global
    te <- ext(habitat_outside); tr <- res(habitat_outside)
    command <- sprintf(
        paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
               "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
        te[1], te[3], te[2], te[4], tr[1], tr[2], 
        fname, file.path(outside_dir, sprintf("%s.tif", species_nm)))
    system(command)
    
    # Clean up
    file.remove(fname)
} else {
    # Split season for birds
    ## resident (1), breeding season (2) and non-breeding season (3)
    plys <- plys %>% filter(sci_name == species)
    species_nm <- gsub(" ", "_", species)
    
    for (season in unique(plys$seasonal)){
        message(sprintf("Season %s for this species.", season))
        
        season_code <- ifelse(season == 1, "R", ifelse(season == 2, "B", "N"))
        ply <- plys %>% filter(seasonal == season)
        fname <- file.path(
            tempdir(), sprintf("%s_%s.geojson", species_nm, season_code))
        st_write(ply, fname)
        
        # For Tanzania
        te <- ext(habitat_tz); tr <- res(habitat_tz)
        command <- sprintf(
            paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
                   "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
            te[1], te[3], te[2], te[4], tr[1], tr[2], fname, 
            file.path(tz_dir, sprintf("%s_%s.tif", species_nm, season_code)))
        system(command)
        
        Sys.sleep(30) # Wait for writing
        
        # For global
        te <- ext(habitat_outside); tr <- res(habitat_outside)
        command <- sprintf(
            paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
                   "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
            te[1], te[3], te[2], te[4], tr[1], tr[2], fname, 
            file.path(outside_dir, sprintf("%s_%s.tif", species_nm, season_code)))
        system(command)
        
        # Clean up
        file.remove(fname)
    }
}
