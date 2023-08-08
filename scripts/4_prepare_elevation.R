# Title     : Prepare elevation dataset
# Objective : Prepare ALOS elevation map for Tanzania and global.
# Created by: Lei Song
# Created on: 07/26/23
# Note: first use GEE to query ALOS map of 100m using:
# var dataset = ee.ImageCollection("JAXA/ALOS/AW3D30/V3_2");
# var elevation = dataset.select('DSM').mosaic();
# 
# Export.image.toDrive({
#     image: elevation,
#     description: 'elevation',
#     folder: 'reconcile',
#     scale: 1000,
#     crs: 'EPSG:4326',
#     maxPixels: 1e13
# });
# Then go to the ALOS DSM official website, register and download the tiles

# Set directories
library(here)
library(sf)
library(stringr)
library(terra)
data_path <- here("data/elevation")
habitat_tz <- here("data/habitats/habitat_tz_refined_final.tif")

################### Get original DSM for Tanzania #####################
fnames <- list.files(file.path(data_path, "DSM"), recursive = TRUE, 
                     pattern = "DSM.tif", full.names = TRUE)
input_string <- paste(fnames, collapse = " ")

command <- sprintf("gdal_merge.py -o %s -co compress=lzw -co BIGTIFF=YES -co TILED=YES %s",
                   file.path(data_path, "elevation.tif"), input_string)
system(command)

# Resample to habitat map
template <- rast(habitat_tz)
te <- ext(template); tr <- res(template); rm(template)
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s -co compress=lzw -co BIGTIFF=YES -co TILED=YES %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    file.path(data_path, "elevation.tif"), 
    file.path(data_path, "elevation_resample.tif"))
system(command)

# This elevation map is ready to use for Tanzania.

############## Get DSM 100m for areas outside of Tanzania ###############
fnames <- list.files(file.path(data_path, "alos_dsm"), full.names = TRUE)
input_string <- paste(fnames, collapse = " ")

command <- sprintf("gdal_merge.py -o %s -co compress=lzw -co BIGTIFF=YES -co TILED=YES %s",
                   file.path(data_path, "elevation_outside.tif"), input_string)
system(command)

# Resample to habitat map
habitat_outside <- here("data/habitats/habitat_outside_lvl2.tif")
template <- rast(habitat_outside)
te <- ext(template); tr <- res(template); rm(template)
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s -co compress=lzw -co BIGTIFF=YES -co TILED=YES %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    file.path(data_path, "elevation_outside.tif"), 
    file.path(data_path, "elevation_outside_resample.tif"))
system(command)

# This elevation map is ready to use for global