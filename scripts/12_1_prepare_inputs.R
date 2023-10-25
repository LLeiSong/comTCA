## --------------------------------------------
## Script name: allocate_cropland
## Purpose of script: allocate cropland for future expansion based on
## different criteria.
## Author: Lei Song
## Date Created: 2023-10-22
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------
## Criteria:
## start with lowest yield one
## start with random crop
## start with highest yield one

# Load libraries
library(terra)
library(sf)
library(dplyr)

# Set directories
agro_dir <- "data/agriculture"
bio_dir <- "data/biodiversity"
carbon_dir <- "data/carbon"
conn_dir <- "data/connectivity"

# Load data
# Prioritization and impact:
## - carbon cost
## - biodiversity
## - landscape connectivity
## - crop yield
## - cost distance
# Masks
## - slope-based farmable area
## - current cropland
## - waterbodies
## - forest/tree

# Slope-based farmable area
## Calculate slope
dem_path <- "data/elevation/elevation.tif"
output_dir <- "data/elevation/slope_4326.tif"
command <- sprintf(
    paste0("gdaldem slope %s %s -p -s 111120 -alg ZevenbergenThorne -compute_edges -co TILED=YES -co compress=lzw -co BIGTIFF=YES"),
    dem_path, output_dir)
system(command)

## Select farmable land
options <- '--NoDataValue=0 --co compress=lzw --type Byte'
eq_string <- sprintf("(A<20)")
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   output_dir, 'data/elevation/farmable.tif',
                   eq_string, options)
system(command)

## The above layer can be obtained from GEE to save resources
## Resample to fit with landcover layer
landcover_geo_fname <- "data/landcover/landcover_geo.tif"
landcover <- rast(landcover_geo_fname)
te <- ext(landcover); tr <- res(landcover)
options_warp <- paste0('-multi -wo NUM_THREADS=ALL_CPUS ',
                       '-co TILED=YES -co compress=lzw -co BIGTIFF=YES')
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s %s %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, 'data/elevation/farmable.tif', 
    'data/elevation/farmable_rsp.tif')
system(command)

## Remove current cropland, forest, waterbodies, and settlement
options_calc <- paste0('--co TILED=YES --co compress=lzw ',
                       '--co BIGTIFF=YES')
crop_msk <- "data/landcover/cropland_binary.tif"
forest_msk <- "data/landcover/treecover_binary.tif"
water_msk <- "data/landcover/water_binary.tif"
stm_msk <- "data/landcover/settlement_binary.tif"
eq_string <- "(A == 1) * (B == 0) * (C == 0) * (D == 0) * (E == 0)"
command <- sprintf('gdal_calc.py -A %s -B %s -C %s -D %s -E %s --outfile=%s --calc="%s" %s',
                   'data/elevation/farmable_rsp.tif', crop_msk, forest_msk, 
                   water_msk, stm_msk, "data/elevation/farmable_calibrated.tif", 
                   eq_string, options_calc)
system(command)

## summarize to 1km grids
template <- "data/landcover/cropland_ratio.tif"
farmable_ratio <- "data/elevation/farmable_ratio.tif"
te <- ext(rast(template)); tr <- res(rast(template))
command <- sprintf(
    paste0("gdalwarp -ot Float32 -te %s %s %s %s -tr %s %s %s %s %s -r average"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, "data/elevation/farmable_calibrated.tif", farmable_ratio)
system(command)

# Mask out not relevant areas
tradeoff_dir <- "data/tradeoff"
if (!dir.exists(tradeoff_dir)) dir.create(tradeoff_dir)
bry <- st_read("data/geoms/mainland_tanzania.geojson")

## Masks
farmable <- rast(farmable_ratio) %>% crop(bry) %>% mask(bry)
writeRaster(farmable, file.path(tradeoff_dir, "farmable_perc.tif"))
cropland <- rast("data/landcover/cropland_ratio.tif") %>% crop(bry) %>% mask(bry)
writeRaster(cropland, file.path(tradeoff_dir, "cropland_perc.tif"))

## Prioritization and impact


