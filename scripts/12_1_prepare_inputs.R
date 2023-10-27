## --------------------------------------------
## Script name: prepare_inputs
## Purpose of script: prepare the inputs for the following analysis.
## Author: Lei Song
## Date Created: 2023-10-22
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------

# Load libraries
library(terra)
library(sf)
library(dplyr)

# Set directories
agro_dir <- "data/agriculture"
bio_dir <- "data/biodiversity"
carbon_dir <- "data/carbon"
conn_dir <- "data/connectivity"
# Mask out not relevant areas
tradeoff_dir <- "data/tradeoff"
if (!dir.exists(tradeoff_dir)) dir.create(tradeoff_dir)

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
    paste0("gdaldem slope %s %s -p -s 111120 -alg ZevenbergenThorne ",
           "-compute_edges -co TILED=YES -co compress=lzw -co BIGTIFF=YES"),
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

# Make the grid template for the record
bry <- st_read("data/geoms/mainland_tanzania.geojson")
template <- rast("data/landcover/cropland_ratio.tif")
template <- mask(template, bry)
template <- trim(template)
values(template)[!is.na(values(template))] <- 1:unlist((ncell(template) - global(template, "isNA")))
writeRaster(template, file.path(tradeoff_dir, "plannint_unit.tif"), datatype = "INT4U",
            gdal = c("COMPRESS=LZW"))

## summarize to 1km grids
template <- rast(file.path(tradeoff_dir, "plannint_unit.tif"))
farmable_ratio <- "data/elevation/farmable_ratio.tif"
te <- ext(template); tr <- res(template)
command <- sprintf(
    paste0("gdalwarp -ot Float32 -te %s %s %s %s -tr %s %s %s %s %s -r average"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, "data/elevation/farmable_calibrated.tif", farmable_ratio)
system(command)

## Areas
farmable <- rast(farmable_ratio) %>% crop(template)
writeRaster(farmable, file.path(tradeoff_dir, "farmable_perc.tif"))
cropland <- rast("data/landcover/cropland_ratio.tif") %>% crop(template)
writeRaster(cropland, file.path(tradeoff_dir, "cropland_perc.tif"))

# Protected areas
pas <- read_sf(
    file.path("data/protected_area", "WDPA_WDOECM_Jan2023_Public_TZA",
              "WDPA_WDOECM_Jan2023_Public_TZA.gdb"),
    layer = "WDPA_WDOECM_poly_Jan2023_TZA")
pas <- pas %>% 
    select(WDPAID, NAME, DESIG) %>% 
    rename(Geometry = SHAPE) %>% 
    filter(!DESIG %in% c("Marine Reserve", "Marine Park"))
pas <- rasterize(pas, template, values = 1, background = 0) %>% mask(template)
writeRaster(pas, file.path(tradeoff_dir, "protected_areas.tif"), 
            datatype = "INT1U", gdal = c("COMPRESS=LZW"))

## Prioritization and impact
# agriculture
current_yield <- file.path(agro_dir, "yield_calibrated_5crops_tz_1km.tif") %>% 
    rast() %>% extend(., ext(template)) %>% crop(template)
writeRaster(current_yield, file.path(tradeoff_dir, "agro_current_yield.tif"))
atn_yield <- file.path(agro_dir, "yield_attainable_5crops_tz_1km.tif") %>%
    rast() %>% extend(., ext(template)) %>% crop(template)
writeRaster(atn_yield, file.path(tradeoff_dir, "agro_attainable_yield.tif"))
travel_time <- rast(file.path(agro_dir, "travel_time.tif"))
travel_time[travel_time== -9999] <- NA
travel_time <- travel_time %>% resample(template) %>% 
    mask(template) %>% crop(template)
writeRaster(travel_time, file.path(tradeoff_dir, "agro_travel_time.tif"))

# biodiversity
bio_index <- rast(file.path(bio_dir, "BIp.tif")) %>% crop(template) %>% mask(template)
writeRaster(bio_index, file.path(tradeoff_dir, "bio_index.tif"))

# connectivity
conn_score <- rast(file.path(conn_dir, "circuit/mean_curmap_robust_norm.tif")) %>% 
    crop(template) %>% mask(template)
writeRaster(conn_score, file.path(tradeoff_dir, "conn_index.tif"))

# carbon
carbon_density <- rast(file.path(carbon_dir, "carbon_density.tif")) %>% resample(template) %>% 
    crop(template) %>% mask(template)
writeRaster(carbon_density, file.path(tradeoff_dir, "carbon_density.tif"))

# check if all layers align
fnames <- list.files(tradeoff_dir, full.names = TRUE)
inputs <- rast(fnames)
