# Title     : Make the image cube of environmental variables
# Objective : Prepare environmental variables for SDM. 
# Created by: Lei Song
# Created on: 09/09/23
# Note.     : Need to prepare the vectors of buildings, big roads, waterbodies
#             rivers and streams from OpenStreetMap and Open Building.

# Load libraries
library(sf)
library(here)
library(terra)
library(dplyr)
library(rgbif)
library(itsdm)
library(corrplot)

sf_use_s2(FALSE)

# Set directories
data_dir <- here("data")
conn_data_dir <- here("data/connectivity")
dst_dir <- file.path(conn_data_dir, 'variables')
if (!dir.exists(dst_dir)) dir.create(dst_dir)

# Get the boundary
bry <- read_sf(here("data/geoms/tanzania_coarse_bry.geojson"))

# Climatic variables
scale <- 0.5
bios <- worldclim2(var = "bio", res = scale,
                   bry = bry, path = tempdir())
write_stars(
    bios, file.path(dst_dir, sprintf("bios_%s.tif", scale)))

vals <- values(bios)
vals <- data.frame(vals) %>% na.omit()
corrplot(corr = cor(vals, method = "spearman"), 
         method = "square", type = "lower",
         diag = FALSE, addCoef.col = "black", 
         number.cex = 0.7, tl.cex = 0.7)

# Bio 1 = Annual Mean temperature
# Bio 4 = Temperature seasonality
# Bio 7 = Temperature Annual Range (BIO5-BIO6)
# Bio 12 = Annual Precipitation
# Bio 14 = Precipitation of Driest Month
# Bio 15 = Precipitation seasonality (coefficient of variation)
# Bio 18 = Precipitation of Warmest Quarter
vals <- values(bios)
vals <- data.frame(vals) %>% na.omit()
vals <- vals %>% select(paste0("bio", c(1, 4, 7, 12, 14, 15, 18)))
corrplot(corr = cor(vals, method = "spearman"), 
         method = "square", type = "lower",
         diag = FALSE, addCoef.col = "black", 
         number.cex = 0.7, tl.cex = 0.7)

bios <- subset(bios, c(1, 4, 7, 12, 14, 15, 18))
writeRaster(bios, file.path(dst_dir, sprintf("bios_%s_selected.tif", scale)))

# NDVI, 13 - 22 average
# Processed in GEE and download
fnames <- list.files(dst_dir, pattern = "lansat8", full.names = TRUE)
fnames <- fnames[c(1, 3, 2)]
ndvis <- rast(fnames)
names(ndvis) <- c("ndvi_dry_season", "ndvi_wet_season", "ndvi_seasonality")
writeRaster(ndvis, file.path(dst_dir, "ndvi_dry_wet_sd.tif"))

# Natural & non-natural land coverage
# - Cropland (one of anthropological impacts)
options_calc <- paste0('--hideNoData ',
                       '--co TILED=YES --co compress=lzw ',
                       '--co BIGTIFF=IF_NEEDED --co NUM_THREADS=ALL_CPUS')
options_warp <- paste0('-multi -wo NUM_THREADS=ALL_CPUS ',
                       '-co TILED=YES -co compress=lzw -co BIGTIFF=IF_NEEDED')
landcover <- here("data/landcover/landcover_geo.tif")
cropland <- here("data/landcover/cropland_binary.tif")
eq_string <- "A==1"
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   landcover, cropland, eq_string, options_calc)
system(command)

cropland_ratio <- here("data/landcover/cropland_ratio.tif")
te <- ext(bios); tr <- res(bios)
command <- sprintf(
    paste0("gdalwarp -ot Float32 -te %s %s %s %s -tr %s %s %s %s %s -r average"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, cropland, cropland_ratio)
system(command)

# - Dense tree/forest
treecover <- here("data/landcover/treecover_binary.tif")
eq_string <- "A==2"
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   landcover, treecover, eq_string, options_calc)
system(command)

treecover_ratio <- here("data/landcover/treecover_ratio.tif")
te <- ext(bios); tr <- res(bios)
command <- sprintf(
    paste0("gdalwarp -ot Float32 -te %s %s %s %s -tr %s %s %s %s %s -r average"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, treecover, treecover_ratio)
system(command)

# - Savanna (Shrub, grassland, and wetland)
savanna <- here("data/landcover/savanna_binary.tif")
eq_string <- "logical_or(logical_or(A==3, A==4), A==8)"
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   landcover, savanna, eq_string, options_calc)
system(command)

savanna_ratio <- here("data/landcover/savanna_ratio.tif")
te <- ext(bios); tr <- res(bios)
command <- sprintf(
    paste0("gdalwarp -ot Float32 -te %s %s %s %s -tr %s %s %s %s %s -r average"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, savanna, savanna_ratio)
system(command)

# - Water
water <- here("data/landcover/water_binary.tif")
eq_string <- "A==5"
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   landcover, water, eq_string, options_calc)
system(command)

water_ratio <- here("data/landcover/water_ratio.tif")
te <- ext(bios); tr <- res(bios)
command <- sprintf(
    paste0("gdalwarp -ot Float32 -te %s %s %s %s -tr %s %s %s %s %s -r average"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, water, water_ratio)
system(command)

# Settlement density
settlement <- here('data/landcover/settlement_binary.tif')
te <- ext(rast(landcover)); tr <- res(rast(landcover))
command <- sprintf(
    paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
           "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    here('data/open_building/buildings_ob.geojson'), settlement)
system(command)

settlement_ratio <- here("data/landcover/settlement_ratio.tif")
te <- ext(bios); tr <- res(bios)
command <- sprintf(
    paste0("gdalwarp -ot Float32 -te %s %s %s %s -tr %s %s %s %s %s -r average ",
           "-srcnodata None"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, settlement, settlement_ratio)
system(command)

fnames <- list.files("data/landcover", pattern = "_ratio.tif", full.names = TRUE)
landcovers <- rast(fnames)
writeRaster(
    landcovers, file.path(dst_dir, "crop_savanna_settlement_tree_water_ratio.tif"))

# Distance to settlement/big roads/waterbodies/rivers/streams
## Prepare the vectors from OpenStreetMap before start.
fnames <- list.files("data/osm", full.names = TRUE)
fnames <- c(fnames, file.path("data/open_building", "buildings_ob.geojson"))

te <- ext(bios); tr <- res(bios)
for (fn in fnames){
    fn_rst <- gsub(".geojson", ".tif", fn)
    
    # Rasterize
    command <- sprintf(
        paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
               "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
        te[1], te[3], te[2], te[4], tr[1], tr[2], fn, fn_rst)
    system(command)
    
    # Distance
    fn_dst <- gsub(".geojson", "_distance.tif", fn)
    command <- sprintf("gdal_proximity.py %s %s -co compress=lzw",
                       fn_rst, fn_dst)
    system(command)
}

distances <- rast(gsub(".geojson", "_distance.tif", fnames))
writeRaster(
    distances, 
    file.path(dst_dir, "distance_to_road_river_stream_waterbody_house.tif"))

# Topographic features: Terrain roughness and slope
dem_raw <- here("data/elevation/elevation.tif")
dem <- here("data/elevation/elevation_1km.tif")
te <- ext(bios); tr <- res(bios)
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s %s %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, dem_raw, dem)
system(command)

item <- "roughness"
fn <- here(sprintf("data/elevation/%s_1km.tif", item))
command <- sprintf("gdaldem %s %s %s -compute_edges -co compress=lzw",
                   item, dem, fn)
system(command)
terrain <- rast(c(dem, fn))
writeRaster(terrain, file.path(dst_dir, "elevation_roughness.tif"))

ndvis <- resample(ndvis, bios[[1]])
variables <- c(bios, ndvis, landcovers, distances, terrain)
variables <- mask(variables, bry)
writeRaster(variables, file.path(dst_dir, "variables.tif"))
