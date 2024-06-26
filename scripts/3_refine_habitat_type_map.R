# Title     : Refine habitat type map
# Objective : Use the fine resolution land cover map to refine the existing
#           : habitat type map. Use gdal functions a lot to speed up.
# Created by: Lei Song
# Created on: 07/26/23
# NOTE      : Take advantage of the efforts for making the original map, 
#             and avoid duplicated work. Do not do any meaningless guess. 
#             The most obvious update is the arable land and urban areas.

# Load libraries
library(here)
library(terra)
library(rgrass)
library(dplyr)
library(sf)
library(tidyverse)

# Set path/dirs
habitat_fname <- file.path(
    here("data/habitats/"), 
    "iucn_habitatclassification_composite_lvl2_ver004.tif")
landcover_fname <- here("data/landcover/landcover.tif")
landcover_geo_fname <- here("data/landcover/landcover_geo.tif")

# Read maps
habitat_lvl2 <- rast(habitat_fname)

# First, re-project land cover map to Geographic projection
options_warp <- paste0('-multi -wo NUM_THREADS=ALL_CPUS ',
                       '-co TILED=YES -co compress=lzw -co BIGTIFF=IF_NEEDED')
command <- sprintf(
    paste0("gdalwarp -t_srs EPSG:4326 %s %s %s"),
    options_warp, landcover_fname, landcover_geo_fname)
system(command)

# Second, clip Tanzania off from the habitat type map
habitat_tz <- here("data/habitats/habitat_tz_lvl2.tif")
habitat_temp <- here("data/habitats/habitat_temp.tif")

## Resample
landcover <- rast(landcover_geo_fname)
te <- ext(landcover); tr <- res(landcover)
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s %s %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, habitat_fname, habitat_temp)
system(command)

## Mask
options_calc <- paste0('--NoDataValue=0 --hideNoData ',
                       '--co TILED=YES --co compress=lzw ',
                       '--co BIGTIFF=IF_NEEDED --co NUM_THREADS=ALL_CPUS')
eq_string <- "(A != 255) * B"
command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                   landcover_geo_fname, habitat_temp, habitat_tz, 
                   eq_string, options_calc)
system(command)
file.remove(habitat_temp)

# Then, get the habitat map outside of Tanzania
## Resample
te <- ext(habitat_lvl2); tr <- res(habitat_lvl2)
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s %s %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, habitat_tz, habitat_temp)
system(command)

## Mask
habitat_outside <- here("data/habitats/habitat_outside_lvl2.tif")
eq_string <- "(A == 0) * B"
command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                   habitat_temp, habitat_fname, habitat_outside, 
                   eq_string, options_calc)
system(command)
file.remove(habitat_temp); rm(habitat_temp)

# After these steps, we got two maps: habitat map of Tanzania and habitat map
# outside of Tanzania. Next step, we will update the habitat map of Tanzania
# (habitat_tz) using the new land cover map (landcover_geo_fname)
################################ Convert table ##############################
# built-up (done)
#   urban areas 1405
#   buffer to get rural garden (20m) 1404
# 
# Cropland (done)
#   Arable land 1401
# 
# wetland + water (done)
#   wetlands 500 - 518
#     
# Forest
#   plantations 1403
#   subtropical/tropical heavily degraded former forest 1406
#   other forest 100 - 109
# 
# shrub + grass + bareland + built-up
#   Savanna 200 - 202
#   Shrubland 300 - 308
#   grassland 400 - 407
#   pastureland 1402
# 
# Desert: should not have any desert in Tanzania

######################## Query Open Buildings dataset ###################
### Follow the tutorial to download Google Open Buildings dataset 
## https://sites.research.google/open-buildings/
bry <- st_read(here("data/geoms/tanzania_coarse_bry.geojson"))
ob_tile <- read_sf(here("data/geoms/tiles.geojson")) %>% 
    slice(unique(unlist(st_intersects(bry, .))))

dir.create(here('data/temp'), showWarnings = F)
options(timeout = 60000)
buildings_ob <- do.call(rbind, lapply(ob_tile$tile_url, function(ul){
    # Download
    dl_nm <- file.path(here('data/temp'), basename(ul))
    if (!file.exists(dl_nm)) download.file(ul, dl_nm)
    
    plys <- read_csv(dl_nm)
    chunks <- split(1:nrow(plys), ceiling(seq_along(1:nrow(plys)) / 10000))
    
    do.call(rbind, pbmclapply(chunks, function(ids){
        st_as_sf(plys[ids, ], wkt = "geometry", crs = 4326) %>% 
            slice(unique(unlist(st_intersects(bry, .)))) %>% 
            st_make_valid()
    }, mc.cores = detectCores()))
}))

dir.create(here('data/open_building'))
st_write(buildings_ob, here('data/open_building/buildings_ob.geojson'))
unlink(here("data/temp"), recursive = TRUE)

###################### Urban areas and Rural gardens #####################
landcover <- rast(landcover_geo_fname)
habitat <- rast(habitat_tz)

# Set up GRASS to use
if (!dir.exists(here("data/grass"))) dir.create(here("data/grass"))
gisBase <- "/opt/local/lib/grass82"
initGRASS(gisBase = gisBase,
          home = here("data/grass"),
          gisDbase = here("data/grass"),
          SG = landcover,
          mapset = "PERMANENT",
          location = "landcover",
          override = TRUE)

## Rasterize Open Buildings
te <- ext(habitat); tr <- res(habitat)
command <- sprintf(
    paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
           "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    here('data/open_building/buildings_ob.geojson'),
    here('data/landcover/buildings.tif'))
system(command)

## Add 5 meter buffer to include residual urban areas,
## and add additional 50 meter buffer as rural gardens.
execGRASS("r.in.gdal",
          flags = c("o", "overwrite"),
          input = here('data/landcover/buildings.tif'),
          band = 1,
          output = "buildings")
execGRASS("g.region", raster = "buildings")

## Make buffering
execGRASS('r.buffer', flags = c("overwrite"),
          parameters = list(input = 'buildings', 
                            output = 'buildings_5',
                            distances = 5))
execGRASS('r.buffer', flags = c("overwrite"),
          parameters = list(input = 'buildings_5', 
                            output = 'buildings_50',
                            distances = 10))

## Save out
execGRASS('r.out.gdal', flags = c("overwrite"),
          parameters = list(
              input = 'buildings_5', type = "Byte",
              createopt = "compress=lzw", 
              output = here('data/landcover/buildings_5_buf.tif')))
execGRASS('r.out.gdal', flags = c("overwrite"),
          parameters = list(
              input = 'buildings_50', type = "Byte",
              createopt = "compress=lzw", 
              output = here('data/landcover/buildings_plus_50_buf.tif')))

## Now put them together
options <- '--NoDataValue=0 --hideNoData --co compress=lzw'
eq_string <- sprintf("logical_or(A==1, A==2)")
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   here('data/landcover/buildings_5_buf.tif'), 
                   here('data/landcover/urban_areas_1405.tif'),
                   eq_string, options)
system(command)

eq_string <- sprintf("(A==2) * logical_or(B==1,B==6)")
command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                   here('data/landcover/buildings_plus_50_buf.tif'), 
                   landcover_geo_fname,
                   here('data/landcover/rural_gardens_1404.tif'),
                   eq_string, options)
system(command)

############################### Cropland ##############################
eq_string <- sprintf("A==1")
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   landcover_geo_fname, 
                   here('data/landcover/arable_land_1401.tif'),
                   eq_string, options)
system(command)

############################### Wetland + water #######################
# Use both landcover map and habitat map
eq_string <- paste0("logical_or(logical_and(A>=500, A<600), B==5) * ",
                    "logical_or(A<900, A>=1400)")
command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                   habitat_tz, landcover_geo_fname, 
                   here('data/landcover/water_wetland_0_500.tif'),
                   eq_string, options)
system(command)

################ Intermediate accumulate the maps #####################
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
eq_string <- paste0("(A==1) * 1401 + (A==0) * (B==1) * 1405 + ",
                    "(A==0) * (B==0) * (C==1) * 1404")
command <- sprintf('gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
                   here('data/landcover/arable_land_1401.tif'), 
                   here('data/landcover/urban_areas_1405.tif'),
                   here('data/landcover/rural_gardens_1404.tif'),
                   here('data/landcover/al_ua_rg_1401_1404_1405.tif'),
                   eq_string, options_int)
system(command)

eq_string <- paste0("logical_and(A>=500, A<600) * A + ",
                    "logical_or(A<500, A>=600) * ((B==1) * 500 + (B==0) * C)")
command <- sprintf('gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
                   habitat_tz, here('data/landcover/water_wetland_0_500.tif'),
                   here('data/landcover/al_ua_rg_1401_1404_1405.tif'),
                   here('data/landcover/wl_al_ua_rg_500_1401_1404_1405.tif'),
                   eq_string, options_int)
system(command)

######################### Forest relevant ##############################
eq_string <- paste0("logical_or(logical_and(A>=100, A<200), ",
                    "A==1406) * (B==2) * A")
command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                   habitat_tz, landcover_geo_fname,
                   here('data/landcover/100_1406.tif'),
                   eq_string, options_int)
system(command)

# Fill with the nearest class
forest_to_fill <- rast(here('data/landcover/100_1406.tif'))
forest_to_fill <- aggregate(forest_to_fill, fact = 20, fun = "modal", 
                         na.rm = TRUE, cores = 11) # back to 100m
for (i in 1:10){
    forest_to_fill <- focal(forest_to_fill, w = 3, fun = "modal", 
                            na.rm = TRUE, na.policy = "only")}

# save out
writeRaster(lcv_to_fill, here('data/landcover/100_1406_filled.tif'))

# Resample the lcv_to_fill and fill the nodata for the first round
template <- rast(landcover_geo_fname)
te <- ext(template); tr <- res(template); rm(template)
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s -co compress=lzw  %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    here('data/landcover/100_1406_filled.tif'), 
    here('data/landcover/100_1406_filled_resample.tif'))
system(command)

# Fill
eq_string <- "(A==2) * ((B!=0) * B + (A==2) * (B==0) * C)"
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
command <- sprintf(
    'gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
    landcover_geo_fname, 
    here('data/landcover/100_1406.tif'), 
    here('data/landcover/100_1406_filled_resample.tif'), 
    here('data/landcover/100_1406_full.tif'),
    eq_string, options_int)
system(command)

## The 1403: plantation class is not trustable at all.
## So update it based on https://doi.org/10.3390/rs13163081.
## NOTE: this mask can always been updated with better available product.

# Update plantation areas
habitat <- rast(habitat_tz)
te <- ext(habitat); tr <- res(habitat)
command <- sprintf(
    paste0("gdal_rasterize -te %s %s %s %s -tr %s %s ",
           "-a_nodata 0 -ot Byte -co compress=lzw -burn 1 %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    here('data/landcover/planatation_area.geojson'),
    here('data/landcover/planatation_area.tif'))
system(command)

eq_string <- paste0("(B==1) * (C==2) * 1403 + (C==2) * (B!=1) * A")
command <- sprintf('gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
                   here('data/landcover/100_1406_full.tif'), 
                   here('data/landcover/planatation_area.tif'),
                   landcover_geo_fname,
                   here('data/landcover/100_1406_1403.tif'),
                   eq_string, options_int)
system(command)

## Fill the rest pixels with the most common forest type: 105:Forest
eq_string <- "(A==2) * (B==0) * 105 + (A==2) * (B!=0) * B"
command <- sprintf(
    'gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
    landcover_geo_fname, 
    here('data/landcover/100_1406_1403.tif'), 
    here('data/landcover/100_1403_1406_full.tif'), 
    eq_string, options_int)
system(command)

################ Intermediate accumulate the maps #########################
eq_string <- "(A!=0) * A + (A==0) * B"
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                   here('data/landcover/wl_al_ua_rg_500_1401_1404_1405.tif'), 
                   here('data/landcover/100_1403_1406_full.tif'), 
                   here('data/landcover/100_500_1400_no_1402.tif'), 
                   eq_string, options_int)
system(command)

####################### Shrub + grass + bareland ##########################
# Respect the habitat type map as much as possible
# Make the landcover mask
eq_string <- "logical_or(logical_or(logical_or(logical_or(A==3,A==4),A==6),A==7),A==8)"
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
command <- sprintf('gdal_calc.py -A %s --outfile=%s --calc="%s" %s',
                   landcover_geo_fname, here('data/landcover/sgb_mask.tif'), 
                   eq_string, options_int)
system(command)

# Apply mask
eq_string <- paste0("logical_or(logical_and(A>=200, A<500), ",
                    "A==1402) * (B==1) * A")
command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                   habitat_tz, here('data/landcover/sgb_mask.tif'),
                   here('data/landcover/200_300_400_1402.tif'),
                   eq_string, options_int)
system(command)

lcv_to_fill <- rast(here('data/landcover/200_300_400_1402.tif'))
lcv_to_fill <- aggregate(lcv_to_fill, fact = 20, fun = "modal", 
                         na.rm = TRUE, cores = 11) # back to 100m

# Use wide neighbor majority to fill the other areas
lcv_to_fill <- focal(
    lcv_to_fill, w = 11, fun = "modal", na.policy = "only", na.rm = T)

# save out
writeRaster(lcv_to_fill, here('data/landcover/200_300_400_1402_filled.tif'))

# Resample the lcv_to_fill and fill the nodata
template <- rast(landcover_geo_fname)
te <- ext(template); tr <- res(template); rm(template)
command <- sprintf(
    paste0("gdalwarp -te %s %s %s %s -tr %s %s -co compress=lzw  %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    here('data/landcover/200_300_400_1402_filled.tif'), 
    here('data/landcover/200_300_400_1402_filled_resample.tif'))
system(command)

# Fill
eq_string <- "A * ((B!=0) * B + (B==0) * C)"
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
command <- sprintf(
    'gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
    here('data/landcover/sgb_mask.tif'), 
    here('data/landcover/200_300_400_1402.tif'), 
    here('data/landcover/200_300_400_1402_filled_resample.tif'), 
    here('data/landcover/200_300_400_1402_full.tif'),
    eq_string, options_int)
system(command)

## Because the majority is Savanna, so all remaining area will be signed as
## dry savanna.
eq_string <- "(D!=255) * ((A!=0) * A + (A==0) * ((C!=0) * C + (C==0) * 201))"
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
command <- sprintf(
    'gdal_calc.py -A %s -B %s -C %s -D %s --outfile=%s --calc="%s" %s',
    here('data/landcover/100_500_1400_no_1402.tif'), # go first
    here('data/landcover/sgb_mask.tif'), 
    here('data/landcover/200_300_400_1402_full.tif'), 
    landcover_geo_fname,
    here('data/landcover/habitat_tz_refined_messy_coastal.tif'), 
    eq_string, options_int)
system(command)

# To mask out the oceanic areas
eq_string <- "logical_and(B!=0, logical_or(B<900, B>=1400)) * A"
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
command <- sprintf(
    'gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
    here('data/landcover/habitat_tz_refined_messy_coastal.tif'),
    habitat_tz,
    here('data/habitats/habitat_tz_refined.tif'), 
    eq_string, options_int)
system(command)

# Finally, add the Rocky areas back, mainly peak of Meru, Kilimanjaro
eq_string <- "(B==600) * B + (B!=600) * A"
options_int <- '--NoDataValue=0 --hideNoData --co compress=lzw --type UInt16'
command <- sprintf(
    'gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
    here('data/habitats/habitat_tz_refined.tif'),
    habitat_tz, here('data/habitats/habitat_tz_refined_final.tif'), 
    eq_string, options_int)
system(command)

# The final habitat type map of Tanzania is fully refined and ready to be used.