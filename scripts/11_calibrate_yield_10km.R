## --------------------------------------------
## Script name: calibrate_yield
## Purpose of script: calibrate the SPAM 2017 gridded yield based on
## NSCA survey database
## Author: Lei Song
## Date Created: 2023-10-11
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------
## Notes: Data: only focus on rainfed.
## https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FSSKBW
## https://www.nbs.go.tz/index.php/en/
## 2021 for Tanzania
## maize: 1.6 tons/ha
## rice: 2.81 tons/ha
## sorghum: 1.04 tons/ha
## beans: 1.3 tons/ha
## cassava: 6.21 tons/ha
## --------------------------------------------

# Load libraries
library(terra)
library(sf)
library(dplyr)
library(readxl)
library(stringr)

# Set directories
data_dir <- "data/agriculture"
spam_ha_dir <- file.path(data_dir, 'SPAM', "spam2017v2r1_ssa_harv_area.geotiff")
spam_yield_dir <- file.path(data_dir, 'SPAM', "spam2017v2r1_ssa_yield.geotiff")

# Read relevant data
## Country boundary
bry <- st_read("data/geoms/mainland_tanzania.geojson") %>% select()

## Selected crops
crops <- c("maize", "paddy", "sorghum", "cassava", "beans")

## Level-1 districts
regions <- st_read("data/geoms/tanzania_level1.geojson") %>% 
    mutate(Region = ifelse(ADM1_EN == "Dar-es-salaam", 
                           "Dar es Salaam", ADM1_EN)) %>% select(Region)

## NSCA agricultural statistics
nsca_yield <- lapply(crops, function(name){
    read_excel(
        file.path(data_dir, "NSCA/region_yield_5crops.xlsx"), 
        sheet = name) %>% 
        mutate(Crop = name)}) %>% bind_rows()

# Match the level-1 admin polygons with NSCA
## Check the region names
if (!identical(intersect(unique(regions$Region), sort(unique(nsca_yield$Region))), 
              sort(unique(nsca_yield$Region)))) {
    message("Not identical region names found!")}

# Subset regions and zones
regions <- regions %>% filter(Region %in% nsca_yield$Region)

# Pre-process gridded harvest area and yield layers
## name conversion
names_cov <- data.frame(
    name = crops, 
    gyga_name = c("Rainfed maize", "Rainfed rice", 
                  "Rainfed sorghum", NA, "Rainfed common bean"),
    code = c("MAIZ", "RICE", "SORG", "CASS", "BEAN"))
practice_cov <- data.frame(
    name = c("All technologies together", "Rainfed high inputs", "Irrigation",
             "Rainfed low inputs", "Rainfed subsistence", "Rainfed"), 
    code = c("A", "H", "I", "L", "R", "S"))

## Yield
yields <- lapply(names_cov$code, function(crp) {
    fnames <- file.path(
        spam_yield_dir, sprintf("spam2017V2r1_SSA_Y_%s_%s.tif", 
                                crp, practice_cov$code))
    rast(fnames) %>% crop(bry) %>% mask(bry)})
names(yields) <- names_cov$name

# Calibrate yield
yields_calib <- do.call(c, lapply(names_cov$name, function(crp){
    # Only use rainfed
    yields_crp <- yields[[crp]]
    # kg/ha to tons/ha
    yield_R_crp <- subset(yields_crp, 5) / 1000
    
    ## Add factor to pixel values to shift the mean to match with surveys.
    do.call(merge, lapply(unique(regions$Region), function(rg){
        ply <- regions %>% filter(Region == rg)
        yield_R_rg <- crop(yield_R_crp, ply) %>% mask(ply)
        
        # Step 2: get real yield
        real_yield <- nsca_yield %>% filter(Region == rg) %>% 
            filter(Season == "both") %>% 
            filter(Crop == crp) %>% pull(Yield)
        
        # use the mean of non-zero values in [0.05, 0.95] to get the factor
        # more robust to outliers
        vals <- values(yield_R_rg)
        vals <- vals[!(is.na(vals))]
        vals <- vals[!is.nan(vals)]
        vals <- vals[vals > 0]
        vals <- vals[vals >= quantile(vals, 0.05) & vals <= quantile(vals, 0.95)]
        fact <- real_yield / mean(vals)
        
        # apply the factor
        yield_R_rg * fact
    }))
}))

# Save out
writeRaster(yields_calib, file.path(data_dir, "SPAM/yield_calibrated_5crops_tz.tif"))
