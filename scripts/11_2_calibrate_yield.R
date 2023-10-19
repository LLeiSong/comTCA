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
## cassava: 6.21 tons/ha
## beans: 1.3 tons/ha
## --------------------------------------------

# Load libraries
library(terra)
library(sf)
library(dplyr)
library(readxl)
library(stringr)

# Set directories
data_dir <- "data/agriculture"

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
    code = c("MAIZ", "RICE", "SORG", "CASS", "BEAN"))

## Yield
yields <- lapply(names_cov$code, function(crp) {
    fname <- file.path(data_dir, sprintf("yield_05_95_%s_R_1km.tif", crp))
    rast(fname) %>% crop(bry) %>% mask(bry)})
names(yields) <- names_cov$name

## Cropland mask
te <- ext(yields[[1]][[1]]); tr <- res(yields[[1]][[1]])
options_warp <- '-co compress=lzw'
command <- sprintf(
    paste0("gdalwarp -t_srs EPSG:4326 -te %s %s %s %s -tr %s %s %s %s %s"),
    te[1], te[3], te[2], te[4], tr[1], tr[2], 
    options_warp, file.path(data_dir, "cropland_tz.tif"), 
    file.path(data_dir, "cropland_tz_1km.tif"))
system(command)
farmland <- rast(file.path(data_dir, "cropland_tz_1km.tif"))

# Calibrate yield
yields_calib <- do.call(c, lapply(names_cov$name, function(crp){
    message(sprintf("Calibdate crop %s", crp))
    # Only use rainfed
    yields_crp <- yields[[crp]]
    
    ## Get the best match quantile value with surveys.
    do.call(merge, lapply(unique(regions$Region), function(rg){
        ply <- regions %>% filter(Region == rg)
        yield_rg <- crop(yields_crp, ply) %>% mask(ply)
        farmland_rg <- crop(farmland, ply) %>% mask(ply)
        
        # Step 2: get real yield
        real_yield <- nsca_yield %>% filter(Region == rg) %>% 
            filter(Season == "both") %>% 
            filter(Crop == crp) %>% pull(Yield)
        
        if (is.na(real_yield)){
            # get the median if no real yield report
            best_quantile <- 'prediction.quantile..0.5'
        } else {
            # Get the best match quantile
            yield_rg_msk <- mask(yield_rg, farmland_rg)
            best_quantile <- global(yield_rg_msk, mean, na.rm = TRUE) %>% 
                data.frame() %>% rename(yield = mean) %>% 
                mutate(quantile = row.names(.)) %>% 
                mutate(residual = abs(yield - real_yield)) %>% 
                filter(residual == min(residual)) %>% slice(1) %>% pull(quantile)
        }
        
        # Extract the best match layer
        subset(yield_rg, best_quantile)
    }))
}))
names(yields_calib) <- names_cov$name

# Check the statistics
global(yields_calib, mean, na.rm = TRUE)

# Save out
fname <- file.path(data_dir, "yield_calibrated_5crops_tz_1km.tif")
writeRaster(yields_calib, fname)

# Use 95% quantile as the attainable yield
yields_atn <- do.call(c, lapply(names_cov$name, function(crp){
    # Only use rainfed
    yields_crp <- yields[[crp]]
    subset(yields_crp, 'prediction.quantile..0.95')
}))
names(yields_atn) <- names_cov$name

# Check the statistics
global(yields_atn, mean, na.rm = TRUE)

# Save out
fname <- file.path(data_dir, "yield_attainable_5crops_tz_1km.tif")
writeRaster(yields_atn, fname)

# Check yield calibrated result
library(ggplot2)

# Collect data
yield_compare <- do.call(rbind, lapply(unique(regions$Region), function(rg){
    ply <- regions %>% filter(Region == rg)
    yield_rg <- crop(yields_calib, ply) %>% mask(ply)
    farmland_rg <- crop(farmland, ply) %>% mask(ply)
    
    # Step 2: get real yield
    real_yield <- nsca_yield %>% filter(Region == rg) %>% 
        filter(Season == "both") %>% select(Region, Crop, Yield)
    
    yield_rg_msk <- mask(yield_rg, farmland_rg)
    yield_cl <- global(yield_rg_msk, mean, na.rm = TRUE) %>% 
        data.frame() %>% mutate(Crop = row.names(.)) %>% 
        rename(Calibrated_Yield = mean) %>% mutate(Region = rg) %>% 
        inner_join(real_yield, by = c('Region', 'Crop'))
}))

# Plot
ggplot(yield_compare) + 
    geom_point(aes(x = Yield, y = Calibrated_Yield)) +
    facet_grid(~Crop, scales = "free") +
    theme_pubclean()
