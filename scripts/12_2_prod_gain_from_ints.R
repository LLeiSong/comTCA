## --------------------------------------------
## Script name: prod_gain_from_ints
## Purpose of script: allocate cropland for future expansion based on
## different criteria.
## Author: Lei Song
## Date Created: 2023-10-22
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------
## 5 crops count for about 90% of food production
## Compare scenarios:
## crop production: only yield, weighted by current cropland area
## biodiversity
## carbon
## connectivity
## hybrid

## then:
## -50% yield gap close
## ---------------------------------------------

# Load libraries
library(terra)
library(parallel)
library(pbapply)

# Set directories
tdf_dir <- "data/tradeoff"

# Mask out the protected areas
pas <- rast(file.path(tdf_dir, "protected_areas.tif"))
pas[pas == 1] <- NA
farmable_area <- file.path(tdf_dir, "farmable_perc.tif") %>% 
    rast() %>% mask(pas) * 100
farmable_area[farmable_area <= 0] <- NA
cropland <- file.path(tdf_dir, "cropland_perc.tif") %>% 
    rast() %>% mask(pas) * 100
# Only these units will be evaluated

# Gather all inputs and standardize them
yields <- rast(file.path(tdf_dir, "agro_current_yield.tif")) %>% 
    mask(farmable_area)

costs <- file.path(
    tdf_dir, c("bio_index.tif", "carbon_density.tif", 
               "conn_index.tif", "agro_travel_time.tif")) %>% 
    rast() %>% mask(farmable_area)
weights <- lapply(costs, function(lyr){
    lyrs <- 1 - stretch(lyr / yields, minv = 0, maxv = 1)
    names(lyrs) <- names(yields)
    lyrs
})

names(costs) <- c("Biodiversity", "Carbon", "Connectivity", "Distance")

yields_as_weights <- 1 - stretch(1 / yields, minv = 0, maxv = 1)
weights <- c(list(yields_as_weights), weights)
names(weights) <- c("Yield", "Biodiversity", "Carbon", "Connectivity", "Distance")

global_weights <- file.path(tdf_dir, "cropland_perc.tif") %>% 
    rast() %>% mask(farmable_area)
names(global_weights) <- c("Cropland coverage")

# Get the production increase from land expansion for each pixel
# intes_gain = (atn_yield - current_yield) * current_land
atn_yields <- rast(file.path(tdf_dir, "agro_attainable_yield.tif")) %>% 
    mask(farmable_area)
atn_yields_150 <- rast(file.path(tdf_dir, "agro_attainable_yield_150.tif")) %>% 
    mask(farmable_area)

# Start the simulations---------------------------------------------------
## production gain from agricultural intensification ---------------------
## weight of planting area for maize, rice, sorghum, cassava and beans are
## 6, 2, 1, 1, 1 (11 in total)
vals_current <- values(c(cropland, yields)) %>% na.omit() %>% data.frame()
vals_atn <- values(c(cropland, atn_yields)) %>% na.omit() %>% data.frame()

# Make the crops before the run
## Current crops, area / 2 * 0.6
current_area <- ceiling(nrow(vals_current) * 0.43 * 0.64)
nums <- floor(current_area * c(6, 2, 1, 1, 1) / 11)
nums[5] <- nums[5] + (current_area - sum(nums))
crops_current <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
                   rep(4, nums[4]), rep(5, nums[5]))

## Attainable crops, area * 0.6, make full use of all land
atn_area <- ceiling(nrow(vals_current) * 0.6)
nums <- floor(atn_area * c(6, 2, 1, 1, 1) / 11)
nums[5] <- nums[5] + (atn_area - sum(nums))
crops_atn <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
               rep(4, nums[4]), rep(5, nums[5]))

prod_gain_from_intense <- do.call(rbind, pblapply(1:100, function(n){
    # Current yield
    set.seed(n)
    simu_plant_current <- vals_current %>% 
        sample_n(size = current_area) %>% 
        mutate(crop = sample(crops_current) + 1)
    
    simu_prod_current <- mclapply(1:nrow(simu_plant_current), function(i){
        simu_plant_current[i, 1] * simu_plant_current[i, simu_plant_current[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    # attainable yield
    set.seed(123 + n)
    simu_plant_atn <- vals_atn %>% 
        sample_n(size = atn_area) %>% 
        mutate(crop = sample(crops_atn) + 1)
    
    simu_prod_atn <- mclapply(1:nrow(simu_plant_atn), function(i){
        simu_plant_atn[i, 1] * simu_plant_atn[i, simu_plant_atn[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    data.frame(Current_production = simu_prod_current,
               Attainable_production = simu_prod_atn,
               Production_gain = simu_prod_atn - simu_prod_current)
}))
write.csv(prod_gain_from_intense, 
          file.path(tdf_dir, "production_gain_from_intensifify.csv"), 
          row.names = TRUE)

# More cassava
vals_current <- values(c(cropland, yields)) %>% na.omit() %>% data.frame()
vals_atn <- values(c(cropland, atn_yields)) %>% na.omit() %>% data.frame()

# Make the crops before the run
## Current crops, area / 2 * 0.6
current_area <- ceiling(nrow(vals_current) * 0.43 * 0.64)
nums <- floor(current_area * c(6, 2, 1, 1, 1) / 11)
nums[5] <- nums[5] + (current_area - sum(nums))
crops_current <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
                   rep(4, nums[4]), rep(5, nums[5]))

## Attainable crops, area * 0.6, make full use of all land
atn_area <- ceiling(nrow(vals_current) * 0.6)
nums <- floor(atn_area * c(5.4, 1.8, 0.9, 2, 0.9) / 11)
nums[5] <- nums[5] + (atn_area - sum(nums))
crops_atn <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
               rep(4, nums[4]), rep(5, nums[5]))

prod_gain_from_intense_mcsv <- do.call(rbind, pblapply(1:100, function(n){
    # Current yield
    set.seed(n)
    simu_plant_current <- vals_current %>% 
        sample_n(size = current_area) %>% 
        mutate(crop = sample(crops_current) + 1)
    
    simu_prod_current <- mclapply(1:nrow(simu_plant_current), function(i){
        simu_plant_current[i, 1] * simu_plant_current[i, simu_plant_current[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    # attainable yield
    set.seed(123 + n)
    simu_plant_atn <- vals_atn %>% 
        sample_n(size = atn_area) %>% 
        mutate(crop = sample(crops_atn) + 1)
    
    simu_prod_atn <- mclapply(1:nrow(simu_plant_atn), function(i){
        simu_plant_atn[i, 1] * simu_plant_atn[i, simu_plant_atn[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    data.frame(Current_production = simu_prod_current,
               Attainable_production = simu_prod_atn,
               Production_gain = simu_prod_atn - simu_prod_current)
}))
write.csv(prod_gain_from_intense_mcsv, 
          file.path(tdf_dir, "production_gain_from_intensifify_more_cassava.csv"), 
          row.names = TRUE)

# 150% yield gap closed
vals_current <- values(c(cropland, yields)) %>% na.omit() %>% data.frame()
vals_atn <- values(c(cropland, atn_yields_150)) %>% na.omit() %>% data.frame()

# Make the crops before the run
## Current crops, area / 2 * 0.6
current_area <- ceiling(nrow(vals_current) * 0.43 * 0.64)
nums <- floor(current_area * c(6, 2, 1, 1, 1) / 11)
nums[5] <- nums[5] + (current_area - sum(nums))
crops_current <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
                   rep(4, nums[4]), rep(5, nums[5]))

## Attainable crops, area * 0.6, make full use of all land
atn_area <- ceiling(nrow(vals_current) * 0.6)
nums <- floor(atn_area * c(6, 2, 1, 1, 1) / 11)
nums[5] <- nums[5] + (atn_area - sum(nums))
crops_atn <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
               rep(4, nums[4]), rep(5, nums[5]))

prod_gain_from_intense_150 <- do.call(rbind, pblapply(1:100, function(n){
    # Current yield
    set.seed(n)
    simu_plant_current <- vals_current %>% 
        sample_n(size = current_area) %>% 
        mutate(crop = sample(crops_current) + 1)
    
    simu_prod_current <- mclapply(1:nrow(simu_plant_current), function(i){
        simu_plant_current[i, 1] * simu_plant_current[i, simu_plant_current[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    # attainable yield
    set.seed(123 + n)
    simu_plant_atn <- vals_atn %>% 
        sample_n(size = atn_area) %>% 
        mutate(crop = sample(crops_atn) + 1)
    
    simu_prod_atn <- mclapply(1:nrow(simu_plant_atn), function(i){
        simu_plant_atn[i, 1] * simu_plant_atn[i, simu_plant_atn[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    data.frame(Current_production = simu_prod_current,
               Attainable_production = simu_prod_atn,
               Production_gain = simu_prod_atn - simu_prod_current)
}))
write.csv(prod_gain_from_intense_150, 
          file.path(tdf_dir, "production_gain_from_intensifify_150.csv"), 
          row.names = TRUE)

# 150% yield gap close and more cassava
vals_current <- values(c(cropland, yields)) %>% na.omit() %>% data.frame()
vals_atn <- values(c(cropland, atn_yields_150)) %>% na.omit() %>% data.frame()

# Make the crops before the run
## Current crops, area / 2 * 0.6
current_area <- ceiling(nrow(vals_current) * 0.43 * 0.64)
nums <- floor(current_area * c(6, 2, 1, 1, 1) / 11)
nums[5] <- nums[5] + (current_area - sum(nums))
crops_current <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
                   rep(4, nums[4]), rep(5, nums[5]))

## Attainable crops, area * 0.6, make full use of all land
atn_area <- ceiling(nrow(vals_current) * 0.6)
nums <- floor(atn_area * c(5.4, 1.8, 0.9, 2, 0.9) / 11)
nums[5] <- nums[5] + (atn_area - sum(nums))
crops_atn <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
               rep(4, nums[4]), rep(5, nums[5]))

prod_gain_from_intense_150_mcsv <- do.call(rbind, pblapply(1:100, function(n){
    # Current yield
    set.seed(n)
    simu_plant_current <- vals_current %>% 
        sample_n(size = current_area) %>% 
        mutate(crop = sample(crops_current) + 1)
    
    simu_prod_current <- mclapply(1:nrow(simu_plant_current), function(i){
        simu_plant_current[i, 1] * simu_plant_current[i, simu_plant_current[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    # attainable yield
    set.seed(123 + n)
    simu_plant_atn <- vals_atn %>% 
        sample_n(size = atn_area) %>% 
        mutate(crop = sample(crops_atn) + 1)
    
    simu_prod_atn <- mclapply(1:nrow(simu_plant_atn), function(i){
        simu_plant_atn[i, 1] * simu_plant_atn[i, simu_plant_atn[i, 7]]
    }, mc.cores = 12) %>% unlist() %>% sum()
    
    data.frame(Current_production = simu_prod_current,
               Attainable_production = simu_prod_atn,
               Production_gain = simu_prod_atn - simu_prod_current)
}))
write.csv(prod_gain_from_intense_150_mcsv, 
          file.path(tdf_dir, "production_gain_from_intensifify_150_more_cassave.csv"), 
          row.names = TRUE)
