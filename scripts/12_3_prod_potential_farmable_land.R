## --------------------------------------------
## Script name: prod_potential_farmable_land
## Purpose of script: production potential from farmable area. 
## Almost the same as prod_gain_from_ints
## Author: Lei Song
## Date Created: 2023-10-22
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------
## 5 crops count for about 80% of food production
## Scenarios to compare:
## - Different land usage on existing cropland
##      - More yield gap close (100 - 150 %)
##      - More cassava to plant (1 - 5 folds)
## ---------------------------------------------

# Load libraries
library(terra)
library(dplyr)
library(parallel)
library(pbmcapply)

# Set directories
tdf_dir <- "data/tradeoff"
dst_dir <- "results/expansion"
if(!dir.exists(dst_dir)) dir.create(dst_dir)

# Define the function
land_potential <- function(tdf_dir, dst_dir, yield=100, cassava=1, land_usage=1.0){
    # yield: one in [100, 110, 120, 130, 140, 150]
    # cassava: one in [1, 2, 3, 4]
    # land_usage: percentage equal or higher than 0.653
    # Mask out the protected areas
    pas <- rast(file.path(tdf_dir, "protected_areas.tif"))
    pas[pas == 1] <- NA
    
    # Get specific cell sizes, the same for all layers
    cellsizes <- cellSize(pas) / 1e6 * 100
    farmable_area <- file.path(tdf_dir, "farmable_perc.tif") %>% 
        rast() %>% mask(pas) * cellsizes
    farmable_area[farmable_area <= 0] <- NA
    # Only these units will be evaluated
    
    # Gather all inputs and standardize them
    yields <- rast(file.path(tdf_dir, "agro_current_yield.tif")) %>% 
        mask(farmable_area)
    
    # Get the production increase from land expansion for each pixel
    # intes_gain = (atn_yield - current_yield) * current_land
    if (yield == 100){
        atn_yields <- rast(file.path(tdf_dir, "agro_attainable_yield.tif")) %>% 
            mask(farmable_area)
    } else {
        atn_yields <- rast(
            file.path(tdf_dir, sprintf("agro_attainable_yield_%s.tif", yield))) %>% 
            mask(farmable_area)
    }
    
    # Start the simulations---------------------------------------------------
    ## production gain from agricultural intensification ---------------------
    ## weight of planting area for maize, rice, sorghum, cassava and beans are
    vals_atn <- values(c(farmable_area, atn_yields)) %>% na.omit() %>% data.frame()
    
    # Convert to planted area to harvested area
    cnt_ratio <- c(0.88, 0.88, 0.87, 0.33, 0.86)
    names(cnt_ratio) <- c("maize", "rice", "sorghum", "cassava", "beans")
    
    # Make the crops before the run
    ## Attainable crops, area * 0.72 with different percentage of usage
    atn_area <- ceiling(nrow(vals_atn) * land_usage * 0.72)
    
    if (cassava == 1){
        nums <- floor(atn_area * c(58, 20, 6, 6, 10) / 100)
    } else if (cassava == 2){
        nums <- floor(atn_area * c(54.3, 18.7, 5.6, 12, 9.4) / 100)
    } else if (cassava == 3){
        nums <- floor(atn_area * c(50.6, 17.5, 5.2, 18, 8.7) / 100)
    } else if (cassava == 4){
        nums <- floor(atn_area * c(46.9, 16.2, 4.8, 24, 8.1) / 100)
    } else if (cassava == 5){
        nums <- floor(atn_area * c(43.2, 14.9, 4.5, 30, 7.4) / 100)
    } else {
        stop("No such value.")
    }
    nums[5] <- nums[5] + (atn_area - sum(nums))
    crops_atn <- c(rep(1, nums[1]), rep(2, nums[2]), rep(3, nums[3]), 
                   rep(4, nums[4]), rep(5, nums[5]))
    
    # Do the calculation
    prod_gain_from_intense <- do.call(rbind, pbmclapply(1:100, function(n){
        # attainable yield
        set.seed(123 + n)
        simu_plant_atn <- vals_atn %>% 
            sample_n(size = atn_area) %>% 
            mutate(crop = sample(crops_atn) + 1) %>% 
            mutate(cnt_ratio = cnt_ratio[crop - 1])
        
        simu_prod_atn <- lapply(1:nrow(simu_plant_atn), function(i){
            simu_plant_atn[i, 1] * simu_plant_atn[i, simu_plant_atn[i, 7]] * 
                simu_plant_atn[i, 8]
        }) %>% unlist() %>% sum()
        
        data.frame(Attainable_production = simu_prod_atn)
    }, mc.cores = 12, ignore.interactive = TRUE))
    fname <- file.path(
        dst_dir, sprintf("prod_gain_CASS%s_Y%s_USG%s.csv", 
                         cassava, yield, round(land_usage, 2) * 100))
    write.csv(prod_gain_from_intense, fname, row.names = TRUE)
    message(sprintf("Save file to %s.", fname))
}

# Yields with full land
for (land_usage in c(0.653, 0.8, 1.0)){
    for (yield in seq(100, 140, 10)){
        land_potential(tdf_dir, dst_dir, yield, 1, land_usage)
    }
    
    # Cassava with same 
    for (cassava in 2:5){
        land_potential(tdf_dir, dst_dir, 100, cassava, land_usage)
    }
}
