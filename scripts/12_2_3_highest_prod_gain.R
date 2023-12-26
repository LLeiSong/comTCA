## --------------------------------------------
## Script name: highest_prod_gain
## Purpose of script: highest possible production gain.
## Author: Lei Song
## Date Created: 2023-10-22
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------
## 5 crops count for about 72% of food production
## Scenarios to compare:
## - Different land usage on existing cropland
##      - More yield gap close (100 - 140 %)
##      - More cassava to plant (1 - 5 folds)
## This is just a run for the best results that
## do not consider other crops and practical issues.
## So in reality, this cannot be reached. This is just
## give us a reference.
## ---------------------------------------------

# Load libraries
library(terra)
library(dplyr)
library(parallel)
library(pbmcapply)

# Set directories
tdf_dir <- "data/tradeoff"
dst_dir <- "results/highest"
if(!dir.exists(dst_dir)) dir.create(dst_dir)

# Define the function
land_intes <- function(tdf_dir, dst_dir, yield=100, cassava=1, land_usage=1.0){
    # yield: one in [100, 110, 120, 130, 140, 150]
    # cassava: one in [1, 2, 3, 4]
    # land_usage: percentage equal or higher than 0.64
    # Mask out the protected areas
    pas <- rast(file.path(tdf_dir, "protected_areas.tif"))
    pas[pas == 1] <- NA
    
    # Get specific cell sizes, the same for all layers
    cellsizes <- cellSize(pas) / 1e6 * 100
    cropland <- file.path(tdf_dir, "cropland_perc.tif") %>% 
        rast() %>% mask(pas) * cellsizes
    cropland[cropland <= 0.9] <- NA
    farmable_area <- file.path(tdf_dir, "farmable_area.tif") %>% 
        rast() %>% mask(pas)
    farmable_area[farmable_area <= 0.9] <- NA
    # Only these units will be evaluated
    
    # Gather all inputs and standardize them
    yields <- rast(file.path(tdf_dir, "agro_current_yield.tif")) %>% 
        mask(cropland)
    
    # Get the production increase from land expansion for each pixel
    # intes_gain = (atn_yield - current_yield) * current_land
    if (yield == 100){
        atn_yields <- rast(file.path(tdf_dir, "agro_attainable_yield.tif")) %>% 
            mask(cropland)
        atn_yields_exp <- rast(file.path(tdf_dir, "agro_attainable_yield.tif")) %>% 
            mask(farmable_area)
    } else {
        atn_yields <- rast(
            file.path(tdf_dir, sprintf("agro_attainable_yield_%s.tif", yield))) %>% 
            mask(cropland)
        atn_yields_exp <- rast(
            file.path(tdf_dir, sprintf("agro_attainable_yield_%s.tif", yield))) %>% 
            mask(farmable_area)
    }
    
    # Start the simulations---------------------------------------------------
    ## production gain from agricultural intensification ---------------------
    ## weight of planting area for maize, rice, sorghum, cassava and beans are
    ## 6, 2, 1, 1, 1 (11 in total)
    vals_current <- values(c(cropland, yields)) %>% data.frame() %>% 
        mutate(cell_id = 1:nrow(.)) %>% na.omit() %>% 
        select(all_of(names(.)[c(7, 1, 2:6)]))
    vals_atn <- values(c(cropland, atn_yields)) %>% data.frame() %>% 
        mutate(cell_id = 1:nrow(.)) %>% na.omit() %>% 
        select(all_of(names(.)[c(7, 1, 2:6)]))
    vals_atn_exp <- values(c(farmable_area, atn_yields_exp)) %>% data.frame() %>% 
        mutate(cell_id = 1:nrow(.)) %>% na.omit() %>% 
        select(all_of(names(.)[c(7, 1, 2:6)]))
    
    # Convert to planted area to harvested area
    cnt_ratio <- c(0.88, 0.88, 0.87, 0.33, 0.86)
    names(cnt_ratio) <- c("maize", "paddy", "sorghum", "cassava", "beans")
    
    # Make the crops before the run
    ## Current crops, area * 0.653 * 0.73
    ## 0.653: planted_area_survey(non_tree_crops) / area_cropland (11765077/18389898) adjust 0.64 to 0.653
    ## 0.72: planted_area_survey(5 crops) / planted_area_survey(non_tree_crops) (8458379 / 11765077)
    ## The reason to adjust them to match with current production.
    current_area <- ceiling(nrow(vals_current) * land_usage * 0.82)
    current_nums <- floor(current_area * c(58, 20, 6, 6, 10) / 100)
    current_nums[5] <- current_nums[5] + (current_area - sum(current_nums))
    names(current_nums) <- names(cnt_ratio)
    
    ## Attainable crops, area * 0.72 with different percentage of usage
    atn_area <- ceiling(nrow(vals_atn) * land_usage * 0.78)
    atn_area_exp <- ceiling(nrow(vals_atn_exp) * land_usage * 0.65)
    
    if (cassava == 1){
        atn_nums <- floor(atn_area * c(58, 20, 6, 6, 10) / 100)
        atn_nums_exp <- floor(atn_area_exp * c(58, 20, 6, 6, 10) / 100)
    } else if (cassava == 2){
        atn_nums <- floor(atn_area * c(54.3, 18.7, 5.6, 12, 9.4) / 100)
        atn_nums_exp <- floor(atn_area_exp * c(54.3, 18.7, 5.6, 12, 9.4) / 100)
    } else if (cassava == 3){
        atn_nums <- floor(atn_area * c(50.6, 17.5, 5.2, 18, 8.7) / 100)
        atn_nums_exp <- floor(atn_area_exp * c(50.6, 17.5, 5.2, 18, 8.7) / 100)
    } else if (cassava == 4){
        atn_nums <- floor(atn_area * c(46.9, 16.2, 4.8, 24, 8.1) / 100)
        atn_nums_exp <- floor(atn_area_exp * c(46.9, 16.2, 4.8, 24, 8.1) / 100)
    } else if (cassava == 5){
        atn_nums <- floor(atn_area * c(43.2, 14.9, 4.5, 30, 7.4) / 100)
        atn_nums_exp <- floor(atn_area_exp * c(43.2, 14.9, 4.5, 30, 7.4) / 100)
    } else {
        stop("No such value.")
    }
    atn_nums[5] <- atn_nums[5] + (atn_area - sum(atn_nums))
    names(atn_nums) <- names(cnt_ratio)
    
    atn_nums_exp[5] <- atn_nums_exp[5] + (atn_area_exp - sum(atn_nums_exp))
    names(atn_nums_exp) <- names(cnt_ratio)
    
    # Make global relative benefits for crops
    # which means applying the convert percentage
    for (crp in names(cnt_ratio)){
        vals_current[[crp]] <- vals_current[[crp]] * cnt_ratio[[crp]]
        vals_atn[[crp]] <- vals_atn[[crp]] * cnt_ratio[[crp]]
        vals_atn_exp[[crp]] <- vals_atn_exp[[crp]] * cnt_ratio[[crp]]
    }
    
    ################## Current condition #####################
    set.seed(123)
    round_num <- 1
    while(TRUE){
        # Make benefit matrix
        if (round_num == 1){
            current_mats <- vals_current %>% 
                mutate(max_ind = max.col(vals_current[, -c(1:2)])) %>% 
                mutate(max_ind = names(cnt_ratio)[max_ind])
        } else {
            current_mats <- vals_current %>% 
                select(-c(unique(current_best_crops$crop))) 
            current_mats <- current_mats %>% 
                mutate(max_ind = max.col(current_mats[, -c(1:2)])) %>% 
                mutate(max_ind = names(current_mats)[-c(1:2)][max_ind]) %>% 
                filter(!cell_id %in% current_best_crops$cell_id)
        }
        
        # Summarize the best pixel for each crop
        current_best <- current_mats %>% 
            group_by(max_ind) %>% summarise(n = n())
        
        # Get the highest production for current situation
        ## Go first for the crops with sufficient best pixels
        current_best <- current_best %>% 
            mutate(n_need = current_nums[current_best$max_ind]) %>% 
            filter(n >= n_need)
        
        best_crops_ss <- do.call(
            rbind, lapply(unique(current_best$max_ind), function(crp){
                num_to_select <- current_best %>% filter(max_ind == crp) %>% 
                    pull(n_need)
                current_mats %>% filter(max_ind == crp) %>% 
                    select(all_of(c("cell_id", "cropland_ratio", crp))) %>% 
                    rename(yield = all_of(crp)) %>% 
                    mutate(prod = cropland_ratio * yield) %>% 
                    arrange(-prod) %>% 
                    slice(1:num_to_select) %>% 
                    mutate(crop = crp)
            }))
        
        # Accumulate the results
        if (round_num == 1){
            current_best_crops <- best_crops_ss
        } else {
            current_best_crops <- rbind(current_best_crops, best_crops_ss)
        }
        
        round_num <- round_num + 1
        
        # Stop if finish
        if (all(names(current_nums) %in% unique(current_best_crops$crop))){
            break
        }
    }
    
    ################## Attainable condition #####################
    set.seed(123)
    round_num <- 1
    while(TRUE){
        # Make benefit matrix
        if (round_num == 1){
            atn_mats <- vals_atn %>% 
                mutate(max_ind = max.col(vals_atn[, -c(1:2)])) %>% 
                mutate(max_ind = names(cnt_ratio)[max_ind])
        } else {
            atn_mats <- vals_atn %>% 
                select(-c(unique(atn_best_crops$crop))) 
            atn_mats <- atn_mats %>% 
                mutate(max_ind = max.col(atn_mats[, -c(1:2)])) %>% 
                mutate(max_ind = names(atn_mats)[-c(1:2)][max_ind]) %>% 
                filter(!cell_id %in% atn_best_crops$cell_id)
        }
        
        # Summarize the best pixel for each crop
        atn_best <- atn_mats %>% 
            group_by(max_ind) %>% summarise(n = n())
        
        # Get the highest production for current situation
        ## Go first for the crops with sufficient best pixels
        atn_best <- atn_best %>% 
            mutate(n_need = atn_nums[atn_best$max_ind]) %>% 
            filter(n >= n_need)
        
        best_crops_ss <- do.call(
            rbind, lapply(unique(atn_best$max_ind), function(crp){
                num_to_select <- atn_best %>% filter(max_ind == crp) %>% 
                    pull(n_need)
                atn_mats %>% filter(max_ind == crp) %>% 
                    select(all_of(c("cell_id", "cropland_ratio", crp))) %>% 
                    rename(yield = all_of(crp)) %>% 
                    mutate(prod = cropland_ratio * yield) %>% 
                    arrange(-prod) %>% 
                    slice(1:num_to_select) %>% 
                    mutate(crop = crp)
            }))
        
        # Accumulate the results
        if (round_num == 1){
            atn_best_crops <- best_crops_ss
        } else {
            atn_best_crops <- rbind(atn_best_crops, best_crops_ss)
        }
        
        round_num <- round_num + 1
        
        # Stop if finish
        if (all(names(atn_nums) %in% unique(atn_best_crops$crop))){
            break
        }
    }
    
    ################## Full expand condition #####################
    set.seed(123)
    round_num <- 1
    while(TRUE){
        # Make benefit matrix
        if (round_num == 1){
            atn_exp_mats <- vals_atn_exp %>% 
                mutate(max_ind = max.col(vals_atn_exp[, -c(1:2)])) %>% 
                mutate(max_ind = names(cnt_ratio)[max_ind])
        } else {
            atn_exp_mats <- vals_atn_exp %>% 
                select(-c(unique(atn_exp_best_crops$crop))) 
            atn_exp_mats <- atn_exp_mats %>% 
                mutate(max_ind = max.col(atn_exp_mats[, -c(1:2)])) %>% 
                mutate(max_ind = names(atn_exp_mats)[-c(1:2)][max_ind]) %>% 
                filter(!cell_id %in% atn_exp_best_crops$cell_id)
        }
        
        # Summarize the best pixel for each crop
        atn_best <- atn_exp_mats %>% 
            group_by(max_ind) %>% summarise(n = n())
        
        # Get the highest production for current situation
        ## Go first for the crops with sufficient best pixels
        atn_best <- atn_best %>% 
            mutate(n_need = atn_nums_exp[atn_best$max_ind]) %>% 
            filter(n >= n_need)
        
        best_crops_ss <- do.call(
            rbind, lapply(unique(atn_best$max_ind), function(crp){
                num_to_select <- atn_best %>% filter(max_ind == crp) %>% 
                    pull(n_need)
                atn_exp_mats %>% filter(max_ind == crp) %>% 
                    select(all_of(c("cell_id", "farmable_ratio", crp))) %>% 
                    rename(yield = all_of(crp)) %>% 
                    mutate(prod = farmable_ratio * yield) %>% 
                    arrange(-prod) %>% 
                    slice(1:num_to_select) %>% 
                    mutate(crop = crp)
            }))
        
        # Accumulate the results
        if (round_num == 1){
            atn_exp_best_crops <- best_crops_ss
        } else {
            atn_exp_best_crops <- rbind(atn_exp_best_crops, best_crops_ss)
        }
        
        round_num <- round_num + 1
        
        # Stop if finish
        if (all(names(atn_nums) %in% unique(atn_exp_best_crops$crop))){
            break
        }
    }
    
    best_prod <- data.frame(
        highest_current_production = sum(current_best_crops$prod),
        highest_attainable_production = sum(atn_best_crops$prod),
        highest_expansion_production = sum(atn_exp_best_crops$prod))
    
    message(sprintf("Scenario: %s, %s, %s.", land_usage, yield, cassava))
    print(best_prod)
    
    fname <- file.path(
        dst_dir, sprintf("highest_prod_gain_CASS%s_Y%s_USG%s.csv", 
                         cassava, yield, round(land_usage, 2) * 100))
    write.csv(best_prod, fname, row.names = TRUE)
    message(sprintf("Save file to %s.", fname))
}

# Yields with intensification and expansion
for (land_usage in c(0.64, 0.8, 1.0)){
    for (yield in seq(100, 140, 10)){
        land_intes(tdf_dir, dst_dir, yield, 1, land_usage)
    }
    
    # Cassava with same 
    for (cassava in 2:5){
        land_intes(tdf_dir, dst_dir, 100, cassava, land_usage)
    }
}
