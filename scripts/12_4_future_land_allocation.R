## --------------------------------------------
## Script name: land_alloc_scenarios
## Purpose of script: function to allocate the
## land for reconciling
## Author: Lei Song
## Date Created: 2023-10-28
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## Compare scenarios:
## yield, cbetas = "1,0,0,0,0"
## biodiversity, cbetas = "0,1,0,0,0"
## carbon, cbetas = "0,0,1,0,0"
## connectivity, cbetas = "0,0,0,1,0"
## cost distance, cbetas = "0,0,0,0,1"
## yield + cost distance (only benefits), cbetas = "0.5,0,0,0,0.5"
## biodiversity + carbon + connectivity (only costs), cbetas = "0,0.33,0.33,0.33,0"
## yield + biodiversity + carbon + cost distance (no connectivity), cbetas = "0.25,0.25,0.25,0,0.25"
## yield + biodiversity + carbon + connectivity + cost distance (everything), cbetas = "0.2,0.2,0.2,0.2,0.2"
## --------------------------------------------

# Load libraries
library(terra)
library(parallel)
library(dplyr)
library(stringr)
library(optparse)

# Function to allocate the land for agriculture
land_allocate <- function(
    # Inputs
    atn_yields, # the attainable yield to calculate production gain
    farmable_area, # all new farmable area
    weights, # weight for other factors except yields
    costs, # the cost by expansion
    # Scenario setting
    name = "Y", # the name for the experiment.
    cbetas, # the prefer factor for each weight, the first one is for yield
    scenario = "Y100", # Y for increased yield gap close, and CASS for increase cassava plant area.
    land_usage = 0.653, # the percentage of annual land use. 0.65, 0.8, or 1.0
    production_need = 0, # production target (e.g. double)
    # Common setting
    spatial = TRUE,
    seed = 123,
    dst_dir){
    
    # Create directories
    num_dir <- file.path(dst_dir, "numbers")
    if (!dir.exists(num_dir)) dir.create(num_dir)
    spatial_dir <- file.path(dst_dir, "spatial")
    if (!dir.exists(spatial_dir)) dir.create(spatial_dir)
    
    # Re-calculate weights 
    calib_weights <- lapply(names(cbetas), function(name){
        weights[[name]] * cbetas[name]})
    names(calib_weights) <- names(weights)
    
    # Get the production gain from expansion
    prod_gain_exp <- atn_yields * farmable_area
    
    # Get global order of the pixels
    glb_weight <- do.call(c, lapply(names(calib_weights$Yield), function(crp){
        # Merge all weights
        weights_crp <- do.call(
            c, lapply(1:length(calib_weights), 
                      function(n) calib_weights[[n]][[crp]])) %>% sum()
    })) %>% max()
    
    # The workflow is like:
    ## First, do an initial run to get a not-bad but not the best result
    ## Starting from there, to converge to the best result
    ## If the calibrated production is higher, then delete the redundant pixels
    ## If the production is lower, add new pixels and continue the converge
    
    ########################## Initial run ##########################
    ## Get decision table for initial run
    ## This initial time, it doesn't matter if there is duplicates for the best
    ## pixels because the pixels are selected one by one.
    decision_mats <- lapply(names(calib_weights$Yield), function(crp){
        # Merge all weights
        weights_crp <- do.call(
            c, lapply(1:length(calib_weights), 
                      function(n) calib_weights[[n]][[crp]])) %>% sum()
        
        vals <- values(c(weights_crp, glb_weight, prod_gain_exp[[crp]], 
                         costs, farmable_area)) %>% 
            data.frame() %>% mutate(cell_id = 1:nrow(.)) %>% 
            na.omit() %>% 
            rename(weight = sum, Production_gain = all_of(crp)) %>% 
            mutate(grp = weight >= max)
        rbind(vals %>% filter(grp == TRUE) %>% 
                  arrange(desc(weight), desc(farmable_ratio)),
              vals %>% filter(grp == FALSE) %>% 
                  mutate(diff = max - weight) %>% 
                  arrange(diff, desc(weight), farmable_ratio) %>% 
                  select(-diff)) %>% select(-c(max))
    })
    names(decision_mats) <- names(calib_weights$Yield)
    
    # Start the allocate
    ## 1. Create a vector to store the selected cell ids. This is important because
    ## cells will not be selected up to down
    accum_cost <- data.frame(
        matrix(ncol = 10, nrow = 0, 
               dimnames=list(NULL, c("cell_id", "Production_gain", 
                                     "Biodiversity", "Carbon", 
                                     "Connectivity", "Distance", 
                                     "farmable_ratio", "weight", "crop", "id"))))
    ## 2. Select pixel by pixel
    set.seed(seed)
    prod_to_catch <- production_need
    pixel_num <- nrow(decision_mats[[1]])
    i <- 1
    
    # Convert to planted area to harvested area
    cnt_ratio <- c(0.88, 0.88, 0.87, 0.33, 0.86)
    names(cnt_ratio) <- c("maize", "paddy", "sorghum", "cassava", "beans")
    
    # Get the ratio of each crop
    cassava <- ifelse(
        str_detect(scenario, "CASS"), str_extract(scenario, "[0-9]+"), 1)
    if (cassava == 1){
        ratios <- c(58, 20, 6, 6, 10) * 10
    } else if (cassava == 2){
        ratios <- c(54.3, 18.7, 5.6, 12, 9.4) * 10
    } else if (cassava == 3){
        ratios <- c(50.6, 17.5, 5.2, 18, 8.7) * 10
    } else if (cassava == 4){
        ratios <- c(46.9, 16.2, 4.8, 24, 8.1) * 10
    } else if (cassava == 5){
        ratios <- c(43.2, 14.9, 4.5, 30, 7.4) * 10
    }
    
    # Do the loop
    while (prod_to_catch > 0) {
        ## Randomly select a crop to grow
        crop_id <- sample(
            c(rep("maize", ratios[1]), rep("paddy", ratios[2]),
              rep("sorghum", ratios[3]), rep("cassava", ratios[4]),
              rep("beans", ratios[5])), 1)
        decision_mat <- decision_mats[[crop_id]]
        
        if (i %% 5000 == 0){
            for (crp in c("maize", "paddy", "sorghum", "cassava", "beans")){
                decision_mats[[crp]] <- decision_mats[[crp]] %>% 
                    filter(!cell_id %in% accum_cost$cell_id)}
            num <- round((production_need - prod_to_catch) / production_need * 100, 2)
            message(sprintf("Finish %s%%", num))
        }
        
        # Remove all selected cell_ids and select the first one
        decision_mat <- decision_mat %>% 
            filter(!cell_id %in% accum_cost$cell_id) %>% slice(1)
        
        # Record the selected row
        accum_cost[i, ] <- decision_mat %>% 
            mutate(crop = crop_id, id = i) %>% 
            select(all_of(names(accum_cost)))
        
        # Move on
        i <- i + 1
        # Apply land usage percentage and the ratio of planted area to 
        # harvested area to get the reasonable production
        prod_to_catch <- prod_to_catch - 
            decision_mat$Production_gain * land_usage * cnt_ratio[crop_id]
    }
    
    # Save out
    fname <- file.path(
        num_dir, sprintf("%s_%s_USG%s_simu_%s.csv", 
                         name, scenario, round(land_usage * 100, 0), seed))
    write.csv(accum_cost, fname, row.names = FALSE)
    message("Finish the initial run. Head to converge.")
    
    ########################## Start the converge ##########################
    # Rebuild the original decision_mats because of the deduction for speeding up
    decision_mats <- lapply(names(calib_weights$Yield), function(crp){
        # Merge all weights
        weights_crp <- do.call(
            c, lapply(1:length(calib_weights), 
                      function(n) calib_weights[[n]][[crp]])) %>% sum()
        
        vals <- values(c(weights_crp, glb_weight, prod_gain_exp[[crp]], 
                         costs, farmable_area)) %>% 
            data.frame() %>% mutate(cell_id = 1:nrow(.)) %>% 
            na.omit() %>% 
            rename(weight = sum, Production_gain = all_of(crp)) %>% 
            mutate(grp = weight >= max)
        rbind(vals %>% filter(grp == TRUE) %>% 
                  arrange(desc(weight), desc(farmable_ratio)),
              vals %>% filter(grp == FALSE) %>% 
                  mutate(diff = max - weight) %>% 
                  arrange(diff, desc(weight), farmable_ratio) %>% 
                  select(-diff)) %>% select(-c(max))
    })
    names(decision_mats) <- names(calib_weights$Yield)
    
    # Remove duplicates
    for (i in 2:length(decision_mats)){
        crp <- names(calib_weights$Yield)[i]
        other_crps <- names(calib_weights$Yield)[1:(i - 1)]
        ids_occup <- unlist(sapply(decision_mats[other_crps], function(x) {
            x %>% filter(grp == TRUE) %>% pull(cell_id)}))
        decision_mats[[i]] <- decision_mats[[i]] %>% 
            mutate(grp = ifelse(cell_id %in% ids_occup, FALSE, grp))
    }
    
    # Start the converging process
    round_num <- 1
    while(TRUE){
        # Record the nrow before the process
        nrow_before <- nrow(accum_cost)
        
        ## Get the initial crop series
        nums_crops <- accum_cost %>% select(crop) %>% 
            group_by(crop) %>% summarise(n = n())
        best_nums_crops <- do.call(rbind, lapply(
            decision_mats, function(tb) tb %>% filter(grp == TRUE) %>% nrow())) %>% 
            data.frame(best_n = .) %>% mutate(crop = row.names(.))
        nums_crops <- left_join(nums_crops, best_nums_crops, by = "crop")
        
        ## Select the best pixels for the self-sufficient crops
        crops_ss <- nums_crops %>% filter(best_n >= n) %>% pull(crop)
        decision_mats_ss <- decision_mats[crops_ss]
        for (i in 1:length(crops_ss)){
            crp <- crops_ss[i]
            if (i == 1) {
                decision_mats_ss[[i]] <- decision_mats_ss[[i]] %>% 
                    slice(1:(nums_crops %>% filter(crop == crp) %>% pull(n)))
            } else {
                decision_mats_ss[[i]] <- decision_mats_ss[[i]] %>% 
                    filter(!cell_id %in% decision_mats_ss[[i - 1]]$cell_id) %>% 
                    slice(1:(nums_crops %>% filter(crop == crp) %>% pull(n)))
            }
        }; names(decision_mats_ss) <- crops_ss
        
        ## Make a record of global best pixels
        best_pixels_crops <- do.call(rbind, lapply(
            crops_ss, function(crp) decision_mats[[crp]] %>% 
                filter(grp == TRUE) %>% mutate(crop = crp))) %>% 
            mutate(benefit = weight * farmable_ratio)
        
        # Now reassign non-self-sufficient crops
        # Use a simple method to speed up
        # WARNING: if the allocation will happen to most of the farmable area,
        # the code may need to be changed to select one by one.
        # I will keep this for now.
        while (TRUE){
            # Selected pixels not to touch
            cell_ids_notouch <- unlist(sapply(decision_mats_ss, function(x) x$cell_id))
            
            # Make weight matrix for non-self-sufficient crops
            crops_not_ss <- setdiff(names(calib_weights[[1]]), crops_ss)
            glb_weight <- do.call(c, lapply(crops_not_ss, function(crp){
                # Merge all weights
                weights_crp <- do.call(
                    c, lapply(1:length(calib_weights), 
                              function(n) calib_weights[[n]][[crp]])) %>% sum()
            })) %>% max()
            
            decision_mats_not_ss <- lapply(crops_not_ss, function(crp){
                # Merge all weights
                weights_crp <- do.call(
                    c, lapply(1:length(calib_weights), 
                              function(n) calib_weights[[n]][[crp]])) %>% sum()
                
                vals <- values(c(weights_crp, glb_weight, prod_gain_exp[[crp]], 
                                 costs, farmable_area)) %>% 
                    data.frame() %>% mutate(cell_id = 1:nrow(.)) %>% 
                    na.omit() %>% 
                    rename(weight = sum, Production_gain = all_of(crp)) %>% 
                    filter(!cell_id %in% cell_ids_notouch) %>% 
                    mutate(grp = weight >= max)
                rbind(vals %>% filter(grp == TRUE) %>% 
                          arrange(desc(weight), desc(farmable_ratio)),
                      vals %>% filter(grp == FALSE) %>% 
                          mutate(diff = max - weight) %>% 
                          arrange(diff, desc(weight), farmable_ratio) %>% 
                          select(-diff)) %>% select(-c(max, grp))
            }); names(decision_mats_not_ss) <- crops_not_ss
            
            # Remove potential duplicates
            for (i in 1:length(crops_not_ss)){
                crp <- crops_not_ss[i]
                if (i == 1) {
                    decision_mats_not_ss[[i]] <- decision_mats_not_ss[[i]] %>% 
                        slice(1:(nums_crops %>% filter(crop == crp) %>% pull(n)))
                } else {
                    decision_mats_not_ss[[i]] <- decision_mats_not_ss[[i]] %>% 
                        filter(!cell_id %in% decision_mats_not_ss[[i - 1]]$cell_id) %>% 
                        slice(1:(nums_crops %>% filter(crop == crp) %>% pull(n)))
                }
            }
            
            # Check if any selected pixels could have higher benefit for ss crops
            benefit_selected <- best_pixels_crops %>% 
                filter(cell_id %in% cell_ids_notouch) %>% arrange(benefit)
            
            ids_to_check <- unlist(sapply(decision_mats_not_ss, function(x) x$cell_id))
            ids_to_check <- ids_to_check[ids_to_check %in% best_pixels_crops$cell_id]
            benefit_abaonded <- best_pixels_crops %>% 
                filter(cell_id %in% ids_to_check) %>% filter(crop %in% crops_ss) %>% 
                arrange(-benefit)
            
            benefit_abaonded <- do.call(
                rbind, lapply(unique(benefit_selected$crop), function(ccrp){
                    ids_has_chance <- c()
                    ids_to_replace <- c()
                    benefit_abaonded_crp <- benefit_abaonded %>% filter(crop == ccrp)
                    benefit_selected_crp <- benefit_selected %>% filter(crop == ccrp)
                    if (nrow(benefit_abaonded_crp) > 0){
                        for (i in 1:nrow(benefit_abaonded_crp)){
                            row_i <- benefit_abaonded_crp[i, ]
                            if (any(row_i$benefit > benefit_selected_crp$benefit)){
                                # Replace one by one
                                ids_has_chance <- c(ids_has_chance, row_i$cell_id)
                                ids_to_replace <- c(ids_to_replace, benefit_selected_crp[1, 'cell_id'])
                                benefit_selected_crp <- benefit_selected_crp %>% slice(-1)
                            }
                        }
                    }
                    
                    benefit_abaonded_crp %>% filter(cell_id %in% ids_has_chance)
                }))
            
            message(sprintf("%s rows to re-allocate.", nrow(benefit_abaonded)))
            
            if (nrow(benefit_abaonded) > 0){
                # Replace these pixels within ss crops
                for (ccrp in unique(benefit_abaonded$crop)){
                    rows_from <- benefit_abaonded %>% filter(crop == ccrp)
                    decision_mats_ss[[ccrp]] <- 
                        rbind(decision_mats_ss[[ccrp]] %>% 
                                  mutate(benefit = weight * farmable_ratio) %>% 
                                  arrange(benefit) %>% slice(-c(1:nrow(rows_from))),
                              rows_from %>% select(-crop)) %>% 
                        arrange(-benefit) %>% select(-benefit)
                }
            } else{
                # All clear, move on
                break
            }
        }
        
        # Now replace these new rows with the old ones
        clmns <- setdiff(names(accum_cost), c("crop", "id"))
        for(crp in crops_ss){
            accum_cost[accum_cost$crop == crp, clmns] <- 
                decision_mats_ss[[crp]] %>% select(all_of(clmns))
        }
        
        for(crp in crops_not_ss){
            accum_cost[accum_cost$crop == crp, clmns] <- 
                decision_mats_not_ss[[crp]] %>% select(all_of(clmns))
        }
        
        # It's the time to calculate how much production for the calibrated plan
        accum_cost <- accum_cost %>% group_by(crop) %>% 
            arrange(-weight, .by_group = TRUE) %>% ungroup()
        accum_cost <- accum_cost %>% 
            mutate(Production_gain_actual = Production_gain * land_usage * cnt_ratio[crop])
        production_actual <- sum(accum_cost$Production_gain_actual)
        
        # More than original run, deduct
        if (production_actual >= production_need){
            message(sprintf("Production: %s. Extra pixels need to be removed.", production_actual))
            
            set.seed(seed)
            prod_to_catch <- production_actual
            while (TRUE){
                ## Randomly select a crop to remove
                crop_id <- sample(
                    c(rep("maize", ratios[1]), rep("paddy", ratios[2]),
                      rep("sorghum", ratios[3]), rep("cassava", ratios[4]),
                      rep("beans", ratios[5])), 1)
                
                row_to_remove <- accum_cost %>% filter(crop == crop_id) %>% 
                    slice(nrow(.))
                
                prod_to_catch <- prod_to_catch - row_to_remove$Production_gain_actual
                
                if (prod_to_catch < production_need){
                    # Work is done.
                    break
                } else {
                    accum_cost <- accum_cost %>% filter(cell_id != row_to_remove$cell_id)
                    
                    # Check if this pixel can be better for other crop
                    benefits_from_row <- do.call(
                        rbind, lapply(setdiff(names(decision_mats), crop_id), 
                                      function(crp){
                        this_w <- decision_mats[[crp]] %>% 
                            filter(cell_id == row_to_remove$cell_id)
                        accum_cost %>% filter(crop == crp) %>% 
                            mutate(weight_gain = this_w$weight - weight) %>% 
                            arrange(-weight_gain) %>% slice(1)
                    })) %>% arrange(-weight_gain) %>% slice(1)
                    
                    # Replace this row
                    if(benefits_from_row$weight_gain > 0){
                        # Rebuild the row
                        row_from <- decision_mats[[benefits_from_row$crop]] %>% 
                            filter(cell_id == row_to_remove$cell_id) %>% 
                            select(all_of(names(accum_cost)[1:8])) %>% 
                            mutate(crop = benefits_from_row$crop, id = benefits_from_row$id,
                                   Production_gain_actual = Production_gain * land_usage * 
                                       cnt_ratio[benefits_from_row$crop])
                        accum_cost <- accum_cost %>% filter(cell_id != benefits_from_row$cell_id)
                        accum_cost <- rbind(accum_cost, row_from)
                        accum_cost <- accum_cost %>% group_by(crop) %>% 
                            arrange(-weight, .by_group = TRUE) %>% ungroup()
                        
                        prod_to_catch <- prod_to_catch + row_from$Production_gain_actual - 
                            benefits_from_row$Production_gain_actual
                    }
                }
            }
            
            # Remove the extra column and re-index.
            accum_cost <- accum_cost %>% select(-Production_gain_actual) %>% 
                mutate(id = 1:nrow(.))
            
            # Post the news
            message(sprintf("Finish round %s for the converge.", round_num))
            message(sprintf("Land area: %s", sum(accum_cost$farmable_ratio)))
            
            # Save out the updated file
            write.csv(accum_cost, fname, row.names = FALSE)
            
            # Stop here
            break
        # Oh, no. More pixels are in need.
        } else {
            message(sprintf("Production: %s. Extra pixels need to be added.", production_actual))
            
            # Remove the extra column because no need it and re-index
            accum_cost <- accum_cost %>% select(-Production_gain_actual) %>% 
                mutate(id = 1:nrow(.))
            
            set.seed(seed)
            prod_to_catch <- production_need - production_actual
            i_cvg <- nrow(accum_cost) + 1
            
            while (prod_to_catch > 0) {
                ## Randomly select a crop to grow
                crop_id <- sample(
                    c(rep("maize", ratios[1]), rep("paddy", ratios[2]),
                      rep("sorghum", ratios[3]), rep("cassava", ratios[4]),
                      rep("beans", ratios[5])), 1)
                decision_mat <- decision_mats[[crop_id]]
                
                # Remove all selected cell_ids and select the first one
                decision_mat <- decision_mat %>% 
                    filter(!cell_id %in% accum_cost$cell_id) %>% slice(1)
                
                # Record the selected row
                accum_cost[i_cvg, ] <- decision_mat %>% 
                    mutate(crop = crop_id, id = i_cvg) %>% 
                    select(all_of(names(accum_cost)))
                
                # Move on
                i_cvg <- i_cvg + 1
                # Apply land usage percentage and the ratio of planted area to 
                # harvested area to get the reasonable production
                prod_to_catch <- prod_to_catch - 
                    decision_mat$Production_gain * land_usage * cnt_ratio[crop_id]
            }
        }
        
        # Post the news
        message(sprintf("Finish round %s for the converge.", round_num))
        message(sprintf("Land area: %s", sum(accum_cost$farmable_ratio)))
        
        # Record the nrow after the process
        nrow_after <- nrow(accum_cost)
        
        if (abs(nrow_before - nrow_after) <= 1){ # avoid endless loop.
            break
        } else {
            # Save out the updated file
            write.csv(accum_cost, fname, row.names = FALSE)
            round_num <- round_num + 1
        }
    }
    
    message("Finish the converge.")
    
    ########################## Save out the spatial ##########################
    if (spatial){
        # Assign the result to spatial map
        agro_land <- farmable_area # as a template
        names(agro_land) <- "expanded_farmland"
        new_vals <- rep(NA, ncell(agro_land))
        crps <- c("maize", "paddy", "sorghum", 'cassava', "beans")
        new_vals[accum_cost$cell_id] <- factor(accum_cost$crop, levels = crps, labels = crps)
        values(agro_land) <- new_vals
        levels(agro_land) <- data.frame(id = 1:5, crop = crps)
        dst_fname <- file.path(
            spatial_dir, sprintf("%s_%s_USG%s_simu_%s.tif", 
                                 name, scenario, round(land_usage * 100, 0), seed))
        writeRaster(agro_land, dst_fname, datatype = "INT1U", overwrite = TRUE)
    }
}

# Set directories
tdf_dir <- "data/tradeoff"
intes_dir <- "results/intensification"
dst_dir <- "results/scenarios"

# Mask out the protected areas
pas <- rast(file.path(tdf_dir, "protected_areas.tif"))
pas[pas == 1] <- NA

# Convert percentage to area in ha
farmable_area <- rast(file.path(tdf_dir, "farmable_area.tif"))
# Only these units will be evaluated

# Command line inputs
option_list <- list(
    make_option(c("-s", "--seed"), 
                action = "store", default = 1, type = 'integer',
                help = "The seed for the iteration [default %default]"),
    make_option(c("-c", "--cbetas"), 
                action = "store", type = 'character',
                help = paste0("The weights for yield, biodiversity, ",
                              "carbon, connectivity and distance in order.")),
    make_option(c("-o", "--scenario"), 
                action = "store", type = 'character',
                help = paste0("The scenario for the simulation, ",
                              "[Y100, Y110, Y120, Y130, Y140, ",
                              "CASS2, CASS3, CASS4, CASS5].")),
    make_option(c("-l", "--land_usage"), 
                action = "store", type = 'numeric',
                help = "The percentage of land usage, [0.653, 0.8, 1.0]."),
    make_option(c("-d", "--dst_dir"), 
                action = "store", type = 'character',
                help = "The sub-path in results directory to save."))
opt <- parse_args(OptionParser(option_list = option_list))
seed <- opt$seed
cbetas <- opt$cbetas
cbetas <- as.numeric(strsplit(cbetas, ",")[[1]])
scenario <- opt$scenario
land_usage <- as.numeric(opt$land_usage)
yield_num <- ifelse(
    str_detect(scenario, "Y"), str_extract(scenario, "[0-9]+"), 100)
cassava <- ifelse(
    str_detect(scenario, "CASS"), str_extract(scenario, "[0-9]+"), 1)

# Make directory for the results
dst_dir <- file.path(dst_dir, opt$dst_dir)
if (!dir.exists(dst_dir)) dir.create(dst_dir)

# Get the right attainable yield layer
atn_yield_fname <- file.path(
    tdf_dir, sprintf("agro_attainable_yield%s.tif", 
                     ifelse(yield_num == 100, "", paste0("_", yield_num))))
atn_yields <- rast(atn_yield_fname) %>% mask(farmable_area)
message(sprintf("Read attainable layer: %s.", basename(atn_yield_fname)))

# Gather all weight inputs and standardize them
cnt_ratio <- c(0.88, 0.88, 0.87, 0.33, 0.86)
names(cnt_ratio) <- c("maize", "paddy", "sorghum", "cassava", "beans")
atn_yields_calib <- do.call(c, lapply(names(atn_yields), function(crp){
    atn_yields[[crp]] * cnt_ratio[[crp]]}))

# Gather all weight inputs and standardize them
## Ecological costs
costs <- file.path(
    tdf_dir, c("bio_index.tif", "carbon_density.tif", 
               "conn_index.tif", "agro_travel_time.tif")) %>% 
    rast() %>% mask(farmable_area)
names(costs) <- c("Biodiversity", "Carbon", "Connectivity", "Distance")

## Get efficiency-based weights for ecological costs
weights <- lapply(costs, function(lyr){
    lyrs <- lyr / atn_yields_calib
    names(lyrs) <- names(atn_yields_calib)
    lyrs
})

## Get yield weights
## factors to prefer to the smallest value.
yields_as_weights <- 1 / atn_yields_calib

## Put them together
weights <- c(list(yields_as_weights), weights)
names(weights) <- c("Yield", "Biodiversity", "Carbon", "Connectivity", "Distance")

## Here to normalize the layers across all layers for each element
weights <- lapply(weights, function(wt_lyrs){
    glb_min <- min(global(wt_lyrs, "min", na.rm = TRUE)[, 1])
    glb_max <- max(global(wt_lyrs, "max", na.rm = TRUE)[, 1])
    
    1 - (wt_lyrs - glb_min) / c(glb_max - glb_min)
})

# Get the production increase on existing cropland
fname <- file.path(
    intes_dir, sprintf("prod_gain_CASS%s_Y%s_USG%s.csv", 
                       cassava, yield_num, round(land_usage * 100, 0)))
prod_intes <- read.csv(fname)
message(sprintf("Read production increase for: %s.", basename(fname)))

production_need <- mean(prod_intes$Current_production) * 2 - 
    mean(prod_intes$Attainable_production)
message(sprintf("Get production remaining to get: %s.", 
                round(production_need, 2)))

# Set the factors
cbetas <- c("Yield" = cbetas[1], "Biodiversity" = cbetas[2], "Carbon" = cbetas[3], 
            "Connectivity" = cbetas[4], "Distance" = cbetas[5])
exp_name <- paste(c("Y", "B", "Ca", "Co", "D")[cbetas != 0], collapse = "")

# Start the run
if (production_need > 0){
    land_allocate(
        atn_yields, farmable_area, weights, costs, 
        exp_name, cbetas, scenario, land_usage, production_need, 
        TRUE, seed, dst_dir)
} else {
    message("No more land is needed.")
}
