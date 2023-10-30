## --------------------------------------------
## Script name: gather_results
## Purpose of script: Gather the results of simulations.
## Author: Lei Song
## Date Created: 2023-10-22
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------

## Compare different scenarios
f_dir <- "/scratch/lsong36/comTCA/data/scenarios"
items <- c("Y", "B", "Ca", "Co", "D", "hybrid")

for (item in items) {
    message(item)
    fnames <- list.files(f_dir, pattern = sprintf("^%s", item), full.names = TRUE)
    
    simu_list <- do.call(rbind, mclapply(fnames, function(fname){
        test <- read.csv(fname)
        # Calculate volume or use weight
        test$Biodiversity <- test$Biodiversity * test$farmable_ratio
        test$Carbon <- test$Carbon * test$farmable_ratio
        test$Connectivity <- test$Connectivity * test$farmable_ratio
        test$Distance <- test$Distance * test$farmable_ratio
        
        # accumulate them
        test$Production_gain <- cumsum(test$Production_gain)
        test$Biodiversity <- cumsum(test$Biodiversity)
        test$Carbon <- cumsum(test$Carbon)
        test$Connectivity <- cumsum(test$Connectivity)
        test$Distance <- cumsum(test$Distance)
        test$farmable_ratio <- cumsum(test$farmable_ratio)
        
        test %>% select('Production_gain', 'Biodiversity', 'Carbon', 
                        'Connectivity',  'Distance', 'farmable_ratio', 'id')
    }, mc.cores = 10))
    
    simu_mean <- simu_list %>% group_by(id) %>% 
        summarise(across(everything(), mean))
    
    simu_sd <- simu_list %>% group_by(id) %>% 
        summarise(across(everything(), sd))
    
    write.csv(simu_mean, file.path(f_dir, sprintf("%s_simulation_mean.csv", item)),
              row.names = FALSE)
    write.csv(simu_sd, file.path(f_dir, sprintf("%s_simulation_sd.csv", item)),
              row.names = FALSE)
}

# Compare different future farming practice
f_dir <- "/scratch/lsong36/comTCA/data/compare"
items <- c("HY", "CS", 'RG')

for (item in items) {
    message(item)
    fnames <- list.files(f_dir, pattern = sprintf("_%s_", item), full.names = TRUE)
    
    simu_list <- do.call(rbind, mclapply(fnames, function(fname){
        test <- read.csv(fname)
        # Calculate volume or use weight
        test$Biodiversity <- test$Biodiversity * test$farmable_ratio
        test$Carbon <- test$Carbon * test$farmable_ratio
        test$Connectivity <- test$Connectivity * test$farmable_ratio
        test$Distance <- test$Distance * test$farmable_ratio
        
        # accumulate them
        test$Production_gain <- cumsum(test$Production_gain)
        test$Biodiversity <- cumsum(test$Biodiversity)
        test$Carbon <- cumsum(test$Carbon)
        test$Connectivity <- cumsum(test$Connectivity)
        test$Distance <- cumsum(test$Distance)
        test$farmable_ratio <- cumsum(test$farmable_ratio)
        
        test %>% select('Production_gain', 'Biodiversity', 'Carbon', 
                        'Connectivity',  'Distance', 'farmable_ratio', 'id')
    }, mc.cores = 10))
    
    simu_mean <- simu_list %>% group_by(id) %>% 
        summarise(across(everything(), mean))
    
    simu_sd <- simu_list %>% group_by(id) %>% 
        summarise(across(everything(), sd))
    
    write.csv(simu_mean, file.path(f_dir, sprintf("%s_hybrid_mean.csv", item)),
              row.names = FALSE)
    write.csv(simu_sd, file.path(f_dir, sprintf("%s_hybrid_sd.csv", item)),
              row.names = FALSE)
}

# Spatial maps
f_dir <- "/scratch/lsong36/comTCA/data/spatial_allocation"
items <- c("HY", "CS", 'RG')

for (item in items) {
    message(item)
    fnames <- list.files(f_dir, pattern = sprintf("_%s_", item), full.names = TRUE)
    lyrs <- rast(fnames)
    
    # For crops
    crps <- do.call(c, lapply(1:5, function(crp){
        lyr <- lyrs
        lyr[lyr != crp] <- NA
        lyr[lyr == crp] <- 1
        sum(lyr, na.rm = TRUE)
    }))
    
    # Overall
    ovr <- lyrs
    ovr[ovr > 0] <- 1
    ovr <- sum(ovr, na.rm = TRUE)
    
    out <- c(ovr, crps)
    names(out) <- c('farmland', "maize", "paddy", "sorghum", "cassava", "beans")
    
    writeRaster(out, file.path(f_dir, sprintf("fut_agro_land_hybrid_%s.tif", item)))
}
