## --------------------------------------------
## Script name: gather_results
## Purpose of script: Gather the results of simulations.
## Author: Lei Song
## Date Created: 2023-10-22
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------

## Library
library(parallel)
library(dplyr)
library(terra)

## Compare different scenarios
data_dir <- "results/USG80" # results/scenarios
usage <- 80
dst_dir <- file.path(data_dir, "summary")
if(!dir.exists(dst_dir)) dir.create(dst_dir)
f_dir <- file.path(data_dir, "numbers")
fnames <- list.files(f_dir, full.names = TRUE)
scenarios <- unique(sapply(strsplit(basename(fnames), "_"), function(x) x[1]))
changes <- unique(sapply(strsplit(basename(fnames), "_"), function(x) x[2]))

## Overall
smry_all_exps <- do.call(rbind, lapply(scenarios, function(item){
    do.call(rbind, lapply(changes, function(change){
        message("Start: ", item, " ", change)
        fnames <- list.files(f_dir, pattern = sprintf("^%s_%s_", item, change), 
                             full.names = TRUE)
        
        if (length(fnames) != 0){
            simu_list <- do.call(rbind, mclapply(fnames, function(fname){
                dt <- read.csv(fname)
                # Calculate volume to use weight
                dt$Biodiversity <- dt$Biodiversity * dt$farmable_ratio
                dt$Carbon <- dt$Carbon * dt$farmable_ratio
                dt$Connectivity <- dt$Connectivity * dt$farmable_ratio
                dt$Distance <- dt$Distance * dt$farmable_ratio
                
                # accumulate them
                data.frame(Biodiversity = sum(dt$Biodiversity),
                           Carbon = sum(dt$Carbon),
                           Connectivity = sum(dt$Connectivity),
                           Distance = sum(dt$Distance),
                           farmable_area = sum(dt$farmable_ratio))
            }, mc.cores = 10))
            
            simu_list %>% summarise(across(everything(), c(mean, sd))) %>% 
                mutate(land_usage = usage, scenario = item, change = change)
        }
    }))
}))

# Tweak a bit
nms <- names(smry_all_exps)
nms <- gsub("1", "mean", nms)
nms <- gsub("2", "sd", nms)
names(smry_all_exps) <- nms

fname <- file.path(dst_dir, sprintf("summary_scenarios_USG%s_mean_sd.csv", usage))
write.csv(smry_all_exps, fname, row.names = FALSE)

## Record everything
for (item in scenarios) {
    for(change in changes){
        message("Start: ", item, " ", change)
        fnames <- list.files(f_dir, pattern = sprintf("^%s_%s_", item, change), 
                             full.names = TRUE)
        
        if (length(fnames) != 0){
            simu_list <- do.call(rbind, mclapply(fnames, function(fname){
                dt <- read.csv(fname)
                # Calculate volume to use weight
                dt$Biodiversity <- dt$Biodiversity * dt$farmable_ratio
                dt$Carbon <- dt$Carbon * dt$farmable_ratio
                dt$Connectivity <- dt$Connectivity * dt$farmable_ratio
                dt$Distance <- dt$Distance * dt$farmable_ratio
                
                # accumulate them
                dt$Production_gain <- cumsum(dt$Production_gain)
                dt$Biodiversity <- cumsum(dt$Biodiversity)
                dt$Carbon <- cumsum(dt$Carbon)
                dt$Connectivity <- cumsum(dt$Connectivity)
                dt$Distance <- cumsum(dt$Distance)
                dt$farmable_ratio <- cumsum(dt$farmable_ratio)
                
                dt %>% select('Production_gain', 'Biodiversity', 'Carbon',
                              'Connectivity',  'Distance', 'farmable_ratio', 'id')
            }, mc.cores = 10))
            
            simu_mean <- simu_list %>% group_by(id) %>% 
                summarise(across(everything(), mean))
            
            simu_sd <- simu_list %>% group_by(id) %>% 
                summarise(across(everything(), sd))
            
            fname <- file.path(dst_dir, sprintf("%s_%s_USG%s_simu_mean.csv", item, change, usage))
            write.csv(simu_mean, fname, row.names = FALSE)
            fname <- file.path(dst_dir, sprintf("%s_%s_USG%s_simu_sd.csv", item, change, usage))
            write.csv(simu_sd, fname, row.names = FALSE)
        }
    }
}

# Spatial maps
f_dir <- file.path(data_dir, "spatial")

for (item in scenarios) {
    for(change in changes){
        message("Start: ", item, " ", change)
        fnames <- list.files(f_dir, pattern = sprintf("^%s_%s_", item, change), 
                             full.names = TRUE)
        if (length(fnames) != 0){
            # Read layers
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
            
            fname <- file.path(dst_dir, sprintf("%s_%s_USG%s_simu.tif", item, change, usage))
            writeRaster(out, fname)
        }
    }
}
