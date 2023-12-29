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
data_dir <- "results/scenarios"
usage <- 64
land_usage <- 0.64

# USG80
data_dir <- "results/USG80"
usage <- 80
land_usage <- 0.8

dst_dir <- file.path(data_dir, "summary")
if(!dir.exists(dst_dir)) dir.create(dst_dir)
f_dir <- file.path(data_dir, "numbers")
fnames <- list.files(f_dir, full.names = TRUE)
scenarios <- unique(sapply(strsplit(basename(fnames), "_"), function(x) x[1]))
changes <- unique(sapply(strsplit(basename(fnames), "_"), function(x) x[2]))

## Overall
## The conversion ratio
cnt_ratio <- c(0.88, 0.88, 0.87, 0.33, 0.86)
names(cnt_ratio) <- c("maize", "paddy", "sorghum", "cassava", "beans")

smry_all_exps <- do.call(rbind, lapply(scenarios, function(item){
    do.call(rbind, lapply(changes, function(change){
        message("Start: ", item, " ", change)
        fnames <- list.files(f_dir, pattern = sprintf("^%s_%s_", item, change), 
                             full.names = TRUE)
        
        if (length(fnames) != 0){
            simu_list <- do.call(rbind, mclapply(fnames, function(fname){
                dt <- read.csv(fname)
                dt <- dt %>% mutate(crop_ratio = cnt_ratio[crop])
                
                # Get iteration index
                iter_n <- strsplit(gsub(".csv", "", basename(fname)), "_")[[1]][5]
                
                # Calculate volume to use weight
                dt$Biodiversity <- dt$Biodiversity * dt$farmable_ratio
                dt$Carbon <- dt$Carbon * dt$farmable_ratio
                dt$Connectivity <- dt$Connectivity * dt$farmable_ratio
                dt$Distance <- dt$Distance * dt$farmable_ratio
                dt$farmable_area_cultv <- dt$farmable_ratio * land_usage * dt$crop_ratio
                
                # accumulate them
                data.frame(Biodiversity = sum(dt$Biodiversity),
                           Carbon = sum(dt$Carbon),
                           Connectivity = sum(dt$Connectivity),
                           Distance = sum(dt$Distance),
                           farmable_area = sum(dt$farmable_ratio),
                           farmable_area_cultv = sum(dt$farmable_area_cultv),
                           num_unit = nrow(dt),
                           weight = mean(dt$weight),
                           iter_n = iter_n)
            }, mc.cores = detectCores() - 1))
            
            simu_list %>% filter(weight == max(weight)) %>% 
                select(-weight) %>% 
                mutate(land_usage = usage, scenario = item, change = change)
        }
    }))
}))

fname <- file.path(dst_dir, sprintf("summary_scenarios_USG%s_best.csv", usage))
write.csv(smry_all_exps, fname, row.names = FALSE)

## Record everything
for (item in scenarios) {
    for(chg in changes){
        message("Start: ", item, " ", chg)
        
        iter_n <- smry_all_exps %>% 
            filter(scenario == item & change == chg) %>% pull(iter_n)
        fname <- file.path(
            f_dir, sprintf("%s_%s_USG%s_simu_%s.csv", 
                           item, chg, usage, iter_n))
        
        if (length(fname) != 0){
            dt <- read.csv(fname)
            
            # Calculate volume to use weight
            dt$Biodiversity <- dt$Biodiversity * dt$farmable_ratio
            dt$Carbon <- dt$Carbon * dt$farmable_ratio
            dt$Connectivity <- dt$Connectivity * dt$farmable_ratio
            dt$Distance <- dt$Distance * dt$farmable_ratio
            dt$Production_gain <- dt$Production_gain * 
                land_usage * cnt_ratio[dt$crop]
            
            # accumulate them
            dt$Production_gain <- cumsum(dt$Production_gain)
            dt$Biodiversity <- cumsum(dt$Biodiversity)
            dt$Carbon <- cumsum(dt$Carbon)
            dt$Connectivity <- cumsum(dt$Connectivity)
            dt$Distance <- cumsum(dt$Distance)
            dt$farmable_ratio <- cumsum(dt$farmable_ratio)
            
            dt %>% select('Production_gain', 'Biodiversity', 'Carbon',
                          'Connectivity',  'Distance', 'farmable_ratio', 
                          'id')
            
            fname <- file.path(dst_dir, sprintf("%s_%s_USG%s_simu_best.csv", item, chg, usage))
            write.csv(dt, fname, row.names = FALSE)
        }
    }
}

# Spatial maps
f_dir <- file.path(data_dir, "spatial")

for (item in scenarios) {
    for(chg in changes){
        message("Start: ", item, " ", chg)
        
        iter_n <- smry_all_exps %>% 
            filter(scenario == item & change == chg) %>% pull(iter_n)
        fname <- file.path(
            f_dir, sprintf("%s_%s_USG%s_simu_%s.tif", 
                           item, chg, usage, iter_n))
        
        if (length(fname) != 0){
            # Read layers
            lyrs <- rast(fname)
            
            # For crops
            crps <- do.call(c, lapply(1:5, function(crp){
                lyr <- lyrs
                lyr[lyr != crp] <- NA
                lyr[lyr == crp] <- 1
                lyr
            }))
            
            # Overall
            ovr <- lyrs
            ovr[ovr > 0] <- 1
            
            out <- c(ovr, crps)
            names(out) <- c('farmland', "maize", "paddy", "sorghum", "cassava", "beans")
            
            fname <- file.path(dst_dir, sprintf("%s_%s_USG%s_best.tif", item, chg, usage))
            writeRaster(out, fname)
        }
    }
}
