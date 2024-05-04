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
library(stringr)

## Directories
data_dir <- "results/scenarios" # scenarios/profit
dst_dir <- file.path(data_dir, "summary")
if(!dir.exists(dst_dir)) dir.create(dst_dir)

## Extract parameters
f_dir <- file.path(data_dir, "numbers")
fnames <- list.files(f_dir, full.names = TRUE)
land_usages <- unique(sapply(strsplit(basename(fnames), "_"), function(x) x[3]))
land_usages <- as.numeric(str_extract(land_usages, "[0-9]{2}")) / 100
scenarios <- unique(sapply(strsplit(basename(fnames), "_"), function(x) x[1]))
changes <- unique(sapply(strsplit(basename(fnames), "_"), function(x) x[2]))

## Overall
smry_all_exps <- do.call(rbind, lapply(land_usages, function(land_usage){
    do.call(rbind, lapply(scenarios, function(item){
        do.call(rbind, lapply(changes, function(change){
            message("Start: ", item, " ", change, " ", land_usage)
            fname <- list.files(
                f_dir, pattern = sprintf(
                    "^%s_%s_USG%2d", item, change, as.integer(land_usage * 100)), 
                full.names = TRUE)
            
            if (length(fname) != 0){
                dt <- read.csv(fname)
                
                # Calculate volume to use weight
                dt$Biodiversity <- dt$Biodiversity * dt$farmable_ratio
                dt$Carbon <- dt$Carbon * dt$farmable_ratio
                dt$Connectivity <- dt$Connectivity * dt$farmable_ratio
                dt$Distance <- dt$Distance * dt$farmable_ratio
                
                # Accumulate them
                data.frame(Biodiversity = sum(dt$Biodiversity),
                           Carbon = sum(dt$Carbon),
                           Connectivity = sum(dt$Connectivity),
                           Distance = sum(dt$Distance),
                           farmable_area = sum(dt$farmable_ratio),
                           num_unit = nrow(dt),
                           weight = mean(dt$weight)) %>% 
                    mutate(land_usage = land_usage, 
                           scenario = item, change = change)
            }
        }))
    }))
}))

# or summary_profit.csv for profit
fname <- file.path(dst_dir, "summary_scenarios_USGs.csv")
write.csv(smry_all_exps, fname, row.names = FALSE)
