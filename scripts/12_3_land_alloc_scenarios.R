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
## crop production: only yield, weighted by current cropland area
## biodiversity
## carbon
## connectivity
## hybrid
## --------------------------------------------

# Load libraries
library(terra)
library(parallel)
library(dplyr)
library(stringr)
library(optparse)

# Function to allocate the land for agriculture
land_allocate_quick <- function(
        atn_yields, # the attainable yield to calculate production gain
        farmable_area, # all new farmable area
        weights, # weight for other factors except yields
        costs, # the cost by expansion
        mod = "hybrid",
        cbetas, # the prefer factor for each weight, the first one is for yield
        shift = "RG", # HY for 150% yield, CS for increase cassava plant area, or HYSC for both, RG for regular.
        production_need = 10943841, # triple the current production
        seed = 123,
        dst_dir){
    
    #  Re-calculate weights
    calib_weights <- lapply(names(cbetas), function(name){
        weights[[name]] * cbetas[name]})
    names(calib_weights) <- names(weights)
    
    # Get the production gain from expansion
    prod_gain_exp <- atn_yields * farmable_area
    
    # Get decision table for each crop
    decision_mats <- lapply(names(calib_weights$Yield), function(crp){
        # Merge all weights
        weights_crp <- do.call(
            c, lapply(1:length(calib_weights), 
                      function(n) calib_weights[[n]][[crp]])) %>% sum()
        
        values(c(weights_crp, prod_gain_exp[[crp]], costs, farmable_area)) %>% 
            data.frame() %>% mutate(cell_id = 1:nrow(.)) %>% 
            na.omit() %>% 
            rename(weight = sum, Production_gain = all_of(crp)) %>% 
            arrange(-weight)
    })
    names(decision_mats) <- names(calib_weights$Yield)
    
    # Start the allocate
    ## 1. Create a vector to store the selected cell ids. This is important because
    ## cells will not be selected up to down
    accum_cost <- data.frame(
        matrix(ncol = 9, nrow = 0, 
               dimnames=list(NULL, c("cell_id", "Production_gain", 
                                     "Biodiversity", "Carbon", 
                                     "Connectivity", "Distance", 
                                     "farmable_ratio", "crop", "id"))))
    ## 2. Select pixel by pixel
    set.seed(seed)
    prod_to_catch <- production_need
    pixel_num <- nrow(decision_mats[[1]])
    i <- 1
    
    # Do the loop
    while (prod_to_catch > 0) {
        ## Randomly select a crop to grow
        if (str_detect(shift, "CS")){
            # (5.4, 1.8, 0.9, 2, 0.9)
            crop_id <- sample(
                c(rep("maize", 54), rep("paddy", 18), 
                  rep("sorghum", 9), rep("cassava", 20), rep("beans", 9)), 1)
        } else {
            crop_id <- sample(c(rep("maize", 6), rep("paddy", 2), 
                                "sorghum", "cassava", "beans"), 1)
        }
        decision_mat <- decision_mats[[crop_id]]
        
        if (i %% 5000 == 0){
            for (crp in c("maize", "paddy", "sorghum", "cassava", "beans")){
                decision_mats[[crp]] <- decision_mats[[crp]] %>% 
                    filter(!cell_id %in% accum_cost$cell_id)
            }
            message(sprintf("Finish %s%%", 
                            round((production_need - prod_to_catch) / production_need * 100, 2)))
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
        prod_to_catch <- prod_to_catch - decision_mat$Production_gain
    }
    
    fname <- file.path(dst_dir, sprintf("%s_%s_simulation_%s.csv", mod, shift, seed))
    write.csv(accum_cost, fname, row.names = FALSE)
}

# Set directories
tdf_dir <- "/scratch/lsong36/comTCA/data/tradeoff"
dst_dir <- "/scratch/lsong36/comTCA/results/scenarios"
if (!dir.exists(dst_dir)) dir.create(dst_dir)

option_list <- list(
    make_option(c("-s", "--seed"), 
                action = "store", default = 1, type = 'integer',
                help = "The seed for the iteration [default %default]"),
    make_option(c("-m", "--mode"), 
                action = "store", type = 'character',
                help = "The mode for the simulation."),
    make_option(c("-f", "--shift"), 
                action = "store", type = 'character',
                help = "The shift for the simulation."))
opt <- parse_args(OptionParser(option_list = option_list))
seed <- opt$seed
mod <- opt$mode
shift <- opt$shift

# Mask out the protected areas
pas <- rast(file.path(tdf_dir, "protected_areas.tif"))
pas[pas == 1] <- NA

# Get specific cell sizes, the same for all layers
cellsizes <- cellSize(pas) / 1e6 * 100

farmable_area <- file.path(tdf_dir, "farmable_perc.tif") %>% 
    rast() %>% mask(pas) * cellsizes
farmable_area[farmable_area <= 0] <- NA
cropland <- file.path(tdf_dir, "cropland_perc.tif") %>% 
    rast() %>% mask(pas) * cellsizes
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
if (shift == "HY"){
    atn_yields <- rast(file.path(tdf_dir, "agro_attainable_yield_150.tif")) %>% 
        mask(farmable_area)
    prod_intes <- read.csv(file.path(tdf_dir, "production_gain_from_intensifify_150.csv"))
} else if (shift == "CS"){
    atn_yields <- rast(file.path(tdf_dir, "agro_attainable_yield.tif")) %>% 
        mask(farmable_area)
    prod_intes <- read.csv(file.path(tdf_dir, "production_gain_from_intensifify_more_cassava.csv"))
} else if (shift == "HYCS"){
    atn_yields <- rast(file.path(tdf_dir, "agro_attainable_yield_150.tif")) %>% 
        mask(farmable_area)
    prod_intes <- read.csv(file.path(tdf_dir, "production_gain_from_intensifify_150_more_cassave.csv"))
} else {
    atn_yields <- rast(file.path(tdf_dir, "agro_attainable_yield.tif")) %>% 
        mask(farmable_area)
    prod_intes <- read.csv(file.path(tdf_dir, "production_gain_from_intensifify.csv"))
}

production_need <- mean(prod_intes$Current_production) * 2 - mean(prod_intes$Attainable_production)

# Set the factors
if (mod == "Y"){
    cbetas <- c("Yield" = 1, "Biodiversity" = 0, "Carbon" = 0, 
                "Connectivity" = 0, "Distance" = 0)
} else if (mod == "B"){
    cbetas <- c("Yield" = 0, "Biodiversity" = 1, "Carbon" = 0, 
                "Connectivity" = 0, "Distance" = 0)
} else if (mod == "Ca"){
    cbetas <- c("Yield" = 0, "Biodiversity" = 0, "Carbon" = 1, 
                "Connectivity" = 0, "Distance" = 0)
} else if (mod == "Co"){
    cbetas <- c("Yield" = 0, "Biodiversity" = 0, "Carbon" = 0, 
                "Connectivity" = 1, "Distance" = 0)
} else if (mod == 'D'){
    cbetas <- c("Yield" = 0, "Biodiversity" = 0, "Carbon" = 0, 
                "Connectivity" = 0, "Distance" = 1)
} else{
    cbetas <- c("Yield" = 1/5, "Biodiversity" = 1/5, "Carbon" = 1/5, 
                "Connectivity" = 1/5, "Distance" = 1/5)
}

# Start the run
land_allocate_quick(atn_yields, farmable_area, weights, 
                    costs, mod, cbetas, shift, production_need, seed, dst_dir)
