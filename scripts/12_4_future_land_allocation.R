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
## yield + biodiversity, cbetas = "0.5,0.5,0,0,0"
## yield + carbon, cbetas = "0.5,0,0.5,0,0"
## yield + connectivity, cbetas = "0.5,0,0,0.5,0"
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
    land_usage = 0.653, # the percentage of annual land use. 0.653, 0.8, or 1.0
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
    
    fname <- file.path(
        num_dir, sprintf("%s_%s_USG%s_simu_%s.csv", 
                         name, scenario, round(land_usage * 100, 0), seed))
    write.csv(accum_cost, fname, row.names = FALSE)
    
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
        writeRaster(agro_land, dst_fname, datatype = "INT1U")  
    }
}

# Set directories
tdf_dir <- "data/tradeoff"
intes_dir <- "results/intensification"
dst_dir <- "results/scenarios"
if (!dir.exists(dst_dir)) dir.create(dst_dir)

# Mask out the protected areas
pas <- rast(file.path(tdf_dir, "protected_areas.tif"))
pas[pas == 1] <- NA

# Get specific cell sizes in ha, the same for all layers
cellsizes <- cellSize(pas) / 1e6 * 100

# Convert percentage to area in ha
farmable_area <- file.path(tdf_dir, "farmable_perc.tif") %>% 
    rast() %>% mask(pas) * cellsizes
farmable_area[farmable_area <= 0] <- NA
cropland <- file.path(tdf_dir, "cropland_perc.tif") %>% 
    rast() %>% mask(pas) * cellsizes
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
                help = "The percentage of land usage, [0.653, 0.8, 1.0]."))
opt <- parse_args(OptionParser(option_list = option_list))
seed <- opt$seed
cbetas <- opt$cbetas
cbetas <- as.numeric(strsplit(cbetas, ",")[[1]])
scenario <- opt$scenario
land_usage <- opt$land_usage
yield_num <- ifelse(
    str_detect(scenario, "Y"), str_extract(scenario, "[0-9]+"), 100)
cassava <- ifelse(
    str_detect(scenario, "CASS"), str_extract(scenario, "[0-9]+"), 1)

# Get the right attainable yield layer
atn_yield_fname <- file.path(
    tdf_dir, sprintf("agro_attainable_yield%s.tif", 
                     ifelse(yield_num == 100, "", paste0("_", yield_num))))
atn_yields <- rast(atn_yield_fname) %>% mask(farmable_area)
message(sprintf("Read attainable layer: %s.", basename(atn_yield_fname)))

# Gather all weight inputs and standardize them
## Ecological costs
costs <- file.path(
    tdf_dir, c("bio_index.tif", "carbon_density.tif", 
               "conn_index.tif", "agro_travel_time.tif")) %>% 
    rast() %>% mask(farmable_area)
names(costs) <- c("Biodiversity", "Carbon", "Connectivity", "Distance")

## Get efficiency-based weights for ecological costs
weights <- lapply(costs, function(lyr){
    lyrs <- 1 - stretch(lyr / atn_yields, minv = 0, maxv = 1)
    names(lyrs) <- names(atn_yields)
    lyrs
})

## Get yield weights
yields_as_weights <- 1 - stretch(1 / atn_yields, minv = 0, maxv = 1)

## Put them together
weights <- c(list(yields_as_weights), weights)
names(weights) <- c("Yield", "Biodiversity", "Carbon", "Connectivity", "Distance")

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
