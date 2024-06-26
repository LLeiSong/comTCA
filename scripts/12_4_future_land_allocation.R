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

## Evaluate the effects of weight on reducing cost
## Y(*), 0% - 100% (20%)

## Usage
## Rscript 12_4_future_land_allocation.R -s 1 -c "0.2,0.2,0.2,0.2,0.2" -o Y100 -l 0.64 -d scenarios
## --------------------------------------------

# Load libraries
library(terra)
library(parallel)
library(dplyr)
library(stringr)
library(optparse)
library(prioritizr)

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
    land_usage = 0.64, # the percentage of annual land use. 0.64, 0.8, or 1.0
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
    cnt_ratio <- c(0.88, 0.88, 0.87, 0.33, 0.86)
    names(cnt_ratio) <- c("maize", "paddy", "sorghum", "cassava", "beans")
    prod_gain_exp <- do.call(c, lapply(names(atn_yields), function(crp){
        atn_yields[[crp]] * (farmable_area * land_usage) * cnt_ratio[[crp]] / 1e6}))
    
    # Get the ratio of each crop
    cassava <- ifelse(
        str_detect(scenario, "CASS"), str_extract(scenario, "[0-9]+"), 1)
    if (cassava == 1){
        ratios <- c(58, 20, 6, 6, 10)
    } else if (cassava == 2){
        ratios <- c(54.3, 18.7, 5.6, 12, 9.4)
    } else if (cassava == 3){
        ratios <- c(50.6, 17.5, 5.2, 18, 8.7)
    } else if (cassava == 4){
        ratios <- c(46.9, 16.2, 4.8, 24, 8.1)
    } else if (cassava == 5){
        ratios <- c(43.2, 14.9, 4.5, 30, 7.4)
    }
    
    # The workflow is like:
    ## First, do an initial run to get a not-bad but not the best result
    ## Starting from there, to converge to the best result
    
    ########################## Initial run ##########################
    set.seed(seed)
    # Merge all weights
    pus <- do.call(c, lapply(names(calib_weights$Yield), function(crp){
        pu <- do.call(c, lapply(1:length(calib_weights), 
                      function(n) calib_weights[[n]][[crp]])) %>% sum()
        pu * farmable_area
    }))
    names(pus) <- names(calib_weights$Yield)
    
    # Calculate the initial targets
    # sum(area * ratios * yields) = production_need
    areas <- production_need / 1e6 / sum(ratios * c(2.77, 11.75, 7.68, 24.07, 2.62))
    prod_targets <- areas * ratios * c(2.77, 11.75, 7.68, 24.07, 2.62)
    names(prod_targets) <- c("maize", "paddy", "sorghum", "cassava", "beans")
    
    targets <- do.call(cbind, lapply(names(prod_targets), function(crp){
        prod_target_crp <- prod_targets
        prod_target_crp[setdiff(names(prod_targets), crp)] <- 0
        prod_target_crp
    }))
    colnames(targets) <- rownames(targets)
    
    p <- problem(
        pus,
        zones(prod_gain_exp, prod_gain_exp,
              prod_gain_exp, prod_gain_exp, prod_gain_exp,
            zone_names = names(pus),
            feature_names = names(prod_gain_exp))) %>%
        add_min_set_objective() %>%
        add_absolute_targets(targets) %>%
        add_binary_decisions() %>% 
        add_gurobi_solver(threads = 8, gap = 0, verbose = FALSE)
    s <- solve(p)
    
    ########################## Start the converge ##########################
    # Good, the problem was solved, now let's converge the solver to match the
    # area requirements
    while (TRUE){
        vals <- (farmable_area * land_usage * s) %>% values() %>% colSums(na.rm = T)
        ratios_run <- vals / sum(vals) * 100
        
        # Stop if converge
        if (all(abs((ratios - ratios_run) / ratios_run) < 0.01)){
            break
        }
        
        prod_targets <- prod_targets * ((ratios - ratios_run) / ratios_run + 1)
        prod_targets <- production_need / 1e6 * (prod_targets / sum(prod_targets))
        
        targets <- do.call(cbind, lapply(names(prod_targets), function(crp){
            prod_target_crp <- prod_targets
            prod_target_crp[setdiff(names(prod_targets), crp)] <- 0
            prod_target_crp
        }))
        colnames(targets) <- rownames(targets)
        
        p <- problem(
            pus,
            zones(prod_gain_exp, prod_gain_exp,
                  prod_gain_exp, prod_gain_exp, prod_gain_exp,
                  zone_names = names(pus),
                  feature_names = names(prod_gain_exp))) %>%
            add_min_set_objective() %>%
            add_absolute_targets(targets) %>%
            add_binary_decisions() %>% 
            add_gurobi_solver(threads = 8, gap = 0, verbose = FALSE)
        s <- solve(p)
    }
    
    ############################ Allocate values ############################
    # Great, the solver has been converged, now allocates all values to pixels
    if (spatial){
        ## Spatial map
        crps <- c("maize", "paddy", "sorghum", 'cassava', "beans")
        agro_land <- category_layer(s)
        names(agro_land) <- "crop"; agro_land[agro_land == 0] <- NA
        levels(agro_land) <- data.frame(id = 1:5, crop = crps)
        dst_fname <- file.path(
            spatial_dir, sprintf("%s_%s_USG%s.tif", name, scenario, round(land_usage * 100, 0)))
        writeRaster(agro_land, dst_fname, datatype = "INT1U", overwrite = TRUE)
    }
    
    ## Values
    lyrs <- c(agro_land, costs, farmable_area, 
              max(pus * s), max(prod_gain_exp * s) * 1e6)
    names(lyrs)[7:8] <- c("weight", "production_gain")
    vals <- values(lyrs) %>% na.omit() %>% data.frame()
    
    dst_fname <- file.path(
        num_dir, sprintf("%s_%s_USG%s.csv", name, scenario, round(land_usage * 100, 0)))
    write.csv(vals, dst_fname, row.names = FALSE)
}

# Set directories
tdf_dir <- "data/tradeoff"
intes_dir <- "results/intensification"
dst_dir <- "results"

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
                help = "The percentage of land usage, [0.64, 0.8, 1.0]."),
    make_option(c("-p", "--production_need"), 
                action = "store", default = 1, type = 'numeric',
                help = "The production to catch from land expansion."),
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
## Ecological costs
costs <- file.path(
    tdf_dir, c("bio_index.tif", "carbon_density.tif", 
               "conn_index.tif", "agro_travel_time.tif")) %>% 
    rast() %>% mask(farmable_area)
names(costs) <- c("Biodiversity", "Carbon", "Connectivity", "Distance")

## Get efficiency-based weights for ecological costs
weights <- lapply(costs, function(lyr){
    lyrs <- lyr / atn_yields
    names(lyrs) <- names(atn_yields)
    lyrs
})

## Get yield weights
## factors to prefer to the smallest value.
yields_as_weights <- 1 / atn_yields

## Put them together
weights <- c(list(yields_as_weights), weights)
names(weights) <- c("Yield", "Biodiversity", "Carbon", "Connectivity", "Distance")

## Here to normalize the layers across all layers for each element
weights <- lapply(weights, function(wt_lyrs){
    glb_min <- min(global(wt_lyrs, "min", na.rm = TRUE)[, 1])
    glb_max <- max(global(wt_lyrs, "max", na.rm = TRUE)[, 1])
    
    (wt_lyrs - glb_min) / c(glb_max - glb_min)
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

# To test the effects of weights
# # Set static parameters
# seed <- 123
# scenario <- 'Y100'
# land_usage <- 0.64
# yield_num <- ifelse(
#     str_detect(scenario, "Y"), str_extract(scenario, "[0-9]+"), 100)
# cassava <- ifelse(
#     str_detect(scenario, "CASS"), str_extract(scenario, "[0-9]+"), 1)
# 
# # Make directory for the results
# dst_dir <- file.path(dst_dir, "profit")
# if (!dir.exists(dst_dir)) dir.create(dst_dir)
# 
# # Loop over different weight values
# wts <- seq(0.1, 0.9, 0.1)
# factors <- c("Yield", "Biodiversity", "Carbon", "Connectivity", "Distance")
# 
# for (ft in factors[-1]){
#     for (wt in wts){
#         # Set weights
#         cbetas <- rep(0, 5)
#         cbetas[1] <- 1 - wt
#         cbetas[which(factors == ft)] <- wt
#         
#         # Set the factors
#         cbetas <- c("Yield" = cbetas[1], "Biodiversity" = cbetas[2], "Carbon" = cbetas[3], 
#                     "Connectivity" = cbetas[4], "Distance" = cbetas[5])
#         idx <- which(cbetas != 0)[-1]
#         exp_name <- paste(
#             c("Y", sprintf("%s%s", c("Y", "B", "Ca", "Co", "D")[idx], 
#                            cbetas[idx] * 10)), collapse = "")
#         message(sprintf("%s, %s, %s", ft, wt, exp_name))
#         print(cbetas)
#         
#         # Run experiments
#         land_allocate(
#             atn_yields, farmable_area, weights, costs,
#             exp_name, cbetas, scenario, land_usage, production_need,
#             TRUE, seed, dst_dir)
#     }
# }