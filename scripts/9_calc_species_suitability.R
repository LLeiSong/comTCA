# Title     : Calculate species-level habitat suitability
# Objective : The main function to calculate the habitat suitability. 
# Created by: Lei Song
# Created on: 09/09/23
# General steps: Step 1: Use the refined range map to make samples, together 
#                with other environmental variables to create SDM.
#                Step 2: Query occurrence from GBIF.
#                Step 3: Ensemble these two models and do prediction.

# Load libraries
## These can be installed from CRAN.
library(sf)
library(here)
library(terra)
library(dplyr)
library(rgbif)
library(lubridate)
library(rgbif)
library(pbmcapply)
library(caret)
library(biomod2)
# remotes::install_github("ropensci/scrubr")
library(scrubr)
sf_use_s2(FALSE)

# Set directories
data_dir <- here("data")
conn_data_dir <- here("data/connectivity")
range_dir <- file.path(
    data_dir, "expert_range_maps/refined_range_rasters/tanzania_100/mammals")
occ_dir <- file.path(conn_data_dir, 'occurrence')
sdm_dir <- file.path(conn_data_dir, "sdm")
dst_dir <- file.path(conn_data_dir, 'suitability')
for (dir_to in c(occ_dir, sdm_dir, dst_dir)){
    if (!dir.exists(dir_to)) dir.create(dir_to)}

# Get the boundary
bry <- read_sf(here("data/geoms/tanzania_coarse_bry.geojson"))

# Define functions
## Query GBIF occurrences
query_gbif <- function(species){
    message(species)
    
    ## Set the time interval for querying on GBIF
    start_year <- 2013
    year <- sprintf('%s,%s',  start_year, year(Sys.Date()))
    
    # Search
    occ <- occ_search(
        scientificName = species,
        country = "TZ",
        hasCoordinate = TRUE,
        limit = 200000,
        year = year,
        hasGeospatialIssue = FALSE)
    
    # Clean the dataset
    occ <- dframe(occ$data) %>%
        coord_impossible() %>%
        coord_incomplete() %>%
        coord_unlikely()
    
    # Convert to sf
    occ <- st_as_sf(
        occ[, c("decimalLongitude", "decimalLatitude")], 
        coords = c("decimalLongitude", "decimalLatitude"),
        crs = 4326)
    rm(start_year, year, nm_search)
    
    fn <- file.path(occ_dir, sprintf("%s.geojson", gsub(" ", "_", species)))
    st_write(occ, fn)
}

## Fit SDM with different parameters
fit_sdm <- function(species){
    message(species)
    
    # Read variables, range map and occurrence
    vars <- rast(
        file.path("data/connectivity", "variables/variables.tif"))
    range_map <- rast(
        file.path("data/expert_range_maps/refined_range_rasters/",
                  "tanzania_100/mammals",
                  sprintf("%s.tif", gsub(" ", "_", species))))
    
    # Make a non-NA mask
    msk <- lapply(1:nlyr(vars), function(n) !is.na(vars[[n]]))
    msk <- sum(do.call(c, msk)) == nlyr(vars)
    msk[msk == 0] <- NA
    # Transform vars
    vars <- st_as_stars(vars)
    vars <- split(vars, "band")
    
    # Monte Carlo simulation
    ## Important parameters: 
    ## prob_pick_pooled_gain (prefer to lower value here) - 0, 0.1, 0.2
    ## ntry (only if use prob_pick_pooled_gain) - 1, 10, 15, 20
    ## sample_size (prefer to higher value here) - 0.8, 0.9, 1.0
    ## max_depth (prefer to slightly higher value) - 15, 20, 30
    ## ndim (recommended to use small value) - 2, 3, 4
    ## scoring_metric: depth, adj_depth
    ## Since the objective is to rank the environmental conditions, not to
    ## really isolate the outliers, so not use prob_pick_pooled_gain and ntry.
    ## max_depth is rather be higher. 
    ## Not use density based methods because the density of elephants in Tanzania
    ## is not linear related to the environmental conditions.
    
    # Define the sets of preferred parameters
    params <- expand.grid(
        sample_size = c(0.8, 0.9, 1.0),
        max_depth = seq(20, 60, 10),
        ndim = c(2, 3, 4),
        scoring_metric = c("depth", "adj_depth"))
    
    models <- lapply(1:nrow(params), function(n) {
        # Parameters
        params_run <- params[n, ]
        
        # Dynamically make pseudo samples of 1000
        set.seed(123 + n) # make sure variability for each parameter pair
        pseudo_occ <- spatSample(
            range_map, 1000, na.rm = TRUE, 
            as.point = TRUE, exhaustive = TRUE) %>% st_as_sf()
        names(pseudo_occ)[1] <- "observation"
        
        message(sprintf("No. of samples: %s.", nrow(pseudo_occ)))
        
        # Split occurrence with 3 folds
        set.seed(234 + n)
        flds <- createFolds(1:nrow(pseudo_occ), 
                            k = 3, list = TRUE, 
                            returnTrain = FALSE)
        
        # Cross validation
        it_sdms <- lapply(flds, function(ids) {
            # Split
            occ_sf <- pseudo_occ[setdiff(1:nrow(pseudo_occ), ids), ]
            occ_test_sf <- pseudo_occ[ids, ]
            
            # Do modeling
            isotree_po(
                obs = occ_sf,
                obs_ind_eval = occ_test_sf,
                variables = vars,
                sample_size = 0.8,
                ndim = params_run$ndim, 
                max_depth = params_run$max_depth,
                scoring_metric = params_run$scoring_metric,
                seed = 10L,
                nthreads = 6,
                response = FALSE,
                spatial_response = FALSE,
                check_variable = FALSE,
                visualize = FALSE)})
        
        # Collect results
        eval_mean <- do.call(rbind, lapply(it_sdms, function(run) {
            eval_test <- run$eval_test
            tibble("cvi25" = eval_test$po_evaluation$cvi$`cvi with 0.25`,
                   "cvi50" = eval_test$po_evaluation$cvi$`cvi with 0.5`,
                   "cvi75" = eval_test$po_evaluation$cvi$`cvi with 0.75`,
                   "cbi" = eval_test$po_evaluation$boyce$cor,
                   "auc_ratio" = eval_test$po_evaluation$roc_ratio$auc_ratio,
                   "sensitivity" = eval_test$pb_evaluation$sensitivity,
                   "specificity" = eval_test$pb_evaluation$specificity,
                   "TSS" = eval_test$pb_evaluation$TSS$`Optimal TSS`,
                   "auc" = eval_test$pb_evaluation$roc$auc,
                   `Jaccard's similarity index` = eval_test$pb_evaluation$`Jaccard's similarity index`,
                   "f-measure" = eval_test$pb_evaluation$`f-measure`,
                   `Overprediction rate` = eval_test$pb_evaluation$`Overprediction rate`,
                   `Underprediction rate` = eval_test$pb_evaluation$`Underprediction rate`)
        })) %>% summarise(across(everything(), mean))
        
        # I like terra more, so convert to rast
        prediction <- do.call(c, lapply(it_sdms, function(run){
            rast(run$prediction)
        })) %>% rast() %>% mean(na.rm = TRUE) %>% st_as_stars()
        
        # Report the progress
        message(sprintf("%.2f %s finished (%s).", 
                        n / nrow(params) * 100, "%", Sys.time()))
        
        # Stack the results
        list(params = params_run,
             eval_mean = eval_mean,
             prediction = prediction)
    })
    
    fn <- file.path(sdm_dir, sprintf("models_%s.rda", gsub(" ", "_", species)))
    save(models, file = fn)
}

# Calculate surface suitability
# and make evaluation using the occurrences
suitability <- function(species){
    # Load SDMs
    fn <- file.path(sdm_dir, sprintf("models_%s.rda", gsub(" ", "_", species)))
    load(fn)
    
    # Get prediction raster
    suit <- do.call(c, lapply(models, function(mod){
        mod$prediction
    })) %>% merge(name = "band") %>% st_apply(1:2, mean)
    
    # Get mean of model evaluation
    eval <- lapply(models, function(mod){
        mod$eval_mean
    }) %>% bind_rows() %>% summarise(across(everything(), mean))
    
    # Make the independent evaluation
    occ <- st_read(
        file.path('data/connectivity/occurrence', 
                  sprintf("%s.geojson", gsub(" ", "_", species))))
    occ <- rasterize(occ, rast(models[[1]]$prediction))
    occ <- as.points(occ) %>% st_as_sf() %>% 
        st_coordinates()
    
    eval_real <- lapply(1:length(models), function(n){
        # Get background points
        mod <- models[[n]]
        pred <- rast(mod$prediction)
        bg_pts <- rasterize(occ, pred, background = 0)
        bg_pts <- mask(bg_pts, pred)
        set.seed(123 + n)
        bg_pts <- spatSample(
            bg_pts, nrow(occ), "stratified", 
            exhaustive = TRUE)
        pred_pa <- extract(pred, bg_pts$cell)[[1]]
        
        boi <- Boyce(obs = occ, pred = pred)$Boyce
        auc <- AUC(obs = bg_pts$last, pred = pred_pa, method = "rank")$AUC
        msus <- threshMeasures(
            obs = bg_pts$last, pred = pred_pa, thresh = "maxTSS", 
            measures = c("Sensitivity", "Specificity", "TSS"), 
            standardize = FALSE)$ThreshMeasures
        as.data.frame(t(rbind(boi, auc, msus)))
    }) %>% bind_rows() %>% summarise(across(everything(), mean))
    
    list(suit = suit, 
         eval = eval,
         eval_real = eval_real)
}

# Get species
species_list <- file.path(conn_data_dir, "select_species.csv") %>% 
    read.csv() %>% pull(sci_name)

# Query occurrence
for (species in species_list){
    query_gbif(species)}

# Fit the final model
for (species in species_list){
    fit_sdm(species)}

# Do prediction
preds_suit <- lapply(species_list, suitability)

fn <- file.path(dst_dir, "preds_suit.rda")
save(preds_suit, file = fn)

suit <- do.call(c, lapply(preds_suit, function(each){
    rast(each$suit)
})) %>% mean(na.rm = TRUE)

eval <- lapply(preds_suit, function(each){
    each$eval
}) %>% bind_rows()

eval_real <- lapply(preds_suit, function(each){
    each$eval_real
}) %>% bind_rows()
