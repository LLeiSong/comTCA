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
fit_sdm_em <- function(species){
    message(species)
    
    # Read variables, range map and occurrence
    vars <- rast(
        file.path("data/connectivity", "variables/variables.tif"))
    
    # Make a non-NA mask
    msk <- lapply(1:nlyr(vars), function(n) !is.na(vars[[n]]))
    msk <- sum(do.call(c, msk)) == nlyr(vars)
    msk[msk == 0] <- NA
    
    # Range map
    range_map <- rast(
        file.path("data/expert_range_maps/refined_range_rasters/",
                  "tanzania_100/mammals",
                  sprintf("%s.tif", gsub(" ", "_", species))))
    range_map <- subst(range_map, NA, 0) %>% 
        resample(msk, method = "near") * msk
    
    # Read occurrences
    occ <- st_read(
        file.path('data/connectivity/occurrence', 
                  sprintf("%s.geojson", gsub(" ", "_", species))))
    occ <- rasterize(occ, vars[[1]]) * msk
    occ <- as.points(occ) %>% st_as_sf() # thinned
    
    # Dynamically make pseudo samples of 1000
    repeats <- 1:10
    em_models <- lapply(repeats, function(n){
        set.seed(123 + n) # make sure variability for each parameter pair
        pseudo_occ <- spatSample(
            range_map, 1000, method = "stratified", na.rm = TRUE, 
            as.point = TRUE, exhaustive = TRUE) %>% st_as_sf()
        names(occ) <- names(pseudo_occ) <- c("observation", "geometry")
        
        obs_data <- rbind(occ, pseudo_occ)
        obs_data <- BIOMOD_FormatingData(
            resp.var = obs_data$observation,
            expl.var = vars,
            dir.name = sdm_dir,
            resp.xy = st_coordinates(obs_data),
            resp.name = species)
        
        models <- BIOMOD_Modeling(
            bm.format = obs_data,
            models = c("FDA", "GAM", "GBM", "GLM", "MARS", "MAXNET", "RF", "XGBOOST"),
            modeling.id = "allModels",
            CV.strategy = "kfold",
            CV.nb.rep = 2,
            CV.perc = 0.7,
            CV.k = 3,
            seed.val = 123,
            nb.cpu = 12)
        
        model_em <- BIOMOD_EnsembleModeling(
            bm.mod = models, # model list from previous step
            models.chosen = "all", # use all of the models
            em.by = "all",
            em.algo = c('EMwmean'),
            metric.select = c('TSS'),
            metric.select.thresh = c(0.7),
            metric.eval = c("KAPPA", "TSS", "ROC"),
            nb.cpu = 12)
        
        model_proj <- BIOMOD_EnsembleForecasting(
            bm.em = model_em,
            proj.name = "Current",
            new.env = vars,
            models.chosen = "all",
            metric.binary = "all",
            metric.filter = "all",
            nb.cpu = 12)
        
        pred <- get_predictions(model_proj)
        eval <- get_evaluations(model_em)
        list(pred = pred, eval = eval)
    })
    
    preds <- do.call(c, lapply(em_models, function(mod){
        mod$pred
    }))
    
    fn <- file.path(sdm_dir, sprintf("suit_%s.tif", gsub(" ", "_", species)))
    writeRaster(preds, fn)
    
    fn <- file.path(sdm_dir, sprintf("eval_%s.rda", gsub(" ", "_", species)))
    save(em_models, file = fn)
}

# Get species
species_list <- file.path(conn_data_dir, "select_species.csv") %>% 
    read.csv() %>% pull(sci_name)

# Query occurrence
for (species in species_list){
    query_gbif(species)}

# Fit the final model
for (species in species_list){
    fit_sdm_em(species)}
