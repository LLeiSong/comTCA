## --------------------------------------------
## Script name: circuit_sim.R
## Purpose of script: script to run circuitscape
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Run in terminal: 
## Rscript scripts/circuit_sim.R -n 50
## --------------------------------------------

# Load libraries
library(here)
library(terra)
library(ggplot2)
library(dplyr)
library(stars)
library(ini)
library(purrr)
library(JuliaCall)
library(parallel)
library(optparse)

# Set directories
ele_data_dir <- here("data/elephant")

# Read PAs
pas <- read_sf(file.path(ele_data_dir, "census.geojson")) %>% select()
suit_fname <- file.path(ele_data_dir, "suitability.asc")

# Parse inline arguments
option_list <- list(
    make_option(c("-n", "--n_iteration"), 
                action = "store", default = 10, type = 'integer',
                help = "No. of time to iteration [default %default]"))
opt <- parse_args(OptionParser(option_list = option_list))
n_iteration <- opt$n_iteration

# set up Julia
# julia_setup(installJulia = TRUE)
julia_install_package_if_needed("Circuitscape")
julia_library("Circuitscape")
config <- read.ini(
    file.path(ele_data_dir, "circuitscape_setting_template.ini"))

# Make folder to save nodes and simulations
src_dir <- file.path(ele_data_dir, "nodes")
dir.create(src_dir)
dst_dir <- file.path(ele_data_dir, "simulations")
dir.create(dst_dir)

# Iterations
suit_map <- rast(suit_fname)
cum_rasters <- lapply(1:n_iteration, function(n) {
    set.seed(100 + n)
    # Get the PAs
    nodes <- pas %>% sample_frac(0.5) %>% 
        mutate(id = 1:nrow(.)) %>% 
        rasterize(suit_map, field = "id")
    nodes[nodes == 0] <- -9999
    
    # Name it
    nm <- sprintf('iter_%s', n)
    message(sprintf("Run experiment for iteration No.%s.", n))
    
    # Save out nodes
    node_name <- file.path(src_dir, sprintf("nodes_%s.asc", nm))
    writeRaster(nodes, node_name, wopt = list(NAflag = -9999))
    
    # Revise the config parameters
    config$`Output options`$write_cur_maps <- 0
    config$`Output options`$write_volt_maps <- 0
    config$`Output options`$write_cum_cur_map_only <- "True"
    config$`Output options`$write_max_cur_maps <- "False"
    config$`Output options`$output_file <-
        file.path(dst_dir, nm)
    config$`Logging Options`$log_file <-
        file.path(dst_dir, sprintf("%s.log", nm))
    config$`Options for pairwise and one-to-all and all-to-one modes`$point_file <-
        node_name
    config$`Habitat raster or graph`$habitat_file <- suit_fname
    fname <- file.path(src_dir, sprintf("%s.ini", nm))
    write.ini(config, fname)
    
    # Run the experiment
    if (file.exists(fname)) {
        julia_command(sprintf('compute("%s")', fname))
    } else {
        warning("Failed to write out config file.")
    }
    
    # Read the files
    cum_raster <- rast(
        file.path(dst_dir, sprintf("%s_cum_curmap.asc", nm)))
    mask(cum_raster, nodes, inverse = TRUE, updatevalue = NA)
})

# Gather results
cum_raster <- do.call(c, cum_rasters) %>% 
    sum(na.rm = TRUE)
mean_raster <- do.call(c, cum_rasters) %>% 
    mean(na.rm = TRUE)
cum_fname <- file.path(ele_data_dir, "cum_curmap.tif")
mean_fname <- file.path(ele_data_dir, "mean_curmap.tif")
writeRaster(cum_raster, cum_fname)
writeRaster(mean_raster, mean_fname)
