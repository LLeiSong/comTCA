# Title     : Calculate multi-species landscape connectivity
# Objective : The main function to calculate the multi-species level habitat 
#             connectivity. Use the mean suitability values as the conductance. 
# Created by: Lei Song
# Created on: 09/09/23
# Note      : Use cg+amg solver to avoid RAM crisis. The user can modify the
#             loop part to an independent script to run on HPC.

# Load libraries
library(sf)
library(here)
library(terra)
library(dplyr)
library(stringr)
library(ini)
library(purrr)
library(JuliaCall)
library(parallel)
library(optparse)
sf_use_s2(FALSE)

# Set directories
conn_data_dir <- here("data/connectivity")
circ_dir <- file.path(conn_data_dir, "circuit")
if (!dir.exists(circ_dir)) dir.create(circ_dir)

# Parse inline arguments
option_list <- list(
    make_option(c("-i", "--initial"), 
                action = "store_true", default = TRUE,
                help = "Initialize inputs [default]"),
    make_option(c("-n", "--n_iteration"), 
                action = "store", default = 10, type = 'integer',
                help = "No. of time to iteration [default %default]"))
opt <- parse_args(OptionParser(option_list = option_list))
initial <- opt$initial
n_iteration <- opt$n_iteration

if (initial){
    # Get the boundary
    bry <- read_sf(here("data/geoms/mainland_tanzania.geojson")) %>% select()
    
    # Prepare inputs for circuitscape
    ## Write out asc file for suitability
    fname <- file.path(conn_data_dir, "multi_species_suitability.tif")
    suitability <- rast(fname)
    suitability <- mask(suitability, bry)
    suitability[is.na(suitability)] <- -9999
    writeRaster(suitability, file.path(circ_dir, 'suitability.asc'),
                wopt = list(NAflag = -9999))
    
    ## Filter the protected areas
    ## Remove culture and forest based PAs by assuming dense forest reserve highly 
    ## impossibly serve as a long-term home range as they are mainly savanna animals.
    ## These areas can still be used by animals if they are important, 
    ## just not be used for nodes in circuitscape.
    pas <- read_sf(
        file.path("data/protected_area", "WDPA_WDOECM_Jan2023_Public_TZA",
                  "WDPA_WDOECM_Jan2023_Public_TZA.gdb"),
        layer = "WDPA_WDOECM_poly_Jan2023_TZA") %>% 
        select(NAME, DESIG_ENG) %>% 
        rename(Geometry = SHAPE)
    pas <- pas %>% 
        filter(!DESIG_ENG %in% c("Marine Reserve", "Marine Park",
                                 "Sanctuary and Closed Forest Reserve")) %>% 
        filter(!str_detect(DESIG_ENG, "forest|Forest"))
    
    ## Double check "World Heritage Site (natural or mixed)" only contains parks
    ## after this step
    
    ## Remove some small ones, smaller than Lake Manyara park
    pas <- pas %>% mutate(area = st_area(.)) %>% 
        filter(area > units::set_units(1e08, "m^2")) %>% 
        select(NAME)
    
    # Could remove some duplicates
    
    st_write(pas, file.path(circ_dir, "pas_selected.geojson"))
}

# Run simulations
## Reload the data
pas <- read_sf(file.path(circ_dir, "pas_selected.geojson")) %>% select()
suit_fname <- file.path(circ_dir, 'suitability.asc')

# set up Julia
# julia_setup(installJulia = TRUE)
julia_install_package_if_needed("Circuitscape")
julia_library("Circuitscape")
config <- read.ini(
    file.path(circ_dir, "circuitscape_setting_template.ini"))

# Make folder to save nodes and simulations
src_dir <- file.path(circ_dir, "nodes")
dst_dir <- file.path(circ_dir, "simulations")
for (dir_to in c(src_dir, dst_dir)){
    if (!dir.exists(dir_to)) dir.create(dir_to)}

# Iterations
suit_map <- rast(suit_fname)
cum_rasters <- lapply(1:n_iteration, function(n) {
    set.seed(100 + n)
    
    # Get the PAs
    ## Use half PAs and randomly generate mini habitat patches as nodes
    nodes <- pas %>% sample_frac(0.5) %>%
        st_sample(size = rep(1, nrow(.))) %>% 
        st_as_sf() %>% mutate(id = 1:nrow(.)) %>% vect() %>% 
        buffer(width = 5000) %>% 
        rasterize(suit_map, field = "id")
    nodes[nodes == 0] <- -9999
    nodes <- mask(nodes, suit_map)
    
    # Name it
    nm <- sprintf('iter_%s', n)
    message(sprintf("Run experiment for iteration No.%s.", n))
    
    # Save out nodes
    node_name <- file.path(src_dir, sprintf("nodes_%s.asc", nm))
    writeRaster(nodes, node_name, wopt = list(NAflag = -9999), overwrite = TRUE)
    
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
    config$`Calculation options`$parallelize <- 'False'
    config$`Calculation options`$solver <- "cg+amg"
    config$`Calculation options`$print_timings <- 'True'
    config$`Calculation options`$print_rusages <- 'True'
    config$`Calculation options`$preemptive_memory_release <- 'True'
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
cum_fname <- file.path(circ_dir, "cum_curmap.tif")
mean_fname <- file.path(circ_dir, "mean_curmap.tif")
writeRaster(cum_raster, cum_fname)
writeRaster(mean_raster, mean_fname)
