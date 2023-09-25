# Title     : The organizing function for refining expert ranges
# Objective : Use the HPC high memory node to refine the rasterized expert range
#           : maps. This script can put in scheduler.
# Created by: Lei Song
# Created on: 07/26/23

# Load libraries
library(sf)
library(parallel)
library(optparse)

# Set directories
data_dir <- "/scratch/lsong36/comTCA/data"

# Parse inline arguments
option_list <- list(
    make_option(c("-t", "--taxon"), 
                action = "store", type = 'character',
                help = "The taxon of the species."),
    make_option(c("-r", "--range"), 
                action = "store", type = 'character',
                help = "The range, global or tanzania."))
opt <- parse_args(OptionParser(option_list = option_list))
taxon <- opt$taxon
range <- opt$range

# Read polygons
plys <- st_read(file.path(data_dir, "expert_range_maps/clean_range_maps",
                          sprintf("%s_tz_relevant.geojson", taxon)))
message(sprintf("%s species to refine.", length(unique(plys$sci_name))))

mclapply(unique(plys$sci_name), function(species){
    cmd <- sprintf('Rscript scripts/refine_expert_ranges.R -s "%s" -t %s -r %s',
                   species, taxon, range)
    system(cmd)
}, mc.cores = 15)