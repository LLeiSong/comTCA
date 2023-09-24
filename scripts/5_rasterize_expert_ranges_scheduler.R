# Title     : The scheduler to rasterize expert range maps
# Objective : Use HPC node as many as possible to speed up the rasterization.
# Created by: Lei Song
# Created on: 07/26/23

# Load libraries
library(sf)
library(optparse)

# Set directories
data_dir <- "/scratch/lsong36/comTCA/data"

# Parse inline arguments
option_list <- list(
    make_option(c("-t", "--taxon"), 
                action = "store", type = 'character',
                help = "The taxon of the species."))
opt <- parse_args(OptionParser(option_list = option_list))
taxon <- opt$taxon

# Read polygons
plys <- st_read(file.path(data_dir, "expert_range_maps/clean_range_maps",
                          sprintf("%s_tz_relevant.geojson", taxon)))
snames <- unique(plys$sci_name)
message(sprintf("%s species to rasterize for %s.", length(snames), taxon))

for (sname in snames){
    system(sprintf("sbatch schedulers/run_rasterize_expert_ranges.sh '%s' %s", 
                   sname, taxon))
}
