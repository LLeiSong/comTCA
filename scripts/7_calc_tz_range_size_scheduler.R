# Title     : Calculate species-level richness index
# Objective : The main function to calculate the species-level richness index
# Created by: Lei Song
# Created on: 07/26/23
# Use HPC nodes to do this

# Load libraries
library(terra)
library(dplyr)
library(parallel)

# Set directories
ref_rg_dir <- file.path("/scratch/lsong36/comTCA/data/expert_range_maps",
                        "refined_range_rasters")
tz_dir <- file.path(ref_rg_dir, "tanzania")
dst_dir <- '/scratch/lsong36/comTCA/data/biodiveristy'
if(!dir.exists(dst_dir)) dir.create(dst_dir)

# Get range sizes for the global refined maps
taxon_group <- c("mammals", "amphibians", "birds", "reptiles")
fnames <- unlist(lapply(taxon_group, function(taxon){
    # Get the list of species
    fnames <- list.files(
        file.path(tz_dir, taxon), full.names = TRUE)
    if (taxon == "birds") fnames <- fnames[!grepl("_[NBR]{1}", fnames)]
    fnames
}))
message(sprintf("Get the fnames: %s.", length(fnames)))

# Save the empty csv file to receive the result
csv_path <- file.path(dst_dir, "tanzania_range_size.csv")
write.csv(data.frame(species = character(), tanzania_range_size = numeric()),
          csv_path, row.names = FALSE)
message("Save out the csv template.")

for (fname in fnames){
    system(sprintf("sbatch schedulers/run_calc_tz_range_size.sh %s", fname))
}