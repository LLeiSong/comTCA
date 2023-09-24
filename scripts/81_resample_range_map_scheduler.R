# Load libraries
library(terra)
library(dplyr)
library(tidyr)
library(stringr)
library(parallel)

# Set directories
hbt_dir <- "/scratch/lsong36/comTCA/data/habitats"
tz_range_dir <- "/scratch/lsong36/comTCA/data/expert_range_maps/refined_range_rasters/tanzania"
tz_range_100_dir <- "/scratch/lsong36/comTCA/data/expert_range_maps/refined_range_rasters/tanzania_100"
if (!dir.exists(tz_range_100_dir)) dir.create(tz_range_100_dir, recursive = TRUE)
sub_dir <- file.path(tz_range_100_dir, c("mammals", "amphibians", "reptiles", "birds"))
for (dr in sub_dir) if (!dir.exists(dr)) dir.create(dr, recursive = TRUE)

# # Make the template, one time thing
# hbt <- file.path(hbt_dir, "iucn_habitatclassification_composite_lvl2_ver004.tif")
# hbt_hight <- file.path(hbt_dir, "habitat_tz_refined_final.tif")
# hbt_low <- file.path(hbt_dir, "habitat_tz_low.tif") # this would be the template
# 
# options_warp <- paste0('-multi -wo NUM_THREADS=ALL_CPUS ',
#                        '-co TILED=YES -co compress=lzw -co BIGTIFF=IF_NEEDED')
# te <- ext(rast(hbt_hight))
# command <- sprintf(
#   paste0("gdalwarp -te %s %s %s %s %s %s %s"),
#   te[1], te[3], te[2], te[4], 
#   options_warp, hbt, hbt_low)
# system(command)

# Get the file list
fnames <- list.files(tz_range_dir, recursive = TRUE, full.names = TRUE)
fnames <- fnames[!grepl("_[NBR]{1}", fnames)]

for(fname in fnames){
    fname_to <- gsub("tanzania", "tanzania_100", fname)
    system(sprintf("sbatch schedulers/run_resample_range_map.sh %s %s", 
                   fname, fname_to))
}