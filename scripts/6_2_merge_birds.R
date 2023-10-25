# Title     : Sum up the bird ranges: resident, breeding and non-breeding
# Objective : The main function to merge the range of birds with different types.
# Created by: Lei Song
# Created on: 07/26/23

# Load libraries
library(dplyr)
library(stringr)
library(parallel)

# Set directories
ref_rg_dir <- file.path("/scratch/lsong36/comTCA/data/expert_range_maps",
                        "refined_range_rasters")
global_dir <- file.path(ref_rg_dir, "global/birds")
outside_dir <- file.path(ref_rg_dir, "outside/birds")
tz_dir <- file.path(ref_rg_dir, "tanzania/birds")

# Get file catalog
global_fnames <- list.files(global_dir, full.names = TRUE) %>% 
    data.frame(fname = .) %>% 
    mutate(species = gsub("_[NRB]{1}.tif", "", basename(fname)))

# Define function
mosaic_layers <- function(fnames, dst_path){
    options_calc <- '--NoDataValue=0 --hideNoData --co compress=lzw --type Byte'
    if (length(fnames) == 2){
        eq_string <- sprintf("logical_or(A==1,B==1)")
        command <- sprintf(
            'gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
            fnames[1], fnames[2], dst_path,
            eq_string, options_calc)
    } else {
        eq_string <- sprintf("logical_or(logical_or(A==1,B==1),C==1)")
        command <- sprintf(
            'gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
            fnames[1], fnames[2], fnames[3], dst_path,
            eq_string, options_calc)
    }
    system(command)
}

# Start to mosaic
mclapply(unique(global_fnames$species), function(sps){
    # Subset
    global_this <- global_fnames %>% filter(species == sps)
    
    # Get the list of fnames
    global_fns <- global_this$fname
    outside_fns <- file.path(outside_dir, basename(global_fns))
    tz_fns <- file.path(tz_dir, basename(global_fns))
    
    # Mosaic
    dst_fname <- sprintf("%s.tif", sps)
    if (length(global_fns) == 1){
        message("Just one layer found, so copy to new one.")
        file.copy(global_fns, file.path(global_dir, dst_fname))
        file.copy(outside_fns, file.path(outside_dir, dst_fname))
        file.copy(tz_fns, file.path(tz_dir, dst_fname))
    } else {
        message("More than one layers found, mosaic them.")
        mosaic_layers(global_fns, file.path(global_dir, dst_fname))
        mosaic_layers(outside_fns, file.path(outside_dir, dst_fname))
        mosaic_layers(tz_fns, file.path(tz_dir, dst_fname))
    }
}, mc.cores = 20)