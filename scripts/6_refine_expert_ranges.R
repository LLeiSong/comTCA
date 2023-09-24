# Title     : The main function for refining expert ranges
# Objective : The main function to actually do the refining. This is the downstream
#           : function for refine_expert_ranges_local.R
# Created by: Lei Song
# Created on: 07/26/23

# Load libraries
library(sf)
library(terra)
library(optparse)
library(dplyr)
library(tidyr)
library(stringr)

# Parse inline arguments
option_list <- list(
    make_option(c("-s", "--species"), 
                action = "store", type = 'character',
                help = "The species to process."),
    make_option(c("-t", "--taxon"), 
                action = "store", type = 'character',
                help = "The taxon of the species."),
    make_option(c("-r", "--range"), 
                action = "store", type = 'character',
                help = "The range, global or tanzania."))
opt <- parse_args(OptionParser(option_list = option_list))
species <- opt$species
taxon <- opt$taxon
range <- opt$range
fname <- gsub(" ", "_", species)

# Set directories
data_dir <- "/scratch/lsong36/comTCA/data"
range_dir <- file.path(
    data_dir, "expert_range_maps/clean_range_rasters", range, taxon)
dst_dir <- file.path(
    data_dir, "expert_range_maps/refined_range_rasters", range, taxon)
if (!dir.exists(dst_dir)) dir.create(dst_dir, recursive = TRUE)

# Save extra layers only outside of Tanzania
## This will be used for comparison in supplementary materials.
if (range == "global"){
    outside_dir <- file.path(
        data_dir, "expert_range_maps/refined_range_rasters", "outside", taxon)
    if (!dir.exists(outside_dir)) dir.create(outside_dir, recursive = TRUE)
}

# Define a small util function
loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])}

# Load layers and species information
if (range == "global"){
    habitats_map <- file.path(
        data_dir, "habitats", "iucn_habitatclassification_composite_lvl2_ver004.tif")
    outside_msk <- file.path(data_dir, "habitats/habitat_outside_lvl2.tif")
    elevation_map <- file.path(data_dir, "elevation/elevation_outside_resample.tif")
} else {
    habitats_map <- file.path(data_dir, "habitats/habitat_tz_refined_final.tif")
    elevation_map <- file.path(data_dir, "elevation/elevation_resample.tif")
}

if (taxon != "birds"){
    range_map <- file.path(range_dir, sprintf("%s.tif", fname))
    info <- loadRData(file.path(data_dir, "expert_range_maps/species_info", 
                                sprintf("%s_info.rda", taxon)))
    id <- which(sapply(info, function(x) x$name) == species)
    info <- info[[id]]
    
    # Load information
    ## No habitat info for some reptiles, then not refine them
    if (length(info$habitats$code) != 0){
        habitats <- info$habitats %>% 
            filter(!str_detect(code, "9.")) %>% # Remove Coral Reef first
            separate(code, c("habitat", "sub_habitat"), fill = "right") %>% 
            mutate(sub_habitat = ifelse(is.na(sub_habitat), 0, sub_habitat)) %>% 
            mutate(code = as.integer(habitat) * 100 + as.integer(sub_habitat)) %>% 
            filter(code < 1500) %>% # remove artificial aquatic and other not mapped
            filter(!(code >= 700 & code < 800)) %>% # remove cave habitats
            filter(!(code >= 900 & code < 1400)) # remove marine
    }
    
    ele_high <- info$elevation_upper
    ele_low <- info$elevation_lower
    if (is.na(ele_high) || is.null(ele_high)) ele_high <- 32767
    if (is.na(ele_low) || is.null(ele_low)) ele_low <- -32767
    
    # Define conditional operation
    ele_string <- sprintf("logical_and(B>=%s, B<=%s)", ele_low, ele_high)
    if (length(info$habitats$code) == 0){
        hbt_string <- "(A!=0)"
    } else if (length(habitats$code) == 0) {
        hbt_string <- "(A!=0)"
    } else if (length(habitats$code) == 1){
        hbt_string <- sprintf("(A==%s)", habitats$code)
    } else if (length(habitats$code) == 2){
        hbt_string <- sprintf("logical_or(A==%s,A==%s)", 
                              habitats$code[1], habitats$code[2])
    } else{
        hbt_string <- sprintf("logical_or(A==%s,A==%s)", 
                              habitats$code[1], habitats$code[2])
        for (i in 3:length(habitats$code)){
            hbt_string <- sprintf("logical_or(%s,A==%s)", hbt_string, habitats$code[i])
        }
    }
    
    eq_string <- sprintf("%s * %s * (C==1)", hbt_string, ele_string)
    options_calc <- '--NoDataValue=0 --hideNoData --co compress=lzw --type Byte'
    command <- sprintf('gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
                       habitats_map, elevation_map, range_map, 
                       file.path(dst_dir, sprintf("%s.tif", fname)),
                       eq_string, options_calc)
    system(command)
    
    if (range == "global") {
        Sys.sleep(30)
        eq_string <- "(A!=0) * B"
        command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                           outside_msk, file.path(dst_dir, sprintf("%s.tif", fname)), 
                           file.path(outside_dir, sprintf("%s.tif", fname)),
                           eq_string, options_calc)
        system(command)
    }
} else {
    range_maps <- list.files(range_dir, pattern = fname, full.names = TRUE)
    info <- loadRData(file.path(data_dir, "expert_range_maps/species_info", 
                                sprintf("%s_info.rda", taxon)))
    id <- which(sapply(info, function(x) x$name) == species)
    info <- info[[id]]
    
    for (range_map in range_maps){
        # Get the season
        this_season <- str_extract(range_map, "[NRB]{1}.")
        this_season <- ifelse(this_season == "R.", "Resident", 
                              ifelse(this_season == "B.", "Breeding Season", 
                                     "Non-Breeding Season"))
        
        # Load information
        ## No habitat info for some reptiles, then not refine them
        if (length(info$habitats$code) != 0){
            habitats <- info$habitats %>% 
                filter(!str_detect(code, "9.")) %>% # Remove Coral Reef first
                separate(code, c("habitat", "sub_habitat"), fill = "right") %>% 
                mutate(sub_habitat = ifelse(is.na(sub_habitat), 0, sub_habitat)) %>% 
                mutate(code = as.integer(habitat) * 100 + as.integer(sub_habitat)) %>% 
                filter(code < 1500) %>% # remove artificial aquatic and other not mapped
                filter(!(code >= 700 & code < 800)) %>% # remove cave habitats
                filter(!(code >= 900 & code < 1400)) %>%  # remove marine
                # subset the season, two species have Unknown, but marginal habitat
                # so ignore.
                filter(season == this_season | is.na(season))
        }
        
        ele_high <- info$elevation_upper
        ele_low <- info$elevation_lower
        if (is.na(ele_high) || is.null(ele_high)) ele_high <- 32767
        if (is.na(ele_low) || is.null(ele_low)) ele_low <- -32767
        
        # Define conditional operation
        ele_string <- sprintf("logical_and(B>=%s, B<=%s)", ele_low, ele_high)
        if (length(info$habitats$code) == 0){
            hbt_string <- "(A!=0)"
        } else if (length(habitats$code) == 0) {
            hbt_string <- "(A!=0)"
        } else if (length(habitats$code) == 1){
            hbt_string <- sprintf("(A==%s)", habitats$code)
        } else if (length(habitats$code) == 2){
            hbt_string <- sprintf("logical_or(A==%s,A==%s)", 
                                  habitats$code[1], habitats$code[2])
        } else{
            hbt_string <- sprintf("logical_or(A==%s,A==%s)", 
                                  habitats$code[1], habitats$code[2])
            for (i in 3:length(habitats$code)){
                hbt_string <- sprintf("logical_or(%s,A==%s)", hbt_string, habitats$code[i])
            }
        }
        
        eq_string <- sprintf("%s * %s * (C==1)", hbt_string, ele_string)
        options_calc <- '--NoDataValue=0 --hideNoData --co compress=lzw --type Byte'
        command <- sprintf('gdal_calc.py -A %s -B %s -C %s --outfile=%s --calc="%s" %s',
                           habitats_map, elevation_map, range_map, 
                           file.path(dst_dir, basename(range_map)),
                           eq_string, options_calc)
        system(command)
        
        if (range == "global") {
            Sys.sleep(30)
            eq_string <- "(A!=0) * B"
            command <- sprintf('gdal_calc.py -A %s -B %s --outfile=%s --calc="%s" %s',
                               outside_msk, file.path(dst_dir, basename(range_map)), 
                               file.path(outside_dir, basename(range_map)),
                               eq_string, options_calc)
            system(command)
        }
    }
}
