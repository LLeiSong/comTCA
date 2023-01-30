## --------------------------------------------
## Script name: species richness calculation
## Purpose of script: calculate species richness
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Notes: Just a code chunk, cannot run outside
## of R markdown.
## --------------------------------------------

# Some reference numbers
# mammal: 170,594.24 km2 (old pre dropping EW EX: 170568.45 km2
# amphibian: 4,443.1651 km2 (old pre dropping EW EX: 4439.1885 km2)
# bird: 471064.67 km2 (even after dropping EW EX)
# reptile: 16,955.5827 km2 (old pre dropping EW EX: 16917.3988 km2)

# Calculate the richness
# Make a raster template
template <- rast(file.path(range_dir, "mammal/Acinonyx jubatus.tif")) %>% 
    crop(bry_buf) %>% terra::mask(bry_buf)
values(template) <- 0
for (txn in taxon_group){
    message(sprintf("Taxonomic group: %s", txn))
    # Read geojson
    threat_weights <- read_sf(
        file.path(species_dir, sprintf("%ss_tz.geojson", txn))) %>% 
        st_drop_geometry() %>% unique() %>% 
        mutate(threat_weight = 
            ifelse(category == "CR", 1,    # if CR then 1, if not, then ->
            ifelse(category == "EN", 1/2,  # if EN then 1/2, if not, then ->
            ifelse(category == "VU", 1/4,  # if VU then 1/4, if not, then ->
            ifelse(category == "NT", 1/8,  # if NT then 1/8, if not, then ->
            ifelse(category == "LC", 1/16, # if LC then 1/16, if not, then ->
            # if DD then 1/4 (a middle ground, conservative estimate).
            ifelse(category == "DD", 1/4, 0)))))))  # If not, then -> 0. 
    
    # Check the threat weights
    message(sprintf("Any invalid category?: %s", 
                    any(threat_weights$threat_weight == 0)))
    
    # Get the list of species
    fnames <- list.files(
        file.path(range_dir, txn), full.names = TRUE)
    if (txn == "bird") fnames <- fnames[!grepl("[1-3]{1}", fnames)]
    species <- tools::file_path_sans_ext(basename(fnames))
    
    # Check file numbers
    message(sprintf("File lengths are the same?: %s", 
                    identical(length(fnames), length(species))))
    
    # Richness
    message("Threat-weighted richness/threat-weighted endemism richness.")
    sr <- rast(template, vals = 0)
    rwsr <- rast(template, vals = 0)
    twsr <- rast(template, vals = 0)
    twrwsr <- rast(template, vals = 0)
    for (fname in fnames){
        ## Species richness----------------------------------------------------
        rg <- rast(fname) %>% crop(bry_buf) %>% terra::mask(bry_buf)
        rg <- extend(rg, ext(sr))
        rg[is.na(rg)] <- 0
        
        ## Endemism_richness---------------------------------------------------
        global_sz <- ranges_sizes %>% 
            filter(species == tools::file_path_sans_ext(basename(fname))) %>% 
            pull(ranges_size)
        rrg <- rg / global_sz
        
        ## Threat-weighted richness--------------------------------------------
        t_weight <- threat_weights %>% 
            filter(sci_name == tools::file_path_sans_ext(basename(fname))) %>% 
            pull(threat_weight)
        twrg <- rg * t_weight
        
        ## Threat-weighted richness--------------------------------------------
        twrrg <- rrg * t_weight
        
        # Accumulate species
        sr <- sr + rg
        rwsr <- rwsr + rrg
        twsr <- twsr + twrg
        twrwsr <- twrwsr + twrrg
        
        # Clean
        rm(rg, rrg, twrg, twrrg)}
    
    # Write out final layers
    writeRaster(
        sr, file.path(dst_dir, sprintf("species_richness_%s.tif", txn)))
    writeRaster(
        rwsr, file.path(dst_dir, sprintf("endemism_richness_%s.tif", txn)))
    writeRaster(
        twsr, file.path(dst_dir, sprintf("threat_weighted_richness_%s.tif", txn)))
    writeRaster(
        twrwsr, file.path(dst_dir, sprintf("threat_weighted_endemism_%s.tif", txn)))
    rm(sr, rwsr, twsr, twrwsr) # Clean up
    
    # Small-ranged species richness
    message("Small-ranged species richness/endemism richness.")
    global_median <- global_medians %>% 
        filter(taxon == txn) %>% pull(refined_median)
    # Check the number
    message(sprintf("Global median: %s", global_median))
    sr_species <- ranges_sizes %>% 
        filter(ranges_size < global_median) %>% 
        pull(species)
    sr_fnames <- fnames[which(species %in% sr_species)]
    
    smsr <- rast(template, vals = 0)
    smrwsr <- rast(template, vals = 0)
    for (fname in sr_fnames){
        rg <- rast(fname) %>% crop(bry_buf) %>% terra::mask(bry_buf)
        rg <- extend(rg, ext(smsr))
        rg[is.na(rg)] <- 0
        
        global_sz <- ranges_sizes %>% 
            filter(species == tools::file_path_sans_ext(basename(fname))) %>% 
            pull(ranges_size)
        rrg <- rg / global_sz
        
        # Accumulate species
        smsr <- smsr + rg
        smrwsr <- smrwsr + rrg
        
        # Clean up
        rm(rg, rrg)}
    
    # Write out final layers
    writeRaster(
        smsr, file.path(dst_dir, sprintf("small_range_species_richness_%s.tif", txn)))
    writeRaster(
        smrwsr, file.path(dst_dir, sprintf("small_range_endemism_richness_%s.tif", txn)))
    rm(smsr, smrwsr) # Clean up
    
    ## Local species richness--------------------------------------------------
    message("Local species richness/endemism richness.")
    lc_species <- local_species[[txn]]
    lc_fnames <- fnames[which(species %in% lc_species)]
    
    lssr <- rast(template, vals = 0)
    lsrwsr <- rast(template, vals = 0)
    for (fname in lc_fnames){
        rg <- rast(fname) %>% crop(bry_buf) %>% terra::mask(bry_buf)
        rg <- extend(rg, ext(lssr))
        rg[is.na(rg)] <- 0
        
        global_sz <- ranges_sizes %>% 
            filter(species == tools::file_path_sans_ext(basename(fname))) %>% 
            pull(ranges_size)
        rrg <- rg / global_sz
        
        # Accumulate species
        lssr <- lssr + rg
        lsrwsr <- lsrwsr + rrg
        
        # Clean up
        rm(rg, rrg)}
    
    # Write out final layers
    writeRaster(
        lssr, file.path(dst_dir, sprintf("local_species_richness_%s.tif", txn)))
    writeRaster(
        lsrwsr, file.path(dst_dir, sprintf("local_species_endemism_richness_%s.tif", txn)))
    rm(lssr, lsrwsr)
}

