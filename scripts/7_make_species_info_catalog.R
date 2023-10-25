# Load libraries
library(sf)
library(dplyr)
library(tidyr)
library(stringr)

# Set directories
data_dir <- "/scratch/lsong36/comTCA/data"
bio_dir <- "/scratch/lsong36/comTCA/data/biodiveristy"

# Read range info
out_range_size <- read.csv(file.path(bio_dir, "outside_range_size.csv"))
tz_range_size <- read.csv(file.path(bio_dir, "tanzania_range_size.csv"))
gl_range_size <- read.csv(file.path(bio_dir, "global_range_size.csv"))

# Gather them together
range_size <- left_join(out_range_size, tz_range_size, by = "species")
range_size <- left_join(range_size, gl_range_size, by = "species")

# Calculate the areas, use degree directly
global_map <- rast(file.path(
    data_dir, "habitats", "habitat_outside_lvl2.tif"))
tz_map <- rast(file.path(data_dir, "habitats/habitat_tz_refined_final.tif"))

range_size <- range_size %>% 
    mutate(outside_range_size = outside_range_size * res(global_map)[1]^2,
           tanzania_range_size = tanzania_range_size * res(tz_map)[1]^2,
           range_size = range_size * res(global_map)[1]^2) %>% 
    mutate(refined_range_size = outside_range_size + tanzania_range_size)

write.csv(range_size, file.path(bio_dir, "range_size.csv"), row.names = FALSE)

# Make the species catalog
fnames <- list.files(file.path(data_dir, "expert_range_maps/clean_range_maps"),
                     pattern = ".geojson", full.names = TRUE)
species_category <- lapply(fnames, function(fname){
    taxon <- str_extract(fname, "mammals|amphibians|reptiles|birds")
    ctgs <- st_read(fname) %>% st_drop_geometry() %>% 
        select(sci_name, category) %>% unique() %>% 
        mutate(taxon = taxon)
    
    # Double check if there are conflicts on category
    abnormal_num <- ctgs %>% 
        group_by(sci_name, category) %>% summarise(n = n()) %>% 
        filter(n > 1) %>% nrow()
    if (abnormal_num != 0) warning("Conflicts found in species category!")
    
    # Return
    ctgs
}) %>% bind_rows()

species_info <- left_join(range_size, species_category, 
                          by = c("species" = "sci_name"))
write.csv(species_info, file.path(bio_dir, "species_info.csv"), row.names = FALSE)
