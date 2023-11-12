# Title     : Select species to calculate species-level habitat suitability
# Objective : The main function to select core species to calculate the habitat 
#             suitability. 
# Created by: Lei Song
# Created on: 09/09/23

# Load libraries
library(sf)
library(terra)
library(here)
library(dplyr)
library(rredlist)
library(dplyr)
library(tidyr)

# Paths
data_dir <- here("data")
conn_data_dir <- here("data/connectivity")
if (!dir.exists(conn_data_dir)) dir.create(conn_data_dir)

bry <- read_sf(here("data/geoms/tanzania_coarse_bry.geojson"))

# select some core grazing species, herbivores
# Some references:
# elephant (Loxodonta africana), buffalo (Syncerus caffer), 
# zebra (Equus quagga), wildebeest (Connochaetes taurinus), 
# hartebeest (Alcelaphus buselaphus), topi (Damaliscus korrigum), 
# impala (Aepyceros melampus) and Thomson's gazelle (Eudorcas thomsonii)
# https://royalsocietypublishing.org/doi/10.1098/rstb.2015.0314
# African elephant (Loxodonta africana), reticulated giraffe (Giraffa reticulata), 
# plains zebra (Equus quagga), and Grevy’s zebra (Equus grevyi)
# https://doi.org/10.1007/s10980-021-01232-8
# Remove small mammals: up to 5 kg. Assume they can't or not necessarily move 
# for a long distance.

# The logic here is first select some qualified species
# Second, check the home range in Tanzania
# Remove the not widely distributed ones.
all_mammals <- st_read(
    file.path("data/expert_range_maps/clean_range_maps",
              "mammals_tz_relevant.geojson")) %>% 
    st_drop_geometry() %>% 
    select(sci_name, category) %>% unique()

# Filter the non-volant, non-immigrant herbivores
## Remove some orders less possibly migrate for long distance
## mainly based on body size (e.g. rodentia) and living style (e.g. Primates)
all_mammals <- lapply(1:nrow(all_mammals), function(n){
    sps <- all_mammals %>% slice(n)
    odr <- rl_search(sps$sci_name)$result$order
    sps %>% mutate(order = odr)
}) %>% bind_rows()
all_mammals <- all_mammals %>% 
    filter(!order %in% c(
        "CHIROPTERA", "CARNIVORA", "RODENTIA", "EULIPOTYPHLA",
        "HYRACOIDEA", "AFROSORICIDA", "TUBULIDENTATA", 
        "MACROSCELIDEA", "LAGOMORPHA", "PHOLIDOTA", "PRIMATES"))

# Filter habitat area, focusing on large habitat range
## Get refined habitat maps
fnames <- file.path("/Users/leisong/Dropbox/research/comTCA/data",
              "expert_range_maps/refined_range_rasters/tanzania_100/mammals",
              sprintf("%s.tif", gsub(" ", "_", all_mammals$sci_name)))
fnames_dc <- list.files(
    file.path("/Users/leisong/Dropbox/research/comTCA/data",
              "expert_range_maps/refined_range_rasters/tanzania_100/mammals"),
    full.names = TRUE)
fnames <- intersect(fnames, fnames_dc); rm(fnames_dc)

## Filter species with desired area
## Assume if a habitat is too small (5%), then better to set up protected area to
## protect them.
tz_area <- rast(fnames[1])
tz_area <- subst(tz_area, NA, 0)
tz_area <- mask(tz_area, bry)
tz_area <- ncell(tz_area) - global(tz_area, "isNA")[[1]]

species_select <- lapply(fnames, function(fn){
    rst <- rast(fn)
    area <- ncell(rst) - global(rst, "isNA")[[1]]
    if (round(area / tz_area, 2) > 0.05){
        gsub("_", " ", gsub(".tif", "", basename(fn)))
    } else {
        NULL
    }
})

species_select <- unlist(species_select[!sapply(species_select, is.null)])
all_mammals <- all_mammals %>% filter(sci_name %in% species_select)

# Remove non-migrant species
species_to_mv <- c(
    "Otolemur garnettii", "Papio cynocephalus", 
    "Galago senegalensis", "Phacochoerus africanus",
    "Chlorocebus pygerythrus", "Otolemur crassicaudatus",
    "Papio anubis")
all_mammals <- all_mammals %>% filter(!sci_name %in% species_to_mv)

# Remove species not in Tanzania (due to the coarse resolution of boundary)
species_select <- lapply(all_mammals$sci_name, function(sps_nm){
    if ("Tanzania, United Republic of" %in% rl_occ_country(sps_nm)$result$country){
        sps_nm
    } else NULL
})
species_select <- unlist(species_select[!sapply(species_select, is.null)])
all_mammals <- all_mammals %>% filter(sci_name %in% species_select)

# Final filter
## Use all available info: habitat range size and threats
fnames <- file.path("/Users/leisong/Dropbox/research/comTCA/data",
                    "expert_range_maps/refined_range_rasters/tanzania_100/mammals",
                    sprintf("%s.tif", gsub(" ", "_", all_mammals$sci_name)))
fnames_dc <- list.files(
    file.path("/Users/leisong/Dropbox/research/comTCA/data",
              "expert_range_maps/refined_range_rasters/tanzania_100/mammals"),
    full.names = TRUE)
fnames <- intersect(fnames, fnames_dc); rm(fnames_dc)

tz_area <- rast(fnames[1])
tz_area <- subst(tz_area, NA, 0)
tz_area <- mask(tz_area, bry)
tz_area <- ncell(tz_area) - global(tz_area, "isNA")[[1]]

species_info <- lapply(fnames, function(fn){
    rst <- rast(fn)
    area <- ncell(rst) - global(rst, "isNA")[[1]]
    sps_nm <- gsub("_", " ", gsub(".tif", "", basename(fn)))
    pup_trend <- rl_search(sps_nm)$result$population_trend
    # 1, 2, 4, 8, 11
    threats <- rl_threats(sps_nm)$result
    threats <- threats %>% separate(code, c("major", "minor")) %>% 
        filter(major %in% c(1, 2, 4, 8, 11))
    data.frame(sci_name = sps_nm,
               abundance = area / tz_area,
               pup_trend = pup_trend,
               threats = ifelse(nrow(threats) > 0, "movement_relevant", "others"))
}) %>% bind_rows()

# Has threats relevant to movement and population is decreasing.
species_select <- species_info %>% 
    filter(threats == "movement_relevant")
all_mammals <- all_mammals %>% filter(sci_name %in% species_select$sci_name)
all_mammals <- left_join(all_mammals, species_select, by = "sci_name") %>% 
    filter(!(category == "LC" & pup_trend != "Decreasing"))

write.csv(all_mammals, "data/connectivity/select_species.csv", row.names = FALSE)
