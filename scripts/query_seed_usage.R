## --------------------------------------------
## Script name: query_seed_usage
## Purpose of script: query seed usage from the
## NSCA survey database.
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Notes: You have to download the survey data.
## --------------------------------------------

# Util function to query seed usage
query_seed_usage <- function(src_dir, season, crop){
    # Load pacakge
    library(tidyr)
    library(stringr)
    library(haven)
    library(dplyr)
    
    # Check inputs
    checkmate::check_choice(season, choices = c("long", "short"))
    # Only because I just test these two types
    checkmate::check_choice(crop, choices = c("maize", "paddy"))
    
    # Target files
    file_index <- ifelse(season == "short", "81", "82")
    
    # Query full table
    ## Assume both local and improved be improved
    seeds <- read_dta(file.path(
        src_dir, sprintf("R0%s_%s_RAINY_SEASON.dta", 
                         file_index, toupper(season)))) %>% 
        select(region, sprintf("q%s1c4", file_index), 
               sprintf("q%s1c11", file_index)) %>% 
        rename(Region = region, 
               Crop = sprintf("q%s1c4", file_index), 
               Seed = sprintf("q%s1c11", file_index)) %>% 
        mutate(Region = as_factor(Region),
               Crop = as_factor(Crop),
               Seed = as_factor(Seed)) %>% 
        mutate(Seed = ifelse(
            Seed == "Local seed", "Local seed",
            "Improved seed"))
    
    # Get Maize and Paddy
    seeds %>% 
        filter(Crop == str_to_title(crop)) %>% 
        group_by(Region, Seed) %>% 
        summarise(value = n()) %>% 
        ungroup() %>% 
        pivot_wider(names_from = Seed,
                    values_fill = 0) %>% 
        mutate(`Seed number` = `Improved seed` + `Local seed`) %>% 
        mutate(`Improved seed` = `Improved seed` / `Seed number`,
               `Local seed` = `Local seed` / `Seed number`,
               Season = season, Crop = crop) %>% 
        select(Region, Season, Crop, `Seed number`, 
               `Improved seed`, `Local seed`)
}
