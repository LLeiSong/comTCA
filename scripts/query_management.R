## --------------------------------------------
## Script name: query_management
## Purpose of script: query management from the
## NSCA survey database
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Notes: You have to download the survey data.
## --------------------------------------------

# Util function to query if irrigation or fertilizer is applied
query_management <- function(src_dir, season, crop, type){
    # Load package
    library(tidyr)
    library(stringr)
    library(haven)
    library(dplyr)
    
    # Check inputs
    checkmate::check_choice(season, choices = c("long", "short"))
    # Only because I just test these types
    checkmate::check_choice(crop, choices = c("maize", "paddy"))
    checkmate::check_choice(type, choices = c("irrigation", "fertilizer"))
    
    # Target files
    file_index <- ifelse(season == "short", "81", "82")
    
    # Column index
    col_index <- ifelse(type == "irrigation", "8a", 16)
    
    # Query full table
    ## Assume both local and improved be improved
    query_table <- read_dta(file.path(
        src_dir, sprintf("R0%s_%s_RAINY_SEASON.dta", 
                         file_index, toupper(season)))) %>% 
        select(region, sprintf("q%s1c4", file_index), 
               sprintf("q%s1c%s", file_index, col_index)) %>% 
        rename(Region = region, 
               Crop = sprintf("q%s1c4", file_index), 
               Management = sprintf(
                   "q%s1c%s", file_index, col_index)) %>% 
        mutate(Region = as_factor(Region),
               Crop = as_factor(Crop),
               Management = as_factor(Management))
    
    # Get Maize and Paddy
    query_table %>% 
        filter(Crop == str_to_title(crop)) %>% 
        group_by(Region, Management) %>% 
        summarise(value = n()) %>% 
        ungroup() %>% 
        pivot_wider(names_from = Management,
                    values_fill = 0) %>% 
        mutate(Percent = Yes / (Yes + No),
               Season = season, Crop = crop, 
               Management = str_to_title(type)) %>% 
        select(Region, Season, Crop, Management, Percent)
}
