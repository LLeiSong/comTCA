## --------------------------------------------
## Script name: query_habitat_elevation_redlist
## Purpose of script: query assessment from IUCN
## Red List
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Notes: A dirty way to query pdf, then parse
## information from pdf file. Please request
## API token to query in a formal way.
## --------------------------------------------

# Load packages
library(tabulizer)
library(stringr)

get_values <- function(a_row){
    # Replace Unknow n
    a_row <- gsub("Unknow n", "Unknown", a_row)
    # Check season
    has_season <- str_detect(a_row, "Breeding|Non- breeding")
    has_season_full <- str_detect(a_row, "Breeding season|Non- breeding season")
    
    # Split with space first
    items <- strsplit(a_row, " ")[[1]]
    if (has_season){
        if (has_season_full){
            season <- str_extract(a_row, "Breeding season|Non- breeding season")
            season <- gsub("- ", "-", season)
            items <- items[!str_detect(items, "Breeding|breeding|season|Non-")]
            items <- c(items[1:(length(items) - 2)], season, 
                       items[(length(items) - 1):length(items)])
        } else {
            season <- str_extract(a_row, "Breeding|Non- breeding")
            season <- gsub("- ", "-", season)
            items <- items[!str_detect(items, "Breeding|Non-")]
            items <- c(items[1:(length(items) - 2)], season, 
                       items[(length(items) - 1):length(items)])
        }
        
        last_three <- items[(length(items) - 2):length(items)]
        
        # For birds/reptiles
        first_one <- str_trim(gsub("Non-|breeding|Breeding|season", "", a_row))
        first_one <- str_trim(gsub(paste(last_three[2:3], collapse = " "), "", first_one))
        
        # first_one <- str_trim(gsub(last_three[1], "", first_one))
        # first_one <- str_trim(gsub(last_three[2], "", first_one))
        # first_one <- str_trim(gsub(last_three[3], "", first_one))
    } else {
        last_three <- items[(length(items) - 2):length(items)]
        # For mammals
        first_one <- str_trim(gsub(paste(last_three, collapse = " "), "", a_row))
    }
    
    lc_class <- strsplit(
        str_extract(first_one, "[0-9]+. [A-Z|a-z]+"), ". ")[[1]][2]
    if (lc_class == "Artificial") {
        lc_class <- strsplit(
            str_extract(first_one, "[0-9]+. [A-Z|a-z]+/[A-Z|a-z]+"), ". ")[[1]][2]
    }
    if (lc_class == "Wetlands") {
        lc_class <- paste(lc_class, "(inland)")
    }
    if (lc_class == "Root"){
        sub_class <- gsub("0. Root -> 6. ", "", first_one)
    } else sub_class <- strsplit(first_one, " - ")[[1]][2]
    
    c(lc_class, sub_class, last_three)
}

# Check if there are any fragmented rows
has_frag_row <- function(tbl){
    base_row <- str_detect(tbl, "^[0-9]+. ")
    
    if (sum(!base_row) > 0) {
        return(TRUE)
    } else return(FALSE)
}

# Merge fragmented rows back to the base row
merge_row <- function(tbl){
    base_row <- str_detect(tbl, "^[0-9]+. ")
    frag_row_id <- which(base_row != TRUE)
    
    for (id in rev(frag_row_id)){
        tbl[id - 1] <- paste(tbl[id - 1], tbl[id], sep = " ")
    }
    
    tbl[-frag_row_id]
}

parse_habitat_elevation <- function(species, data_dir, pdf_fname = "temp"){
    # Read the info
    assessment_id <- species[[1]]
    taxon_id <- species[[2]]
    sci_name <- species[[3]]
    category <- species[[4]]
    
    # Download the pdf file
    url <- sprintf("https://www.iucnredlist.org/species/pdf/%s", assessment_id)
    
    # For some species included lately, no report yet
    if (RCurl::url.exists(url)) {
        download.file(url = url, destfile = file.path(
            data_dir, sprintf("%s.pdf", pdf_fname)))
        
        # Parse the text in the pdf
        text <- extract_text(file.path(data_dir, sprintf("%s.pdf", pdf_fname)))
        text <- strsplit(text, "\n")[[1]]
    } else text <- NA
    
    # Collect table of habitats
    ids <- str_detect(
        text, paste0(
            "^[0-9]{1}. Forest|",
            "^[0-9]{1}. Savanna|",
            "^[0-9]{1}. Shrubland|",
            "^[0-9]{1}. Grassland|",
            "^0. Root -> 6. Rocky areas (eg. inland cliffs, mountain peaks)|",
            "^0. Root -> 18. Unknown|",
            "^0. Root -> 18. |",
            "^0. Root -> 6. |",
            "^5. Wetlands|",
            "^[0-9]{1}. Caves and Subterranean Habitats ",
            "(non-aquatic)|",
            "^[0-9]{2}. Artificial/Terrestrial|",
            "^[0-9]{2}. Artificial/Aquatic|",
            "^13. Marine Coastal/Supratidal|",
            "^12. Marine Intertidal|",
            "^8. Desert|",
            "^9. Marine Neritic|",
            "^10. Marine Oceanic|",
            "^11. Marine Deep Ocean Floor"))
    ids <- str_detect(
        text, "Suitable|Marginal") | ids
    ids <- which(ids)
    
    # Fix the column names
    col_names <- c("Habitat", "Sub_habitat", "Season",
                   "Suitability", "Major_Importance")
    
    if ((!identical(text, NA)) & (sum(str_detect(
        text, "Habitat Season Suitability MajorImportance?")) > 0)){
        # Move below the habitat table
        id_lower <- min(which(str_detect(
            text, "Habitat Season Suitability MajorImportance?")))
        
        ids <- ids[ids > id_lower]
        
        # Adjust the last id
        if (any(str_detect(text[max(ids) + 1:4],
                           "Suitable|Unknown|Marginal")) &
            (!any(str_detect(text[max(ids) + 1:4], "[0-9]+")))){
            ids_ext <- max(ids) + 1:4
            id_max <- ids_ext[str_detect(text[ids_ext], "Suitable|Unknown|Marginal")]
            # If the habitat table exist
            ids <- seq(min(ids), id_max)
        } else ids <- seq(min(ids), max(ids))
        
        hbt_table <- text[ids]
        
        # Merge multiple-line rows
        if (has_frag_row(hbt_table)) {
            hbt_table <- merge_row(hbt_table)
        }
        
        # Check footnote
        ftn_id <- which(str_detect(hbt_table, " © The IUCN.+MajorImportance\\?"))
        hbt_table[ftn_id] <- gsub(" © The IUCN.+MajorImportance\\?", "",
                                  hbt_table[ftn_id])
        
        # Collect columns
        col_vls <- do.call(rbind, lapply(hbt_table, get_values))
        col_vls <- data.frame(col_vls); names(col_vls) <- col_names
    } else {
        # If the habitat table does not exist
        col_vls <- data.frame(t(rep(NA, 5))); names(col_vls) <- col_names
    }
    
    # Find out elevation limits
    low_ele <- grep("Lower elevation limit", text)
    if (length(low_ele) > 0){
        low_ele <- text[low_ele]
        low_ele <- strsplit(low_ele, ": ")[[1]][2]
    } else low_ele <- NA
    high_ele <- grep("Upper elevation limit", text)
    if (length(high_ele) > 0){
        high_ele <- text[high_ele]
        high_ele <- strsplit(high_ele, ": ")[[1]][2]
    } else high_ele <- NA
    
    # Delete the file
    file.remove(file.path(data_dir, sprintf("%s.pdf", pdf_fname)))
    
    # Return
    return(list(assessment_id = assessment_id,
         taxon_id = taxon_id,
         scientific_name = sci_name,
         category = category,
         habitats = col_vls,
         lower_elevation = low_ele,
         upper_elevation = high_ele))
}
