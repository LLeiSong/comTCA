## --------------------------------------------
## Script name: query_dsm
## Purpose of script: query ALOS DSM elevation
## dataset from the official website
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Notes: Make sure to register first.
## --------------------------------------------

# Load libraries
library(stringr)

query_dsm <- function(bbox, dst_path){
    # Adjust the ranges
    ranges <- bbox
    ranges[1:2] <- floor(ranges[1:2])
    ranges[3:4] <- ceiling(ranges[3:4])
    
    # Get the tiles for downloading
    ## longitude
    if (ranges[1] > 0 & ranges[3] >= 0) {
        message("The whole area is in East side.")
        
        start_loc <- ranges[1] %/% 5 * 5
        end_loc <- (ranges[3] %/% 5 + 1) * 5
        lon_zones <- seq(start_loc, end_loc, 5)
        lon_zones <- paste0("E", str_pad(lon_zones, 3, pad = "0"))
    } else if (ranges[1] < 0 & ranges[3] >= 0) {
        message("The whole area is across West and East side.")
        
        start_loc <- ranges[1] %/% 5 * 5
        end_loc <- (ranges[3] %/% 5 + 1) * 5
        lon_zones <- seq(start_loc, end_loc, 5)
        lon_signs <- ifelse(lon_zones >= 0, "E", "W")
        lon_zones <- paste0(lon_signs, str_pad(abs(lon_zones), 3, pad = "0"))
        rm(lon_signs)
    } else if (ranges[1] < 0 & ranges[3] < 0) {
        message("The whole area is in West side.")
        
        start_loc <- ranges[1] %/% 5 * 5
        end_loc <- (ranges[3] %/% 5 + 1) * 5
        lon_zones <- seq(start_loc, end_loc, 5)
        lon_zones <- paste0("W", str_pad(abs(lon_zones), 3, pad = "0"))
    }
    ## latitude
    if (ranges[2] > 0 & ranges[4] >= 0) {
        message("The whole area is in North side.")
        
        start_loc <- ranges[2] %/% 5 * 5
        end_loc <- (ranges[4] %/% 5 + 1) * 5
        lat_zones <- seq(start_loc, end_loc, 5)
        lat_zones <- paste0("N", str_pad(lat_zones, 3, pad = "0"))
    } else if (ranges[2] < 0 & ranges[4] >= 0) {
        message("The whole area is across South and North side.")
        
        start_loc <- ranges[2] %/% 5 * 5
        end_loc <- (ranges[4] %/% 5 + 1) * 5
        lat_zones <- seq(start_loc, end_loc, 5)
        lat_signs <- ifelse(lat_zones >= 0, "N", "S")
        lat_zones <- paste0(lat_signs, str_pad(abs(lat_zones), 3, pad = "0"))
        rm(lat_signs)
    } else if (ranges[2] < 0 & ranges[4] < 0) {
        message("The whole area is in South side.")
        
        start_loc <- ranges[2] %/% 5 * 5
        end_loc <- (ranges[4] %/% 5 + 1) * 5
        lat_zones <- seq(start_loc, end_loc, 5)
        lat_zones <- paste0("S", str_pad(abs(lat_zones), 3, pad = "0"))
    }; rm(ranges, start_loc, end_loc)
    
    zones <- sapply(1:(length(lat_zones) - 1), function(n) {
        sprintf("%s%s_%s%s", lat_zones[n], 
                lon_zones[1:(length(lon_zones) - 1)],
                lat_zones[n + 1], 
                lon_zones[2:length(lon_zones)])})
    
    # Make a list of urls
    base_url <- "https://www.eorc.jaxa.jp/ALOS/aw3d30/data/release_v2012"
    urls <- file.path(base_url, sprintf("%s.zip", zones))
    rm(zones, base_url, lat_zones, lon_zones)
    
    # Start to download
    opts <- options()
    options(timeout = 60 * 120)
    lapply(urls, function (url_dl) {
        message(sprintf("Download %s to %s", url_dl, dst_path))
        
        # Download
        fname <- file.path(dst_path, basename(url_dl))
        if (!dir.exists(gsub(".zip", "", fname))){
            tryCatch({
                download.file(url_dl, 
                              destfile = fname)
            }, error = function(e) {
                warning("Failed to download, try again.")
                file.remove(fname)
                tryCatch({
                    download.file(url_dl, 
                                  destfile = fname)
                }, error = function(e) {
                    warning("Failed to download again.")
                    file.remove(fname)})})
            
            # unzip
            if (file.exists(fname)) {
                cmd <- sprintf("unzip %s -d %s", fname, dst_path)
                try(system(cmd, intern = TRUE, ignore.stdout = TRUE))
                file.remove(fname)}
        }
    }); options(opts)
    rm(bopts, cmd, data_path, fname, url_dl, urls)
}
