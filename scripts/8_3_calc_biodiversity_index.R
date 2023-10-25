# Title     : Calculate biodiversity index
# Objective : The main function to calculate biodiversity index that will be
#             used for the final calculation.
# Created by: Lei Song
# Created on: 09/09/23
# Usage: python scripts/calc_species_richness.py --txn birds

library(sf)
library(here)
library(terra)
sf_use_s2(FALSE)

data_dir <- here("data")
bio_data_dir <- here("data/biodiversity")

bry <- read_sf(here("data/geoms/tanzania_coarse_bry.geojson"))
bbox_ll <- bry %>% st_bbox() %>% 
st_as_sfc() %>% st_buffer(1, endCapStyle = "SQUARE")

# Read metrics
bhi <- rast(file.path(bio_data_dir, "BILBI_P_BHIv2_Habitat_2020.tif")) %>% 
    crop(bbox_ll)
msa <- rast(file.path(bio_data_dir, "TerrestrialMSA_2015_World.tif")) %>% 
    crop(bbox_ll) %>% resample(bhi)
bii <- rast(file.path(bio_data_dir, "lbii.asc")) %>% 
    crop(bbox_ll)

# save out
bhi_msa_bii <- c(bhi, msa, bii)
names(bhi_msa_bii) <- c("BHI", "MSA", "BII")
writeRaster(bhi_msa_bii, file.path(bio_data_dir, "bhi_msa_bii.tif"))

# Read files
bhi_msa_bii <- rast(file.path(bio_data_dir, "bhi_msa_bii.tif"))

# Calculate species component
normalize <- function(x, robust = TRUE) {
    if (robust) {
        stretch(x, minv = 0, maxv = 1, minq = 0.01, maxq = 0.99)
    } else {
        (x - minmax(x)[1]) / (minmax(x)[2] - minmax(x)[1])
    }
}

fnames <- list.files(bio_data_dir, pattern = "*weighted*", full.names = TRUE)
richness_fns <- fnames[!grepl("rarity", fnames)]

richness <- do.call(
    c, lapply(richness_fns, function(fn){
        normalize(rast(fn) %>% mask(bry))})) %>% 
    mean() %>% normalize(FALSE)

endemism_fns <- fnames[grepl("rarity", fnames)]
# Use robust scaling to avoid the impacts of outliers
endemism <- do.call(
    c, lapply(endemism_fns, function(fn){
        normalize(sqrt(sqrt(rast(fn) %>% mask(bry))))})) %>% 
    mean() %>% normalize(FALSE)

# Geometric mean for species richness and endemism richness
# (the Nth root of the product of N numbers - 
# in the case of combining two layers x and y, 
# a pixels value would be the square root of x*y)
S <- sqrt(richness * endemism)

# Calculate ecosystems component
# Arithmetic mean for MSA and BII
# (the Nth root of the product of N numbers - 
# in the case of combining two layers x and y, 
# a pixels value would be the square root of x*y)
E <- mean(bhi_msa_bii[["MSA"]], bhi_msa_bii[["BII"]])

# Calculate biodiversity indicators
b <- resample(S, E) + E
c <- bhi_msa_bii[["BHI"]]
BIp <- b * c

# Save out
# NOTE: the values outside of Tanzania is not valid
# because species components does not cover these areas.
names(BIp) <- "BIp"
writeRaster(BIp, file.path(bio_data_dir, "BIp.tif"))
