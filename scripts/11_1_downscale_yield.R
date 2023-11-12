## --------------------------------------------
## Script name: downscale_yield
## Purpose of script: downscale the calibrated SPAM 2017 gridded yield to 1 km
## Author: Lei Song
## Date Created: 2023-10-11
##
## Copyright (c) Lei Song, 2023
## Email: lsong@ucsb.edu
## --------------------------------------------
## Variables:
## MODIS LST, WorldClim prec, and soil
## 
## --------------------------------------------

# Libraries
library(glmnet)
library(ranger)
library(sf)
library(terra)
library(dplyr)

# Set directories
data_dir <- file.path("data/agriculture")
spam_yield_dir <- file.path(data_dir, 'SPAM', "spam2017v2r1_ssa_yield.geotiff")

## Country boundary
bry <- st_read("data/geoms/mainland_tanzania.geojson") %>% select()

## Rainfed Yield from SPAM
crops <- c("MAIZ", "RICE", "SORG", "CASS", "BEAN")
yields <- do.call(c, lapply(crops, function(crp) {
    fnames <- file.path(
        spam_yield_dir, sprintf("spam2017V2r1_SSA_Y_%s_R.tif", crp))
    rast(fnames) %>% crop(bry) %>% mask(bry)}))
names(yields) <- crops

# convert: kg/ha to tons/ha and set 0 to NA
yields <- yields / 1000
yields[yields == 0] <- NA
writeRaster(yields, file.path(data_dir, "yields_10km.tif"), overwrite = TRUE)

# MODIS ET and LST -----------------
## Process in GEE
### 10 km grid
fnames <- file.path(
    data_dir, "modis",
    sprintf("modis_%s_monthly_mean_2010_2020_10km.tif",
            c("et", "tmean", "tmax", "tmin", "gdd")))
modis_10km <- do.call(c, lapply(fnames, rast))
names(modis_10km) <- lapply(
    c("et", "tmean", "tmax", "tmin", "gdd"), function(item){
    sprintf("%s%s", item, 1:12)}) %>% unlist()

### 1 km grid
fnames <- file.path(
    data_dir, "modis",
    sprintf("modis_%s_monthly_mean_2010_2020_1km.tif",
            c("et", "tmean", "tmax", "tmin", "gdd")))
modis_1km <- do.call(c, lapply(fnames, rast))
names(modis_1km) <- lapply(
    c("et", "tmean", "tmax", "tmin", "gdd"), function(item){
        sprintf("%s%s", item, 1:12)}) %>% unlist()

# WorldClim precipitation related --------------
## Download from WorldClim website
fnames <- file.path(
    data_dir, "WorldClim",
    sprintf("%s_current.tif", c("prec", "srad", "vapr")))
wc_1km <- do.call(c, lapply(fnames, rast))

# Soil----------------
fnames <- list.files(file.path(data_dir, "SoilGrid250m"),
                     pattern = ".tif", full.names = TRUE, recursive = TRUE)
soils <- do.call(c, lapply(fnames, rast))

# Save out
vars_10km <- c(resample(modis_10km, yields), resample(wc_1km, yields),
               resample(soils, yields))
writeRaster(vars_10km, file.path(data_dir, "variables_10km.tif"), overwrite = TRUE)
vars_1km <- c(resample(modis_1km, wc_1km), wc_1km,
              resample(soils, wc_1km))
writeRaster(vars_1km, file.path(data_dir, "variables_1km.tif"), overwrite = TRUE)

# Clean
rm(wc_1km, soils, fnames); gc()

# Fit crop yield downscaling models
## Use Lasso to select features
## Fit random forest with overfitting set.
downscale_yield <- function(yield, variables, variables_to){
    # Get values
    vals <- values(c(yield, variables)) %>% 
        data.frame() %>% na.omit()
    names(vals)[1] <- "yield"
    
    # Fit lasso to select features
    lambdas <- 10^seq(2, -3, by = -.1)
    lasso_reg <- cv.glmnet(as.matrix(vals[, -1]), vals[, 1], alpha = 1, 
                           lambda = lambdas, standardize = TRUE, nfolds = 5)
    lasso_model <- glmnet(vals[, -1], vals[, 1], alpha = 1, 
                          lambda = lasso_reg$lambda.min, standardize = TRUE)
    
    coefs <- lasso_model$beta[, 1]
    features <- names(coefs)[abs(coefs) > 0.001]
    message("Select features using Lasso")
    
    # Fit quantile regression forest
    vals <- vals %>% select(all_of(c(names(vals)[1], features)))
    set.seed(123)
    rf_model <- ranger(yield ~ ., data = vals, quantreg = TRUE)
    save(rf_model, file = file.path(data_dir, sprintf("qrf_%s.rda", names(yield))))
    message("Fit the Quantile Regression Forest model")
    
    # Make prediction for training variables
    pred_qt <- predict(rf_model, vals[, -1], type = "quantiles",
                       quantiles = seq(0.05, 0.95, 0.01))
    pred_qt <- cbind(vals['yield'], pred_qt$predictions) %>% data.frame()
    fname <- file.path(data_dir, sprintf("yield_05_95_%s_R_1km.csv", names(yield)))
    write.csv(pred_qt, fname, row.names = FALSE)
    
    names(variables_to) <- gsub("-", ".", names(variables_to))
    yield_qt <- predict(
        subset(variables_to, features), rf_model, 
        type = "quantiles", quantiles = seq(0.05, 0.95, 0.01), na.rm = TRUE)
    
    fname <- file.path(data_dir, sprintf("yield_05_95_%s_R_1km.tif", names(yield)))
    writeRaster(yield_qt, fname)
    message("Write out the downscaled yield")
}

# Loop over crops
for (crp in crops){
    message(sprintf("Downscale %s yield", crp))
    yield <- yields[crp]
    downscale_yield(yield, vars_10km, vars_1km)
}

# Visualize yield downscaling models
library(ggplot2)
library(ggpubr)
fig_yield_qrf <- function(crp){
    fname <- file.path(data_dir, sprintf("yield_05_95_%s_R_1km.csv", crp))
    vals <- read.csv(fname)
    load(file.path(data_dir, sprintf("qrf_%s.rda", crp)))
    
    # Make figure
    mse_text <- sprintf("OOB prediction error (MSE): %s", round(rf_model$prediction.error, 2))
    r_text <- sprintf("R squared (OOB): %s", round(rf_model$r.squared, 2))
    pd <- position_dodge(0.1)
    ggplot(vals, aes(x = yield)) + 
        geom_errorbar(
            aes(ymin = `quantile..0.25`, ymax = `quantile..0.75`, 
                color = "Predicted 50% interval"), width = 0, position = pd) +
        geom_abline(color = "grey", linetype = 'dashed') +
        geom_point(aes(y = `quantile..0.5`, color = "Predicted median"), 
                   position = pd, size = 0.3) +
        scale_color_manual(
            "", values = c("Predicted median" = 'black',
                           'Predicted 50% interval' = 'lightblue')) +
        annotate(geom = "text", x = 0, 
                 y = max(apply(vals, MARGIN = 2, FUN = max)[-1]) * 0.95, 
                 label = mse_text, hjust = 0, size = 3) +
        annotate(geom = "text", x = 0, 
                 y = max(apply(vals, MARGIN = 2, FUN = max)[-1]) * 0.90, 
                 label = r_text, hjust = 0, size = 3) +
        xlab("Actual yield (t/ha)") +
        ylab("Predicted yield (t/ha)") +
        theme_pubclean() +
        theme(legend.position = "top",
              legend.text = element_text(size = 12),
              text = element_text(size = 12))
}

fig_list <- lapply(crops, fig_yield_qrf)
ggarrange(
    plotlist = fig_list,
    labels = c("A", "B", "C", "D", "E"),
    ncol = 3, nrow = 2, common.legend = TRUE, legend = "bottom")
ggsave("figures/S2_eval_yield_downscale.png", width = 8, height = 6, bg = "white")
