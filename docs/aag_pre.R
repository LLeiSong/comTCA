library(terra)
library(sf)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tidyterra)
library(stringr)
library(dplyr)
library(colorspace)
library(stars)
library(rnaturalearth)
library(ggfx)
library(stringr)
library(tidyverse)
library(ggsci) # for paper color palette
library(ggnewscale) # help to have multiple same scales
library(ggtext) # for rich text
library(ggsflabel)
library(ggthemes)
library(gganimate)
library(tmap)
sf_use_s2(FALSE)

# Paths
agro_dir <- "data/agriculture"
lc_dir <- "data/landcover"
tdf_dir <- "data/tradeoff"
c1_dir <- "/Users/pinot/dropbox/research/hrlcm"
c3_dir <- "/Users/pinot/dropbox/research/eleDistribution"
c4_dir <- "/Users/pinot/dropbox/research/comTCA"

# Study area
bry <- read_sf("data/geoms/mainland_tanzania.geojson") %>% 
    select() %>% mutate(name = "Tanzania")

ecozones <- read_sf(file.path(c1_dir, "data/ecozones/ecozones_shp.shp")) %>% 
    st_intersection(bry)

# Mountains
mnts <- read_sf(file.path(c3_dir, "data/mountains/GMBA_Inventory_v2",
                          "GMBA_Inventory_v2.0_standard.shp")) %>% 
    filter(MapUnit == "Basic") %>% 
    filter(str_detect(Countries, "Tanzania")) %>% 
    st_intersection(bry)
mnts <- mnts %>% mutate(height = Elev_High - Elev_Low) %>% 
    filter(height > 1000)

# Protected areas
pas <- read_sf(
    file.path("data/protected_area", "WDPA_WDOECM_Jan2023_Public_TZA",
              "WDPA_WDOECM_Jan2023_Public_TZA.gdb"),
    layer = "WDPA_WDOECM_poly_Jan2023_TZA")
pas <- pas %>% 
    select(WDPAID, NAME, DESIG) %>% 
    rename(Geometry = SHAPE) %>% 
    filter(!DESIG %in% c("Marine Reserve", "Marine Park"))
pas <- pas %>% 
    slice(unique(unlist(st_intersects(bry, pas)))) %>% 
    st_simplify()

cropland <- rast(file.path(c4_dir, "data/tradeoff/cropland_perc.tif"))
cropland[cropland < 0.1] <- NA

# Make the figure
## Overview
countries <- ne_countries(
    type = 'countries', scale = 'large',
    continent = "Africa", returnclass = 'sf')
east_africa <- read_sf(file.path("data/geoms/east_africa.sqlite"))

# Make overview
overview <- ggplotGrob(ggplot() +
                           with_shadow(
                               geom_sf(data = countries, fill = 'lightgrey', 
                                       color = 'black', linewidth = 0.3),
                               sigma = 1, x_offset = 10, y_offset = 10) +
                           geom_sf(data = bry, fill = 'coral', 
                                   color = 'black', linewidth = 0.3) +
                           coord_sf(xlim = c(-17.892, 51.892), 
                                    ylim = c(-35.321, 37.715)) +
                           theme_transparent() +
                           theme(axis.title = element_blank(),
                                 axis.text = element_blank(),
                                 axis.ticks = element_blank(),
                                 panel.grid.major = element_blank(),
                                 plot.margin = unit(rep(-0.6, 4), "cm")))

## Main figure
g <- ggplot() +
    geom_sf(data = ecozones, aes(fill = PHY_REG, color = PHY_REG),
            linewidth = 0.2) +
    scale_fill_hypso_d(
        palette = "wiki-2.0_hypso",
        guide = guide_legend(
            title.position = "top", 
            title.hjust = 0.5, direction = "horizontal", ncol = 1)) +
    scale_color_hypso_d(
        palette = "wiki-2.0_hypso",
        guide = guide_legend(
            title.position = "top", 
            title.hjust = 0.5, direction = "horizontal", ncol = 1)) +
    labs(fill = "Agro-ecological zone", color = "Agro-ecological zone") +
    geom_sf(data = bry, color = "black", 
            fill = "transparent", linewidth = 0.8) +
    annotation_custom(grob = overview, 
                      xmin = 37.2, xmax = 40.1,
                      ymin = - 3.8, ymax = -1.4) +
    coord_sf() + theme_void() +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 10, face = "bold"),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(rep(0, 4)), 
          legend.key.height = unit(0.4, 'cm'),
          legend.key.width = unit(0.4, 'cm'),
          plot.margin = unit(c(0, 0, 0, 0), "cm"),
          legend.position = "bottom",
          legend.box.just = "left",
          legend.box = "vertical")

g
ggsave("figures/aag_fig1-1.png",
       bg = "transparent",
       width = 4, height = 6, dpi = 500)

g <- g + 
    geom_sf(data = pas, color = "white", 
            fill = "#00968f", linewidth = 0.4)
g
ggsave("figures/aag_fig1-2.png",
       bg = "transparent",
       width = 4, height = 6, dpi = 500)

g <- g +
    new_scale_fill() +
    geom_spatraster(data = cropland, show.legend = FALSE) +
    scale_fill_material("amber", na.value = NA)

g
ggsave("figures/aag_fig1-3.png",
       bg = "transparent",
       width = 4, height = 6, dpi = 500)

# TOP 5 country in Africa
## FAO
pops <- data.frame(
    year = seq(2025, 2050, 5),
    tanzania = c(71.4, 81.9, 93.1, 105.0, 117.3, 130.0), # million
    world = c(8.2, 8.5, 8.8, 9.2, 9.4, 9.7)) # billion

pops <- pops %>% 
    mutate(label = paste0(round(tanzania, 0), "M"))

m <- png::readPNG("figures/community.png")
g <- matrix(rgb(m[,,1],m[,,2],m[,,3], m[,,4] * 0.5), nrow=dim(m)[1])
g <- grid::rasterGrob(g, interpolate = TRUE)

ani <- ggplot(pops, aes(x = year, y = tanzania, label = tanzania)) + 
    annotation_custom(g, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_line(linewidth = 2) +
    geom_point(shape = 21, fill = "white", size = 4, stroke = 2,
               aes(group = seq_along(year))) +
    geom_text(aes(group = seq_along(year)), hjust = 0.5, vjust = -1.1,
              size = 7, fontface = 'bold') +
    geom_point(shape = 21, fill = "coral", size = 4,
               stroke = 2.5, color = "coral") +
    ylim(70, 140) +
    labs(x = "Year", y = "Population in millions") +
    theme_fivethirtyeight() +
    theme(text = element_text(size = 18, face = "bold"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    transition_reveal(year)

anim_save("figures/aag_fig2.gif", animation = ani, duration = 3,
          renderer = gifski_renderer(loop = FALSE),
          width = 500, height = 300)

# Simplify figure 2 in the manuscript
f_dir <- "results/scenarios/summary"
items <- c("Y", "I", "B", "Ca", "Co", "D", "YBCaCoD")

# Load dataset
vals_cmpr <- read.csv(file.path(f_dir, "summary_scenarios_USGs.csv")) %>% 
    filter(land_usage == 0.64) %>% 
    filter(change == "Y100") %>% 
    filter(scenario %in% items) %>% 
    mutate(scenario = factor(
        scenario, levels = items, labels = items)) %>% 
    mutate(farmable_area_intense = farmable_area / num_unit) %>% 
    select(-c(land_usage, change, num_unit, weight))
names(vals_cmpr) <- c("B", "Ca", "Co", "D", "Y", "scenario", "I")

# Figure
fig_list <- lapply(c("Y", "I", "B", "Ca", "Co", "D"), function(item){
    # Clean data
    # item_order <- c(item, setdiff(items, item))
    item_order <- items
    vals <- vals_cmpr %>% 
        select(all_of(c(item, "scenario"))) %>% 
        mutate(scenario  = factor(
            scenario, levels = item_order, labels = item_order)) %>% 
        rename(value = all_of(item))
    item_mod <- ifelse(item == 'I', "Y", item)
    vals <- vals %>% 
        mutate(group = ifelse(
            scenario == item_mod, "Optimal", 
            ifelse(scenario %in% setdiff(c("B", "Ca", "Co", "D", "Y"), item_mod), 
                   "Single factor", "Combined factors"))) %>% 
        mutate(group = factor(
            group, levels = c("Single factor", "Optimal", "Combined factors"),
            labels = c("Single factor (", "Optimal for each factor)", "Hybrid")))
    
    # Separate different conditions
    if (str_detect(item, "Y")){
        vals$value <- vals$value / 1e6
        xlab_nm <- expression(bold(paste("New cropland area (ha * ", bold("10")^bold("6"), ")      ")))
    } else if (str_detect(item, "I")){
        vals$value <- vals$value
        xlab_nm <- "New cropland area / No. of unit (ha)                  "
    } else if (str_detect(item, "B")) {
        vals$value <- vals$value / 1e6
        xlab_nm <- expression(bold(paste("Biodiversity loss (", bold("10")^bold("6"), ")")))
    } else if (str_detect(item, "Ca")){
        vals$value <- vals$value / 1e8
        xlab_nm <- expression(bold(paste("Carbon loss (tons * ", bold("10")^bold("8"), ")     ")))
    } else if (str_detect(item, "Co")){
        vals$value <- vals$value / 1e6
        xlab_nm <- expression(bold(paste("Connectivity loss (", bold("10")^bold("5"), ")")))
    } else if (str_detect(item, "D")){
        vals$value <- vals$value / 60 / 1e6
        xlab_nm <- expression(bold(paste("Travel time (hours * ", bold("10")^bold("6"), ")            ")))}
    
    # Make figure
    ggplot(data = vals) + 
        geom_segment(aes(
            y = scenario, yend = scenario, 
            x = 0, xend = value, color = group), linewidth = 0.8) +
        geom_point(aes(y = scenario, x = value, color = group), size = 2) +
        scale_color_manual(
            name = "Solution Type", 
            values = c("grey", pal_npg("nrc")(1), 
                       pal_npg("nrc")(3)[3])) + 
        scale_y_discrete(labels=c("B" = "B", "Ca" = "Ca", "Co" = "Co", "D" = "D", "Y" = "Y",
            "YBCaCoD" = expression(bold(YBCaCoD)), parse=TRUE)) +
        theme_light() +
        theme(panel.grid.major.y = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank()) +
        ylab("") +
        xlab(xlab_nm) + theme_light() +
        theme(panel.grid.major.y = element_blank(),
              panel.grid.major.x = element_line(
                  color = "lightgrey", linetype = "dashed"),
              panel.grid.minor.x = element_blank(),
              panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_text(face = "bold"),
              text = element_text(size = 10, color = "black"),
              axis.title = element_text(color = "black", size = 10),
              axis.text = element_text(color = "black", size = 10),
              legend.text = element_text(size = 12, color = "black"),
              legend.title = element_text(size = 12, color = "black"),
              legend.position = "top",
              plot.margin = unit(c(0.1, 0.2, 0.1, 0.2), "cm"))
})

ggarrange(plotlist = fig_list, nrow = 2, ncol = 3, vjust = 5,
          common.legend = TRUE, legend = "top")

ggsave("figures/aag_fig3.png",
       bg = "white",
       width = 7.5, height = 3.5, dpi = 500)

# Summarize the maps
f_dir <- "results/scenarios/spatial"
fnames <- list.files(f_dir, full.names = TRUE, pattern = "Y100_USG64.tif$")
fnames <- fnames[-c(2, 8, 9)]
maps <- lapply(fnames, function(fname){
    rst <- rast(fname)
    rst[is.na(rst)] <- 0
    lyrs <- do.call(c, lapply(0:5, function(x) {
        if (x == 0){
            rst > 0
        } else rst == x
    }))
    names(lyrs) <- c("farmland", "Maize", "Paddy", "Sorghum", "Cassava", "Beans")
    lyrs
})

layers <- do.call(c, lapply(names(maps[[1]]), function(lyr){
    do.call(c, lapply(maps, function(x) x[lyr])) %>% sum()
})) %>% crop(bry) %>% mask(bry)
names(layers) <- names(maps[[1]])
layers[layers == 0] <- NA

# Make the figure
# Protected areas
pas <- read_sf(
    file.path("data/protected_area", "WDPA_WDOECM_Jan2023_Public_TZA",
              "WDPA_WDOECM_Jan2023_Public_TZA.gdb"),
    layer = "WDPA_WDOECM_poly_Jan2023_TZA")
pas <- pas %>% 
    select(WDPAID, NAME, DESIG) %>% 
    rename(Geometry = SHAPE) %>% 
    filter(!DESIG %in% c("Marine Reserve", "Marine Park"))
pas <- pas %>% 
    slice(unique(unlist(st_intersects(bry, pas)))) %>% 
    st_simplify()

# Highlight areas
areas <- st_read('data/osm/areas_fig3.geojson') %>% 
    st_intersection(., bry)
areas$name[2] <- "Lake\nVictoria"

# Overall agreement
ggplot() +
    geom_sf(data = pas, color = "#E5F5E0", 
            fill = "#E5F5E0", linewidth = 0) +
    geom_spatraster(data = layers$farmland, na.rm = TRUE) +
    scale_fill_continuous_sequential(
        name = "Agreement level (1 - 6)", breaks = c(1, 3, 5, 7, 9),
        palette = "Plasma", na.value = NA,
        guide = guide_colourbar(
            title.position = "top", title.hjust = 0.5)) +
    geom_sf(data = areas, color = "black",
            fill = "transparent", linewidth = 0.3) +
    geom_sf_text(data = areas[2, ], fontface='bold',
                 aes(label = "Lake"), nudge_y = 0.8, nudge_x = 0.5,
                 colour = "black", size = 3) +
    geom_sf_text(data = areas[2, ], fontface='bold',
                 aes(label = "Victoria"), nudge_y = 0.5, nudge_x = 0.5,
                 colour = "black", size = 3) +
    geom_sf_text_repel(data = areas[3, ], fontface='bold',
                       aes(label = name), colour = "black", size = 3, 
                       nudge_y = 1.2, nudge_x = 0.3, force = 100) +
    geom_sf_text_repel(data = areas[1, ], fontface='bold',
                       aes(label = name), colour = "black", size = 3, 
                       nudge_y = 0.9, nudge_x = -0.2, force = 100) +
    geom_sf_text_repel(data = areas[4, ], fontface='bold',
                       aes(label = name), colour = "black", size = 3, 
                       nudge_y = 1.2, nudge_x = -1.4, force = 100) +
    geom_sf_text_repel(data = areas[5, ], fontface='bold',
                       aes(label = name), colour = "black", size = 3, 
                       nudge_y = -1.5, nudge_x = 0, force = 100) +
    geom_sf(data = bry, color = "black", 
            fill = "transparent", linewidth = 0.6) +
    coord_sf() + theme_void() +
    theme(text = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(0.1, 'cm'),
          legend.margin = margin(rep(0, 4)), 
          legend.key.height = unit(0.4, 'cm'),
          legend.key.width = unit(0.8, 'cm'),
          plot.margin = unit(c(0, -1, 0, -1), "cm"),
          legend.position = "bottom",
          legend.box = "vertical")

ggsave("figures/aag_fig3_2.png",
       bg = "white",
       width = 4, height = 4.2, dpi = 500)

# Cropland explode figure
cropland <- rast(file.path(c4_dir, "data/tradeoff/cropland_perc.tif"))
msk <- cropland
cropland[cropland < 0.1] <- 0
cropland <- mask(cropland, pas, inverse = TRUE)

# Get values
crps <- cropland
for (n in 2:10){
    set.seed(n)
    vals <- values(cropland) %>% data.frame()
    ids_to_select <- which(!is.na(vals$cropland_ratio) & vals$cropland_ratio < 0.9)
    ids_to_select <- sample(ids_to_select, 64483)
    
    vals$cropland_ratio[ids_to_select] <- 1
    if (n == 10){
        vals$cropland_ratio[!is.na(vals$cropland_ratio)] <- 1
        vals$cropland_ratio[1] <- 0
    }
    new_cropland <- cropland
    values(new_cropland) <- vals
    
    cropland <- new_cropland
    crps <- c(crps, cropland)
}

crps_to_pa <- crps[[10]]
cropland <- crps_to_pa

for (n in 2:7){
    new_cropland <- buffer(cropland, 20000)
    new_cropland <- cover(cropland, new_cropland)
    new_cropland[new_cropland == 0] <- NA
    new_cropland <- mask(new_cropland, msk)
    
    cropland <- new_cropland
    crps_to_pa <- c(crps_to_pa, cropland)
}

crps <- c(crps, crps_to_pa)
names(crps) <- 1:17

tmap_options(output.dpi.animation = 300)
anim <- tm_shape(pas) + tm_polygons(col = "#00968f", border.col = "#00968f") +
    tm_shape(crps, raster.downsample = FALSE,
                 title = "Cropland coverage") + 
    tm_raster(palette = pal_material("amber")(7)) + 
    tm_shape(bry) + tm_borders(lwd = 2) +
    tm_facets(nrow = 1, ncol = 1) +
    tm_layout(frame = FALSE, panel.show = FALSE, legend.show = FALSE)

tmap_animation(anim, "figures/aag_fig3.gif", loop = TRUE, fps = 2)


# Summarize the maps
fnames <- list.files(
    "results/scenarios/spatial", 
    full.names = TRUE, pattern = "^YBCaCoD.+.tif$")
maps <- lapply(fnames, function(fname){
    rst <- rast(fname)
    rst[is.na(rst)] <- 0
    lyrs <- do.call(c, lapply(0:5, function(x) {
        if (x == 0){
            rst > 0
        } else rst == x
    }))
    names(lyrs) <- c("farmland", "Maize", "Paddy", "Sorghum", "Cassava", "Beans")
    lyrs
})

# Make the figure
# Protected areas
pas <- read_sf(
    file.path("data/protected_area", "WDPA_WDOECM_Jan2023_Public_TZA",
              "WDPA_WDOECM_Jan2023_Public_TZA.gdb"),
    layer = "WDPA_WDOECM_poly_Jan2023_TZA")
pas <- pas %>% 
    select(WDPAID, NAME, DESIG) %>% 
    rename(Geometry = SHAPE) %>% 
    filter(!DESIG %in% c("Marine Reserve", "Marine Park"))
pas <- pas %>% 
    slice(unique(unlist(st_intersects(bry, pas)))) %>% 
    st_simplify()

# Spatial land allocation for each crop
practices <- str_extract(fnames, "Y100|CASS3|Y140")
usages <- str_extract(fnames, "USG64|USG80")
crop_cols <- c(
    "#BF812D", "#543005",  "#e9a3c9", "#c51b7d", "#fdae6b", 
    "#e6550d", "#9ecae1", "#3182bd",  "#bcbddc", "#756bb1")

# Make legend
dat <- data.frame(
    rst = c("Maize", "Maize2", "Paddy rice", "Paddy rice2", 
            "Sorghum", "Sorghum2", "Cassava", "Cassava2", 
            "Common beans", "Common beans2"), 
    x = 1:10, y = 1:10) %>% 
    mutate(rst = factor(
        rst, levels = c("Maize", "Maize2", "Paddy rice", "Paddy rice2", 
                        "Sorghum", "Sorghum2", "Cassava", "Cassava2", 
                        "Common beans", "Common beans2"),
        labels = c("Maize", "Maize2", "Paddy rice", "Paddy rice2", 
                   "Sorghum", "Sorghum2", "Cassava", "Cassava2", 
                   "Common beans", "Common beans2")))
dat2 <- data.frame(x = 1:2, y = 1:2, z = c("64", "80")) %>% 
    mutate(z = factor(z, levels = c("64", "80"), labels = c("64", "80")))
legend <-  ggplot() + 
    geom_raster(data = dat2, aes(x = x, y = y, fill = z)) +
    scale_fill_manual('Cropland usage rate',
                      values = c("white", "white"), 
                      labels = c("                                      64%", 
                                 "                                      80%"),
                      guide = guide_legend(
                          title.position = "top", 
                          title.hjust = 0.4, title.vjust = -11, 
                          direction = "vertical", 
                          order = 1)) +
    new_scale_fill() +
    geom_raster(data = dat,  aes(x = x, y = y, fill = rst)) +
    scale_fill_manual('Suggested crop',
                      values = crop_cols,
                      labels = c("Maize", "", "Paddy rice", "", "Sorghum", "", 
                                 "Cassava", "", "Common beans", ""),
                      guide = guide_legend(
                          title.position = "top", 
                          title.hjust = 0.5, direction = "horizontal",
                          order = 2)) +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position = "top", 
          legend.box = "horizontal")
legend_lc <- get_legend(legend)

fig_list <- lapply(c("Y100", "CASS3", "Y140"), function(practice){
    if (practice %in% c("CASS3", "Y140")){
        ggplot() +
            geom_sf(data = bry, color = "black", 
                    fill = "transparent", linewidth = 0) +
            coord_sf() + theme_void() +
            theme(text = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.spacing.y = unit(0.1, 'cm'),
                  legend.margin = margin(rep(0, 4)), 
                  legend.key.height = unit(0.4, 'cm'),
                  legend.key.width = unit(0.8, 'cm'),
                  plot.margin = unit(c(0, -1, 0, -1), "cm"),
                  legend.position = "none")
    } else {
        map <- maps[[c(1:6)[practices == practice & usages == "USG64"]]]
        map <- map %>% crop(bry) %>% mask(bry)
        map[map == 0] <- NA
        
        # Put them together into the same layer
        map <- sum(map[[2]], (map[[3]] + 1), (map[[4]] + 2), (map[[5]] + 3), 
                   (map[[6]] + 4), na.rm = TRUE)
        
        # Set to categorical
        map <- as.factor(map)
        levels(map) <- data.frame(
            value = 1:5, 
            usage = c("Maize", "Paddy rice", "Sorghum", 
                      "Cassava", "Common beans"))
        
        # Get sub_map
        sub_ply <- ext(c(xmin = 35.1, xmax = 36.1, 
                         ymin = -4.2, ymax = -3.2)) %>% 
            vect(crs = "EPSG:4326") %>% st_as_sf(crs = 4326) %>% 
            st_centroid() %>% st_buffer(0.5)
        sub_map <- crop(map, sub_ply) %>% mask(sub_ply)
        cols_sub <- crop_cols[c(1,3,5,7,9)][na.omit(freq(map)[["value"]]) %in% freq(sub_map)[["value"]]]
        zoomin <- ggplotGrob(
            ggplot() +
                geom_spatraster(data = sub_map, na.rm = TRUE) +
                scale_fill_manual(
                    values = cols_sub,
                    na.translate = FALSE) +
                geom_sf(data = sub_ply, color = "black", 
                        fill = "transparent", linewidth = 0.6) +
                coord_sf() + theme_void() +
                theme(text = element_text(size = 14),
                      legend.text = element_text(size = 14),
                      legend.spacing.y = unit(0.1, 'cm'),
                      legend.margin = margin(rep(0, 4)),
                      legend.key.height = unit(0.4, 'cm'),
                      legend.key.width = unit(0.8, 'cm'),
                      plot.margin = unit(c(0, -1, 0, -1), "cm"),
                      legend.position = "none"))
        
        # Figure
        ggplot() +
            geom_spatraster(data = map, na.rm = TRUE) +
            scale_fill_manual(
                values = crop_cols[c(1,3,5,7,9)],
                na.translate = FALSE) +
            geom_sf(data = pas, color = "#f0f0f0", 
                    fill = "#f0f0f0", linewidth = 0) +
            geom_sf(data = bry, color = "black", 
                    fill = "transparent", linewidth = 0.6) +
            geom_sf(data = sub_ply, color = "black", 
                    fill = "transparent", linewidth = 0.4) +
            annotation_custom(grob = zoomin,
                              xmin = 29.9, xmax = 33.1,
                              ymin = - 12.1, ymax = - 8.9) +
            coord_sf() + theme_void() +
            theme(text = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.spacing.y = unit(0.1, 'cm'),
                  legend.margin = margin(rep(0, 4)), 
                  legend.key.height = unit(0.4, 'cm'),
                  legend.key.width = unit(0.8, 'cm'),
                  plot.margin = unit(c(0, -1, 0, -1), "cm"),
                  legend.position = "none")
    }
})

g1 <- ggarrange(legend_lc, ggarrange(plotlist = fig_list, nrow = 1, ncol = 3,
                                     labels = c("BL", "Sce1", "Sce2"),
                                     font.label = list(size = 14),
                                     vjust = 3, hjust = 0), 
                nrow = 2, ncol = 1, heights = c(3, 10))

# New area and ecological costs change
vals <- read.csv(
    file.path("results/scenarios/summary", 
              "summary_scenarios_USGs.csv")) %>% 
    filter(scenario == "YBCaCoD") %>% 
    select(-c(scenario, num_unit, weight))
names(vals) <- c("B", "Ca", "Co", "D", "Y", "USG", "Practice")

vals <- vals %>% pivot_longer(
    c("B", "Ca", "Co", "D", "Y"), names_to = "Element")
vals <- vals %>% group_by(Element) %>% 
    reframe(USG = USG, Practice = Practice,
            value = value / max(value) * 100) %>% 
    mutate(USG = factor(USG, labels = c("64%", "80%"), levels = c(0.64, 0.8))) %>% 
    mutate(Practice = factor(Practice, labels = c("BL", "Sce1", "Sce2"),
                             levels = c("Y100", "CASS3", "Y140"))) %>% 
    mutate(Element = factor(
        Element, levels = c("Y", "B", "Ca", "Co", "D"),
        labels = c("New cropland area", "Biodiveristy loss", "Carbon loss",
                   "Connectivity loss", "Travel time")))
vals <- vals %>% mutate(group = ifelse(USG == "64%" & Practice == "BL", "B", "P"))
vals <- vals %>% mutate(value = ifelse(Practice %in% c("Sce1", "Sce2"), NA, value))
vals <- vals %>% mutate(value = ifelse(USG == "80%", NA, value))

g2 <- ggplot(data = vals) + 
    geom_linerange(aes(
        x = Practice, ymin = 0, ymax = value, linetype = USG, color = group), 
        linewidth = 0.8, position = position_dodge2(width = 0.7)) +
    geom_point(aes(x = Practice, y = value, color = group), size = 2,
               position = position_dodge2(width = 0.7)) +
    geom_text(data = vals %>% mutate(value = ifelse(value == 100, NA, value)), 
              aes(x = Practice, y = value + 8, label = round(value, 1)), 
              size = 3, color = "black", fontface = "bold",
              position = position_dodge2(width = 0.7)) + 
    scale_color_manual(values = c("#d73027", "black")) +
    xlab("") + ylab("Percentage") + labs(linetype = "Cropland usage rate") +
    facet_grid(~Element) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(
              color = "lightgrey", linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank(),
          text = element_text(size = 14, color = "black"),
          axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(color = "black", size = 14),
          axis.text.x = element_text(face = "bold"),
          legend.text = element_text(size = 14, color = "black"),
          legend.title = element_text(size = 14, color = "black"),
          legend.position = "top",
          legend.margin = margin(rep(0, 4)),
          plot.margin = unit(rep(0, 4), "cm")) +
    guides(color = "none")

ggarrange(g1, g2, nrow = 2, ncol = 1, heights = c(6, 5))

ggsave("figures/aag_fig4_1.png",
       bg = "white",
       width = 9.6, height = 6.4, dpi = 500)

fig_list <- lapply(c("Y100", "CASS3", "Y140"), function(practice){
    if (practice %in% c("Y140")){
        ggplot() +
            geom_sf(data = bry, color = "black", 
                    fill = "transparent", linewidth = 0) +
            coord_sf() + theme_void() +
            theme(text = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.spacing.y = unit(0.1, 'cm'),
                  legend.margin = margin(rep(0, 4)), 
                  legend.key.height = unit(0.4, 'cm'),
                  legend.key.width = unit(0.8, 'cm'),
                  plot.margin = unit(c(0, -1, 0, -1), "cm"),
                  legend.position = "none")
    } else {
        map <- maps[[c(1:6)[practices == practice & usages == "USG64"]]]
        map <- map %>% crop(bry) %>% mask(bry)
        map[map == 0] <- NA
        
        # Put them together into the same layer
        map <- sum(map[[2]], (map[[3]] + 1), (map[[4]] + 2), (map[[5]] + 3), 
                   (map[[6]] + 4), na.rm = TRUE)
        
        # Set to categorical
        map <- as.factor(map)
        levels(map) <- data.frame(
            value = 1:5, 
            usage = c("Maize", "Paddy rice", "Sorghum", 
                      "Cassava", "Common beans"))
        
        # Get sub_map
        sub_ply <- ext(c(xmin = 35.1, xmax = 36.1, 
                         ymin = -4.2, ymax = -3.2)) %>% 
            vect(crs = "EPSG:4326") %>% st_as_sf(crs = 4326) %>% 
            st_centroid() %>% st_buffer(0.5)
        sub_map <- crop(map, sub_ply) %>% mask(sub_ply)
        cols_sub <- crop_cols[c(1,3,5,7,9)][na.omit(freq(map)[["value"]]) %in% freq(sub_map)[["value"]]]
        zoomin <- ggplotGrob(
            ggplot() +
                geom_spatraster(data = sub_map, na.rm = TRUE) +
                scale_fill_manual(
                    values = cols_sub,
                    na.translate = FALSE) +
                geom_sf(data = sub_ply, color = "black", 
                        fill = "transparent", linewidth = 0.6) +
                coord_sf() + theme_void() +
                theme(text = element_text(size = 14),
                      legend.text = element_text(size = 14),
                      legend.spacing.y = unit(0.1, 'cm'),
                      legend.margin = margin(rep(0, 4)),
                      legend.key.height = unit(0.4, 'cm'),
                      legend.key.width = unit(0.8, 'cm'),
                      plot.margin = unit(c(0, -1, 0, -1), "cm"),
                      legend.position = "none"))
        
        # Figure
        ggplot() +
            geom_spatraster(data = map, na.rm = TRUE) +
            scale_fill_manual(
                values = crop_cols[c(1,3,5,7,9)],
                na.translate = FALSE) +
            geom_sf(data = pas, color = "#f0f0f0", 
                    fill = "#f0f0f0", linewidth = 0) +
            geom_sf(data = bry, color = "black", 
                    fill = "transparent", linewidth = 0.6) +
            geom_sf(data = sub_ply, color = "black", 
                    fill = "transparent", linewidth = 0.4) +
            annotation_custom(grob = zoomin,
                              xmin = 29.9, xmax = 33.1,
                              ymin = - 12.1, ymax = - 8.9) +
            coord_sf() + theme_void() +
            theme(text = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.spacing.y = unit(0.1, 'cm'),
                  legend.margin = margin(rep(0, 4)), 
                  legend.key.height = unit(0.4, 'cm'),
                  legend.key.width = unit(0.8, 'cm'),
                  plot.margin = unit(c(0, -1, 0, -1), "cm"),
                  legend.position = "none")
    }
})

g1 <- ggarrange(legend_lc, ggarrange(plotlist = fig_list, nrow = 1, ncol = 3,
                                     labels = c("BL", "Sce1", "Sce2"),
                                     font.label = list(size = 14),
                                     vjust = 3, hjust = 0), 
                nrow = 2, ncol = 1, heights = c(3, 10))

# New area and ecological costs change
vals <- read.csv(
    file.path("results/scenarios/summary", 
              "summary_scenarios_USGs.csv")) %>% 
    filter(scenario == "YBCaCoD") %>% 
    select(-c(scenario, num_unit, weight))
names(vals) <- c("B", "Ca", "Co", "D", "Y", "USG", "Practice")

vals <- vals %>% pivot_longer(
    c("B", "Ca", "Co", "D", "Y"), names_to = "Element")
vals <- vals %>% group_by(Element) %>% 
    reframe(USG = USG, Practice = Practice,
            value = value / max(value) * 100) %>% 
    mutate(USG = factor(USG, labels = c("64%", "80%"), levels = c(0.64, 0.8))) %>% 
    mutate(Practice = factor(Practice, labels = c("BL", "Sce1", "Sce2"),
                             levels = c("Y100", "CASS3", "Y140"))) %>% 
    mutate(Element = factor(
        Element, levels = c("Y", "B", "Ca", "Co", "D"),
        labels = c("New cropland area", "Biodiveristy loss", "Carbon loss",
                   "Connectivity loss", "Travel time")))
vals <- vals %>% mutate(group = ifelse(USG == "64%" & Practice == "BL", "B", "P"))
vals <- vals %>% mutate(value = ifelse(Practice %in% c("Sce2"), NA, value))
vals <- vals %>% mutate(value = ifelse(USG == "80%", NA, value))

g2 <- ggplot(data = vals) + 
    geom_linerange(aes(
        x = Practice, ymin = 0, ymax = value, linetype = USG, color = group), 
        linewidth = 0.8, position = position_dodge2(width = 0.7)) +
    geom_point(aes(x = Practice, y = value, color = group), size = 2,
               position = position_dodge2(width = 0.7)) +
    geom_text(data = vals %>% mutate(value = ifelse(value == 100, NA, value)), 
              aes(x = Practice, y = value + 8, label = round(value, 1)), 
              size = 3, color = "black", fontface = "bold",
              position = position_dodge2(width = 0.7)) + 
    scale_color_manual(values = c("#d73027", "black")) +
    xlab("") + ylab("Percentage") + labs(linetype = "Cropland usage rate") +
    facet_grid(~Element) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(
              color = "lightgrey", linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank(),
          text = element_text(size = 14, color = "black"),
          axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(color = "black", size = 14),
          axis.text.x = element_text(face = "bold"),
          legend.text = element_text(size = 14, color = "black"),
          legend.title = element_text(size = 14, color = "black"),
          legend.position = "top",
          legend.margin = margin(rep(0, 4)),
          plot.margin = unit(rep(0, 4), "cm")) +
    guides(color = "none")

ggarrange(g1, g2, nrow = 2, ncol = 1, heights = c(6, 5))

ggsave("figures/aag_fig4_2.png",
       bg = "white",
       width = 9.6, height = 6.4, dpi = 500)

fig_list <- lapply(c("Y100", "CASS3", "Y140"), function(practice){
    map <- maps[[c(1:6)[practices == practice & usages == "USG64"]]]
    map <- map %>% crop(bry) %>% mask(bry)
    map[map == 0] <- NA
    
    # Put them together into the same layer
    map <- sum(map[[2]], (map[[3]] + 1), (map[[4]] + 2), (map[[5]] + 3), 
               (map[[6]] + 4), na.rm = TRUE)
    
    # Set to categorical
    map <- as.factor(map)
    levels(map) <- data.frame(
        value = 1:5, 
        usage = c("Maize", "Paddy rice", "Sorghum", 
                  "Cassava", "Common beans"))
    
    # Get sub_map
    sub_ply <- ext(c(xmin = 35.1, xmax = 36.1, 
                     ymin = -4.2, ymax = -3.2)) %>% 
        vect(crs = "EPSG:4326") %>% st_as_sf(crs = 4326) %>% 
        st_centroid() %>% st_buffer(0.5)
    sub_map <- crop(map, sub_ply) %>% mask(sub_ply)
    cols_sub <- crop_cols[c(1,3,5,7,9)][na.omit(freq(map)[["value"]]) %in% freq(sub_map)[["value"]]]
    zoomin <- ggplotGrob(
        ggplot() +
            geom_spatraster(data = sub_map, na.rm = TRUE) +
            scale_fill_manual(
                values = cols_sub,
                na.translate = FALSE) +
            geom_sf(data = sub_ply, color = "black", 
                    fill = "transparent", linewidth = 0.6) +
            coord_sf() + theme_void() +
            theme(text = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.spacing.y = unit(0.1, 'cm'),
                  legend.margin = margin(rep(0, 4)),
                  legend.key.height = unit(0.4, 'cm'),
                  legend.key.width = unit(0.8, 'cm'),
                  plot.margin = unit(c(0, -1, 0, -1), "cm"),
                  legend.position = "none"))
    
    # Figure
    ggplot() +
        geom_spatraster(data = map, na.rm = TRUE) +
        scale_fill_manual(
            values = crop_cols[c(1,3,5,7,9)],
            na.translate = FALSE) +
        geom_sf(data = pas, color = "#f0f0f0", 
                fill = "#f0f0f0", linewidth = 0) +
        geom_sf(data = bry, color = "black", 
                fill = "transparent", linewidth = 0.6) +
        geom_sf(data = sub_ply, color = "black", 
                fill = "transparent", linewidth = 0.4) +
        annotation_custom(grob = zoomin,
                          xmin = 29.9, xmax = 33.1,
                          ymin = - 12.1, ymax = - 8.9) +
        coord_sf() + theme_void() +
        theme(text = element_text(size = 14),
              legend.text = element_text(size = 14),
              legend.spacing.y = unit(0.1, 'cm'),
              legend.margin = margin(rep(0, 4)), 
              legend.key.height = unit(0.4, 'cm'),
              legend.key.width = unit(0.8, 'cm'),
              plot.margin = unit(c(0, -1, 0, -1), "cm"),
              legend.position = "none")
})

g1 <- ggarrange(legend_lc, ggarrange(plotlist = fig_list, nrow = 1, ncol = 3,
                                     labels = c("BL", "Sce1", "Sce2"),
                                     font.label = list(size = 14),
                                     vjust = 3, hjust = 0), 
                nrow = 2, ncol = 1, heights = c(3, 10))

# New area and ecological costs change
vals <- read.csv(
    file.path("results/scenarios/summary", 
              "summary_scenarios_USGs.csv")) %>% 
    filter(scenario == "YBCaCoD") %>% 
    select(-c(scenario, num_unit, weight))
names(vals) <- c("B", "Ca", "Co", "D", "Y", "USG", "Practice")

vals <- vals %>% pivot_longer(
    c("B", "Ca", "Co", "D", "Y"), names_to = "Element")
vals <- vals %>% group_by(Element) %>% 
    reframe(USG = USG, Practice = Practice,
            value = value / max(value) * 100) %>% 
    mutate(USG = factor(USG, labels = c("64%", "80%"), levels = c(0.64, 0.8))) %>% 
    mutate(Practice = factor(Practice, labels = c("BL", "Sce1", "Sce2"),
                             levels = c("Y100", "CASS3", "Y140"))) %>% 
    mutate(Element = factor(
        Element, levels = c("Y", "B", "Ca", "Co", "D"),
        labels = c("New cropland area", "Biodiveristy loss", "Carbon loss",
                   "Connectivity loss", "Travel time")))
vals <- vals %>% mutate(group = ifelse(USG == "64%" & Practice == "BL", "B", "P"))
vals <- vals %>% mutate(value = ifelse(USG == "80%", NA, value))

g2 <- ggplot(data = vals) + 
    geom_linerange(aes(
        x = Practice, ymin = 0, ymax = value, linetype = USG, color = group), 
        linewidth = 0.8, position = position_dodge2(width = 0.7)) +
    geom_point(aes(x = Practice, y = value, color = group), size = 2,
               position = position_dodge2(width = 0.7)) +
    geom_text(data = vals %>% mutate(value = ifelse(value == 100, NA, value)), 
              aes(x = Practice, y = value + 8, label = round(value, 1)), 
              size = 3, color = "black", fontface = "bold",
              position = position_dodge2(width = 0.7)) + 
    scale_color_manual(values = c("#d73027", "black")) +
    xlab("") + ylab("Percentage") + labs(linetype = "Cropland usage rate") +
    facet_grid(~Element) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(
              color = "lightgrey", linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank(),
          text = element_text(size = 14, color = "black"),
          axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(color = "black", size = 14),
          axis.text.x = element_text(face = "bold"),
          legend.text = element_text(size = 14, color = "black"),
          legend.title = element_text(size = 14, color = "black"),
          legend.position = "top",
          legend.margin = margin(rep(0, 4)),
          plot.margin = unit(rep(0, 4), "cm")) +
    guides(color = "none")

ggarrange(g1, g2, nrow = 2, ncol = 1, heights = c(6, 5))

ggsave("figures/aag_fig4_3.png",
       bg = "white",
       width = 9.6, height = 6.4, dpi = 500)

fig_list <- lapply(c("Y100", "CASS3", "Y140"), function(practice){
    map_64 <- maps[[c(1:6)[practices == practice & usages == "USG64"]]]
    map_80 <- maps[[c(1:6)[practices == practice & usages == "USG80"]]]
    map <- (map_64 + map_80) %>% crop(bry) %>% mask(bry)
    map[map == 0] <- NA
    
    # Put them together into the same layer
    map <- sum(map[[2]], (map[[3]] + 2), (map[[4]] + 4), (map[[5]] + 6), 
               (map[[6]] + 8), na.rm = TRUE)
    
    # Set to categorical
    map <- as.factor(map)
    levels(map) <- data.frame(
        value = 1:10, 
        usage = as.vector(
            sapply(c("Maize", "Paddy rice", "Sorghum", 
                     "Cassava", "Common beans"), 
                   function(x) sprintf("%s_%s", x, c("64", "80")))))
    
    # Get sub_map
    sub_ply <- ext(c(xmin = 35.1, xmax = 36.1, 
                     ymin = -4.2, ymax = -3.2)) %>% 
        vect(crs = "EPSG:4326") %>% st_as_sf(crs = 4326) %>% 
        st_centroid() %>% st_buffer(0.5)
    sub_map <- crop(map, sub_ply) %>% mask(sub_ply)
    cols_sub <- crop_cols[na.omit(freq(map)[["value"]]) %in% freq(sub_map)[["value"]]]
    zoomin <- ggplotGrob(
        ggplot() +
            geom_spatraster(data = sub_map, na.rm = TRUE) +
            scale_fill_manual(
                values = cols_sub,
                na.translate = FALSE) +
            geom_sf(data = sub_ply, color = "black", 
                    fill = "transparent", linewidth = 0.6) +
            coord_sf() + theme_void() +
            theme(text = element_text(size = 14),
                  legend.text = element_text(size = 14),
                  legend.spacing.y = unit(0.1, 'cm'),
                  legend.margin = margin(rep(0, 4)),
                  legend.key.height = unit(0.4, 'cm'),
                  legend.key.width = unit(0.8, 'cm'),
                  plot.margin = unit(c(0, -1, 0, -1), "cm"),
                  legend.position = "none"))
    
    # Figure
    ggplot() +
        geom_spatraster(data = map, na.rm = TRUE) +
        scale_fill_manual(
            values = crop_cols,
            na.translate = FALSE) +
        geom_sf(data = pas, color = "#f0f0f0", 
                fill = "#f0f0f0", linewidth = 0) +
        geom_sf(data = bry, color = "black", 
                fill = "transparent", linewidth = 0.6) +
        geom_sf(data = sub_ply, color = "black", 
                fill = "transparent", linewidth = 0.4) +
        annotation_custom(grob = zoomin,
                          xmin = 29.9, xmax = 33.1,
                          ymin = - 12.1, ymax = - 8.9) +
        coord_sf() + theme_void() +
        theme(text = element_text(size = 14),
              legend.text = element_text(size = 14),
              legend.spacing.y = unit(0.1, 'cm'),
              legend.margin = margin(rep(0, 4)), 
              legend.key.height = unit(0.4, 'cm'),
              legend.key.width = unit(0.8, 'cm'),
              plot.margin = unit(c(0, -1, 0, -1), "cm"),
              legend.position = "none")
})

g1 <- ggarrange(legend_lc, ggarrange(plotlist = fig_list, nrow = 1, ncol = 3,
                                     labels = c("BL", "Sce1", "Sce2"),
                                     font.label = list(size = 14),
                                     vjust = 3, hjust = 0), 
                nrow = 2, ncol = 1, heights = c(3, 10))

# New area and ecological costs change
vals <- read.csv(
    file.path("results/scenarios/summary", 
              "summary_scenarios_USGs.csv")) %>% 
    filter(scenario == "YBCaCoD") %>% 
    select(-c(scenario, num_unit, weight))
names(vals) <- c("B", "Ca", "Co", "D", "Y", "USG", "Practice")

vals <- vals %>% pivot_longer(
    c("B", "Ca", "Co", "D", "Y"), names_to = "Element")
vals <- vals %>% group_by(Element) %>% 
    reframe(USG = USG, Practice = Practice,
            value = value / max(value) * 100) %>% 
    mutate(USG = factor(USG, labels = c("64%", "80%"), levels = c(0.64, 0.8))) %>% 
    mutate(Practice = factor(Practice, labels = c("BL", "Sce1", "Sce2"),
                             levels = c("Y100", "CASS3", "Y140"))) %>% 
    mutate(Element = factor(
        Element, levels = c("Y", "B", "Ca", "Co", "D"),
        labels = c("New cropland area", "Biodiveristy loss", "Carbon loss",
                   "Connectivity loss", "Travel time")))
vals <- vals %>% mutate(group = ifelse(USG == "64%" & Practice == "BL", "B", "P"))

g2 <- ggplot(data = vals) + 
    geom_linerange(aes(
        x = Practice, ymin = 0, ymax = value, linetype = USG, color = group), 
        linewidth = 0.8, position = position_dodge2(width = 0.7)) +
    geom_point(aes(x = Practice, y = value, color = group), size = 2,
               position = position_dodge2(width = 0.7)) +
    geom_text(data = vals %>% mutate(value = ifelse(value == 100, NA, value)), 
              aes(x = Practice, y = value + 8, label = round(value, 1)), 
              size = 3, color = "black", fontface = "bold",
              position = position_dodge2(width = 0.7)) + 
    scale_color_manual(values = c("#d73027", "black")) +
    xlab("") + ylab("Percentage") + labs(linetype = "Cropland usage rate") +
    facet_grid(~Element) +
    theme_light() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(
              color = "lightgrey", linetype = "dashed"),
          panel.grid.minor.y = element_blank(),
          panel.border = element_blank(),
          strip.text.x = element_text(size = 14, color = "black"),
          axis.ticks.x = element_blank(),
          text = element_text(size = 14, color = "black"),
          axis.title = element_text(color = "black", size = 14),
          axis.text = element_text(color = "black", size = 14),
          axis.text.x = element_text(face = "bold"),
          legend.text = element_text(size = 14, color = "black"),
          legend.title = element_text(size = 14, color = "black"),
          legend.position = "top",
          legend.margin = margin(rep(0, 4)),
          plot.margin = unit(rep(0, 4), "cm")) +
    guides(color = "none")

ggarrange(g1, g2, nrow = 2, ncol = 1, heights = c(6, 5))

ggsave("figures/aag_fig4_4.png",
       bg = "white",
       width = 9.6, height = 6.4, dpi = 500)
