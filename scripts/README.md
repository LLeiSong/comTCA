## Prepare species and habitat

- Run `1_subset_expert_range_maps.R` to select species by taxon that are involved with Tanzania.
- Run `2_query_habitat_elevation.R` to query the habitat details and preferred elevation from ICUN Redlist.
- Run `3_refine_habitat_type_map.R` to refine the habitat type map by the high resolution land cover map.
- Run `4_prepare_elevation.R` to prepare the elevation layer from ALOS DSM layer.
- Run `5_rasterize_expert_ranges_scheduler.R` to rasterize all expert range polygons using HPC nodes. The script `5_rasterize_expert_ranges.R` can be used independently to rasterize the polygons.
- Run `6_refine_expert_ranges_scheduler.R` to refine the rasterized expert range maps using the refined habitat type map. Similarly, the script `6_refine_expert_ranges.R` can be used independently.
- Run `7_calc_global_range_size.py`, `7_calc_outside_range_size.py` and `7_calc_tz_range_size_scheduler.R` to calculate the habitat size inside and outside of Tanzania and global for area weighted biodiversity index.

## Biodiversity

- Run `8_1_resample_range_map_scheduler.R` to resample the Tanzania range maps because the original resolution is too computational resource consuming.
- Run `8_2_calc_species_richness.py` to calculate the different types of species richness index.
- Run `8_3_calc_biodiversity_index.R` to calculate the proactive biodiversity index.

## Habitat connectivity

- Run `9_1_select_core_species.R` to select core and representative species for connectivity simulation.
- Run `9_2_env_variables_cooking.R` to prepare image cube of environmental variables.
- Run `9_3_calc_species_suitability.R` to fit ensemble SDM using `biomod2` for the selected species and get the average suitability.
- Run `9_4_calc_connectivity.R` to simulate landscape connectivity using Circuitscape.

## Carbon cost

- Run the Rmarkdown `10_carbon_cost.qmd` to calculate the carbon density.

## Crop yield and gaps

- Run `11_downscale_yield.R` to use quantile regression forest to downscale SPAM 2017 yield maps of 10 km to 1 km.
- Run `11_calibrate_yield.R` to calibrate the downscaled yield to match with district surveys and then use predicted 95% quantile as the attainable yield.