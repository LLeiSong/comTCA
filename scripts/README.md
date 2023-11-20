## Refine species habitat, expert range maps of species, and calculate related range sizes

1. Run `1_subset_expert_range_maps.R` to select species by taxon that are involved with Tanzania.
2. Run `2_query_habitat_elevation.R` to query the habitat details and preferred elevation from ICUN Redlist.
3. Run `3_refine_habitat_type_map.R` to refine the habitat type map by the high resolution land cover map.
4. Run `4_prepare_elevation.R` to prepare the elevation layer from ALOS DSM layer.
5. Run `5_rasterize_expert_ranges_scheduler.R` to rasterize all expert range polygons using HPC nodes. The script `5_rasterize_expert_ranges.R` can be used independently to rasterize the polygons.
6. Run `6_refine_expert_ranges_scheduler.R` to refine the rasterized expert range maps using the refined habitat type map. Similarly, the script `6_refine_expert_ranges.R` can be used independently.
7. Run `6_2_merge_birds.R` to merge the ranges of birds; resident, breeding, and non-breeding.
8. Run `7_calc_global_range_size.py`, `7_calc_outside_range_size.py` and `7_calc_tz_range_size_scheduler.R` to calculate the habitat size inside and outside of Tanzania and global for area weighted biodiversity index.
9. Finally run `7_make_species_info_catalog.R` to collect all numbers together into one file.

Note: the range sizes calculated from step 8 are all in pixels. Step 9 convert pixels into "square degree". So converting these values to square kilometers may have slight distortion.

## Biodiversity

1. Run `8_1_resample_range_map_scheduler.R` to resample the Tanzania range maps because the original resolution is too computational resource consuming.
2. Run `8_2_calc_species_richness.py` to calculate the different types of species richness index.
3. Run `8_3_calc_biodiversity_index.R` to calculate the proactive biodiversity index.

## Habitat connectivity

1. Run `9_1_select_core_species.R` to select core and representative species for connectivity simulation.
2. Run `9_2_env_variables_cooking.R` to prepare image cube of environmental variables.
3. Run `9_3_calc_species_suitability.R` to fit ensemble SDM using `biomod2` for the selected species and get the average suitability.
4. Run `9_4_calc_connectivity.R` to simulate landscape connectivity using Circuitscape.

## Carbon cost

- Run the Rmarkdown `10_carbon_cost.qmd` to calculate the carbon density.

## Crop yield and gaps

1. Run `11_downscale_yield.R` to use quantile regression forest to downscale SPAM 2017 yield maps of 10 km to 1 km.
2. Run `11_calibrate_yield.R` to calibrate the down-scaled yield to match with district surveys for current yields and then to match with FAO projections for attainable yields.

## Spatial tradeoff modeling for new cropland allocation

1. Run `12_1_prepare_inputs.R` to pick relevant layers from different directories and save to one path for convenience.
2. Run `12_2_prod_gain_from_ints.R` to first get the production growth on existing cropland, which will determine the needed area of extra cropland.
3. Run `12_3_prod_potential_farmable_land.R` (Optional) to get the production potential for all existing and farmable area.
4. Run `12_4_future_land_allocation.R` and various setting parameters to allocate land for future agriculture.
5. Run `12_5_gather_results.R` to collect and summarize the simulation results.

Note: Some scripts are designed to take command line inputs to conveniently submit to HPC for parallel computation. The slurm schedulers can be found in folder `~schedulers`. 

## Analysis

After all these steps, the users should get the results based on how they set. The tables/maps should locate in `~/results` with structure like:

- expansion (optional)
- intensification
- other folder named by you for simulations

The document to make figures for our study locates in `~docs` as `figures_in_manuscript.qmd`.
