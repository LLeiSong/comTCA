# Function to allocate the land for agriculture
land_allocate <- function(
        potential_production, 
        inputs_std, 
        production_need, 
        cbetas,
        dst_dir,
        scenario){
    #  Calculate weights
    weights <- do.call(c, lapply(names(cbetas), function(name){
        inputs_std[[name]] * cbetas[name]
    })) %>% sum()
    
    # Get the decision table
    decision_mat <- values(c(weights, potential_production)) %>% 
        data.frame() %>% mutate(cell_id = 1:nrow(.)) %>% 
        na.omit() %>% arrange(-sum)
    
    # Match the target
    prod_to_catch <- production_need
    row_id <- 1
    while(prod_to_catch > 0){
        prod_to_catch <- prod_to_catch - decision_mat[row_id, 2]
        row_id <- row_id + 1
    }
    
    # Locate the land
    cell_ids <- decision_mat %>% slice(1:(row_id - 1)) %>% pull(cell_id)
    
    # Set pixels
    agro_land <- potential_production # as a template
    new_vals <- rep(0, ncell(potential_production))
    new_vals[cell_ids] <- 1
    values(agro_land) <- new_vals
    dst_fname <- file.path(
        dst_dir, sprintf("fut_agro_land_%s_equal.tif", scenario))
    writeRaster(agro_land, dst_fname)
}