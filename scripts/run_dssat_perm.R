## --------------------------------------------
## Script name: run_dssat_perm
## Purpose of script: function to run DSSAT
## experiments.
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Notes: Make sure you install DSSAT and set 
## up properly.
## --------------------------------------------

# Load libraries
library(Dasst)
library(DSSAT)
library(dplyr)

exec_csm <- function(projdir, csmdir, rundir, btype = "B", 
                     csm = "dscsm048", 
                     bname = "R_DSSBatch.v48") {
    bcmd <- sprintf("%s/%s %s %s", csmdir, csm, btype, bname)
    setwd(rundir)
    system(bcmd, ignore.stdout = TRUE)
    setwd(projdir)
}

batch_file <- function(xl, xfiles, outdir, btype, bname = "R_DSSBatch", RP = 1, 
                       SQ = 0, OP = 0, CO = 0) {
    outname <- paste0(bname, ".v48")
    # i <- 1
    # xfiles <- unlist(xfiles)
    header <- rbind(paste0("$BATCH(", btype, ")"), "!", 
                    sprintf("%6s %92s %6s %6s %6s %6s", "@FILEX", "TRTNO",
                            "RP", "SQ", "OP", "CO"))
    guts <- do.call(rbind, lapply(1:length(xl), function(i) {
        trt <- xl[[i]]$N  # pull out treatment numbers
        cbind(sprintf("%6s %86s %6i %6i %6i %6i", xfiles[i], trt, RP, SQ, OP, CO))
    }))  
    outbatch <- rbind(header, guts)
    # outbatch <- rbind(rbind(paste0("$BATCH(", btype, ")"), "!",
    #                       sprintf("%6s %92s %6s %6s %6s %6s", "@FILEX", "TRTNO",
    #                                "RP", "SQ", "OP", "CO")),
    #                          cbind(sprintf("%6s %86s %6i %6i %6i %6i", xfiles, 1,
    #                                      RP, SQ, OP, CO)))
    
    # Write the batch file to the selected folder  
    write(outbatch, fp(outdir, outname), append = FALSE)
    return(outname)
}

run_dssat_perm <- function(crop = "MZ", 
                           plant_density = 3.7, 
                           row_spacing = 90, 
                           N_kg = 5, h20_s = 0.05, 
                           soilid, 
                           weatherid,
                           start_year, 
                           end_year, 
                           planting_date,
                           coords, 
                           cultivars,
                           dssat_path,
                           proj_path) {
    # Set initial date 30 days before the planting date
    initial_date <- planting_date - 7
    
    # Extract soil info
    profdat <- read_sol_hor(
        solfile = file.path(dssat_path, "Soil/SOIL.SOL"),
        profiles = soilid)
    
    
    ## Make fields for x_tab function
    # combine WTH name, coordinates, soil grid id, and profile data
    # Create unique field ID (FID) and name for output X file
    # Fixed parameters
    message("Step1 -- Generate field file.")
    field <- cbind("WTH" = weatherid, 
                   lat = coords$lat, lon = coords$lon, 
                   soilid, profdat) %>%
        mutate(ID_FIELD = fid(WTH), XNAME = xname(WTH, length(WTH))) %>%
        mutate(CR = crop, PPOP = plant_density, PDATE = planting_date, 
               SDATE = initial_date, ICDAT = initial_date, PLRS = row_spacing, 
               FAMN = N_kg, H20_s = h20_s, NYERS = (end_year - start_year + 1))
    xtab <- x_tab(fields = field)
    
    message("Step2 -- Build cultivars.")
    cult <- as.vector(cultivars)
    
    message("Step3 -- Build tcomb.")
    tcomb <- expand.grid(list("PDATE" = toString(planting_date), 
                              "INGENO" = cult),
                         stringsAsFactors = FALSE)
    
    message("Step4 -- Build ttab.")
    ttab <- cbind("N" = 1:nrow(tcomb), tcomb,
                  t_tab(tvars = c("PDATE", "INGENO"), 
                        topts = c("MP", "CU"),
                        ttab = tcomb))
    
    message("Step5 -- Build xtabl.")
    xtabl <- lapply(1:nrow(xtab), function(j) { # j <- 1
        d <- xtab[j, ]
        d <- do.call(rbind, lapply(1:nrow(ttab), function(x) d))  # expand rows
        upd_col <- intersect(colnames(d), colnames(ttab))
        d[, c(upd_col) := ttab[, upd_col]]  # update columns having variable values
        xt <- cbind(data.table(ttab), d[, !colnames(d) %in% upd_col, with = FALSE])
        xt %>% mutate(TNAME = INGENO, VBOSE = "A")
    })
    
    message("Step6 -- Build xrun.")
    rundir <- file.path(dssat_path, "Maize")
    
    do.call(rbind, lapply(xtabl, function(xf) {
        xfnm <- x_file(xtab = xf, outdir = rundir, 
                       z = "01", xtype = ".MZX")
        bname <- batch_file(xl = list(xf), xfiles = xfnm,
                            outdir = rundir, btype = "MAIZE")
        exec_csm(projdir = proj_path, csm = "dscsm048", 
                 csmdir = dssat_path, rundir = rundir, 
                 bname = bname)
         read_output(file.path(dssat_path, "Maize", "Summary.OUT"))
    })) %>% select(c("EXNAME", "HWAH"))
}

