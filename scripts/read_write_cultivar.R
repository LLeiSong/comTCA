## --------------------------------------------
## Script name: dssatr_write_cultigen and
## dssat_read_cultigen
## Purpose of script: write cultivar file
## for DSSAT.
## Author: Lei Song
## Date Created: 2023-01-24
##
## Copyright (c) Lei Song, 2023
## Email: lsong@clarku.edu
## --------------------------------------------
## Notes: Modified based on R package
## https://github.com/agroimpacts/rcropmod.
## Author: Lyndon Estes
## --------------------------------------------

# Read cultivar file
dssat_read_cultigen <- function(cultigen_file) {
    temp <- tempfile()
    on.exit(unlink(temp))
    
    # cat("\n",cultigen_file)
    readr::read_lines(cultigen_file) %>%
        gsub("!.*", "", .) %>%
        trimws("right") %>%
        magrittr::extract(. != "") %>%
        magrittr::extract(-1:-2) %>%
        readr::write_lines(temp)
    
    readr::read_fwf(
        temp, col_positions = readr::fwf_positions(
            start = c(1, 8, 25, 31, 38, 44, 50, 56, 62, 68),
            end =   c(6, 23, 29, 36, 42, 48, 54, 60, 66, 72),
            col_names = c("@VAR#",
                          "VRNAME..........",
                          "EXPNO",
                          "ECO#",
                          "P1",
                          "P2",
                          "P5",
                          "G2",
                          "G3",
                          "PHINT")),
        col_types = readr::cols(
            `@VAR#` = readr::col_character(),
            `VRNAME..........` = readr::col_character(),
            `EXPNO` = readr::col_character(),
            `ECO#` = readr::col_character(),
            P1 = readr::col_double(),
            P2 = readr::col_double(),
            P5 = readr::col_double(),
            G2 = readr::col_double(),
            G3 = readr::col_double(),
            PHINT = readr::col_double()
        ))
    
}

# Write cultivar file
dssatr_write_cultigen <- function(cultivars, 
                                  file.name = "NAMELESS", 
                                  dssat_path = NULL, 
                                  output_dir = ".") {
    # Must set dssat Data path
    if(is.null(dssat_path)) {
        stop("Must have a valid cultivar file for the reference.")
    }
    
    # Create dir if not exists
    if (! dir.exists(output_dir)) {
        dir.create(output_dir,
                   showWarnings = FALSE,
                   recursive = TRUE)}

    # Set file name
    file.name <- paste0(sprintf("%-8s", file.name))
    out.file <- file.path(output_dir, file.name)
    
    # Create template for the file
    template <- file.path(dssat_path, "Genotype/MZCER048.CUL")
    template <- readLines(template)
    
    headers <- template[1:50]
    
    lines <- paste0(
        format(cultivars$'@VAR#', width = 6)," ",
        format(cultivars$VRNAME.........., width = 16)," ",
        format(cultivars$EXPNO, width = 5)," ",
        format(cultivars$'ECO#', width = 6)," ",
        format(round(cultivars$P1 %>% as.numeric(), digits=1), width=5, digits=1, nsmall=1), " ",
        format(round(cultivars$P2 %>% as.numeric(), digits=1), width=5, digits=1, nsmall=1), " ",
        format(round(cultivars$P5 %>% as.numeric(), digits=1), width=5, digits=1, nsmall=1), " ",
        format(round(cultivars$G2 %>% as.numeric(), digits=1), width=5, digits=1, nsmall=1), " ",
        format(round(cultivars$G3 %>% as.numeric(), digits=1), width=5, digits=1, nsmall=1), " ",
        format(round(cultivars$PHINT %>% as.numeric(), digits=1), width=5, digits=1, nsmall=1))
    
    readr::write_lines(c(headers, lines), file = out.file)
    return(output_dir)
}
