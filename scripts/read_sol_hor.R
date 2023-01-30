## --------------------------------------------
## Script name: read_sol_hor
## Purpose of script: read soil layer statistics
## from SOL file used in DSSAT
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

library(data.table)

read_lines2 <- function(fname) {
    s <- file.info(fname)$size 
    buf <- readChar(fname, s, useBytes = TRUE)
    strsplit(buf,"\r\n", fixed = TRUE, useBytes = TRUE)[[1]]
}

read_sol_hor <- function(solfile, profiles) {
    sdat <- read_lines2(solfile)
    ndx <- which(nchar(sdat) == 1)
    idx <- which(Reduce("+", lapply(profiles, grepl, sdat, fixed = TRUE)) == 1)
    ndx <- ndx[ndx > min(idx)]
    ndx <- sapply(idx, function(x) ndx[ndx > x][1] - 1)
    pdat <- do.call(rbind.data.frame, lapply(1:length(idx), function(i) {
        p <- sdat[idx[i]:ndx[i]]
        ntier <- grep("@", p)  # index of SOL headers
        hind <- c(ntier[3] + 1, ifelse(length(ntier) == 4, ntier[4] - 1, 
                                       length(p)))
        hcol <- rbind("slb" = c(2, 6), "slll" = c(14, 18), "sdul" = c(20, 24))
        dul_ll0 <- sapply(1:nrow(hcol), function(x) {
            sapply(hind[1]:hind[2], function(y) {
                as.numeric(substr(p[y], hcol[x, 1], hcol[x, 2]))
            })
        })
        
        dul_ll <- matrix(dul_ll0,length(dul_ll0)/3, 3)
        colnames(dul_ll) <- c("SLB", "SDLL", "SDUL")
        dul_llv <- unlist(lapply(colnames(dul_ll), function(x) {
            d <- dul_ll[, x]
            c(d, rep(NA, 10 - length(d)))
        }))
        return(dul_llv)
    }))
    colnames(pdat) <- unlist(lapply(c("SLB", "SLLL", "SDUL"), paste0, 1:10))
    pdat <- as.data.table(pdat)
    pdat <- cbind("prof" = gsub("\\*", "", gsub(" .*.", "", profiles)), pdat)
    pdat[, prof := as.character(prof)]
    return(pdat)
}
