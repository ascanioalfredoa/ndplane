#### Functions to simulate I, D, and NDP using virtual species ####

#### Normalize suitability values from a virtualspecies object
normalize_suitability <- function(vsp) {
    # Where vsp = virtualspecies object
    vsp_ras <- vsp$suitab.raster
    norm_vsp <- vsp_ras/cellStats(vsp_ras, stat = sum)
    norm_vsp
}

#### Calculate Schoener's D and I from normalized suitability raster ####
calculate_similarity <- function(sp1, sp2, method = c("D", "I")) {
    if(!all(method %in% c("D", "I"))) {
        stop("Method for niche similarity calculation is not implemented or does not exist")
    }
    
    if("D" %in% method) {
        D <- 1 - 0.5*(cellStats(abs(sp1 - sp2), stat = sum))
    }
    
    if("I" %in% method) {
        I <- 1 - 0.5*sqrt(cellStats((sqrt(sp1) - sqrt(sp2))^2, stat = sum))
    }
    
    data.frame(D, I)

}


