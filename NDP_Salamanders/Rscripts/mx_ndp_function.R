maxent_ndp <- function(mx_sp1, mx_sp2, cut_off = 0.95) {
    library(pracma)
    library(zoo)
    # Save variable name
    var <- unique(names(mx_sp1$response.plots)[[1]], names(mx_sp1$response.plots)[[1]])
    
    # Save species responses in different objects
    A <- mx_sp1$response.plots[[1]]$data %>% filter(source == "Suitability") #%>% filter(grepl(taxa_comb[j, 1], category))
    B <- mx_sp2$response.plots[[1]]$data %>% filter(source == "Suitability") #%>% filter(grepl(taxa_comb[j, 2], category))
    
    A$source <- mx_sp1$species.name
    B$source <- mx_sp2$species.name
    
    # Recalculate original probability value to assess quantiles
    A$cs <- cumsum(A$value/sum(A$value))
    B$cs <- cumsum(B$value/sum(B$value))
    
    # Cut responses given the cut_off threshold
    lco <- (1-cut_off)/2
    uco <- 1-((1-cut_off)/2)
    
    A <- A[A$cs >= lco & A$cs <= uco, ]
    B <- B[B$cs >= lco & B$cs <= uco, ]
    #### Integration by determining the lower function ####
    # Where the values of X are shared between distribution
    # And where the values of Y are the lowest between the two possibilities
    
    C <- 
        #A[A[, 1] %in% B[, 1], ] %>%
        A[A[, 1] >= min(B[, 1]) & A[, 1] <= max(B[, 1]), 1:3] %>%
        #left_join(B[B[, 1] %in% A[, 1], ], by = "layer", suffix = c("_a", "_b")) %>%
        full_join(B[B[, 1] >= min(A[, 1]) & B[, 1] <= max(A[, 1]), 1:3], by = "layer") %>% 
        #rbind(B[B[, 1] > min(A[, 1]) & B[, 1] < max(A[, 1]), ]) %>% 
        arrange(layer) %>% 
        mutate(value.x = na.approx(value.x, x = layer, na.rm = F)) %>%
        mutate(value.y = na.approx(value.y, x = layer, na.rm = F)) %>% 
        dplyr::select(1, 2, 4) %>%
        filter(complete.cases(.)) %>%
        #rename(!!(.[["Source_a"]] %>% unique()) := y_a, !!(.[["Source_b"]] %>% unique()) := y_b) %>%
        #dplyr::select(1, 2, 5) %>% 
        pivot_longer(cols = 2:3, names_to = "source_curve", values_to = "y") %>%
        group_by(layer) %>%
        summarise(y = min(y))
    
    # Calculate area
    NicheDiss <- 1 - (trapz(C[[1]], C[[2]])/trapz(A[[1]], A[[2]]) + trapz(C[[1]], C[[2]])/trapz(B[[1]], B[[2]]))/2
    
    NicheEx <- 1 - (min(c(max(A[[1]]), max(B[[1]]))) - max(c(min(A[[1]]), min(B[[1]]))))/(max(c(max(A[[1]]), max(B[[1]]))) - min(c(min(A[[1]]), min(B[[1]]))))
    if(NicheEx > 1) NicheEx <- 1
    
    res <- list()
    res[["indices"]] <- data.frame(dissimilarity = NicheDiss, exclusivity = NicheEx, var = var)
    res[["response_curves"]] <- rbind(A, B)
    res[["response_curves"]]$var <- var
    
    return(res)
}
