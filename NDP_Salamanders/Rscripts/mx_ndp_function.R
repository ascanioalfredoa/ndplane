maxent_ndp <- function(mx_sp1, mx_sp2, cut_off = 0.95) {
    
}
# Save species responses in different objects
A <- mar_res #%>% filter(grepl(taxa_comb[j, 1], category))
B <- spo_res #%>% filter(grepl(taxa_comb[j, 2], category))

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
    A[A[, 1] %in% B[, 1], ] %>%
    left_join(B[B[, 1] %in% A[, 1], ], by = "layer", suffix = c("_a", "_b")) %>%
    #rename(!!(.[["Source_a"]] %>% unique()) := y_a, !!(.[["Source_b"]] %>% unique()) := y_b) %>%
    dplyr::select(1, 2, 6) %>% 
    pivot_longer(cols = 2:3, names_to = "source_curve", values_to = "y") %>%
    group_by(layer) %>%
    summarise(y = min(y))

# Calculate area
NicheDiss <- 1 - (trapz(C[[1]], C[[2]])/trapz(A[[1]], A[[2]]) + trapz(C[[1]], C[[2]])/trapz(B[[1]], B[[2]]))/2

RNO <- 1 - (min(c(max(A[[1]]), max(B[[1]]))) - max(c(min(A[[1]]), min(B[[1]]))))/(max(c(max(A[[1]]), max(B[[1]]))) - min(c(min(A[[1]]), min(B[[1]]))))
if(RNO > 1) RNO <- 1
