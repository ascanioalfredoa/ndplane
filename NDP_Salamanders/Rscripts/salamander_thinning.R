#### Load packages ####
library(terra)
library(tidyverse)
library(tidyterra)
source("Rscripts/thin_records.R")

#### Load Data ####
sal <- vect("Data/Ambystoma/Full_WGS84.shp")

#### Filter A. opacum and A. maculatum ####
sal <- 
    sal %>%
    filter(Cmmn_Nm %in% c("Marbled Salamander", "Spotted Salamander"))

ggplot(sal) +
    geom_spatvector(aes(color = Cmmn_Nm))

#### Thin occurrences to 10km ####
spatdat2 <- 
    sal %>% 
    terra::split(f = "Cmmn_Nm") %>%
    lapply(FUN = function(x) thin_records(x, thin.par = 10, reps = 10)) %>%
    do.call(rbind, .)

writeVector(spatdat2, "Data/Ambystoma/Amb_macmar_thin.shp")

#### Plot resulting occurrences ####
ggplot(spatdat2) +
    geom_spatvector(aes(color = Cmmn_Nm))
