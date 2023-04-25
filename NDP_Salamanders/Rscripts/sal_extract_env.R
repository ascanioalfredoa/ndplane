#### Load packages ####
library(terra)
library(tidyverse)
library(tidyterra)


##### Load thinned occurrences ####
sal_thin <- vect("Data/Ambystoma/Amb_macmar_thin.shp")

#### Load environmental chelsa data ####
chelsa_files <- list.files("Data/CHELSA_v2_1/", full.names = T)
chelsa_bioclim <- rast(chelsa_files)

sal_bioclim <- terra::extract(chelsa_bioclim, sal_thin)


#### Add environmental values to salamander data ####
sal <- 
    sal_thin %>%
    mutate(sal_bioclim[, -1])

as_tibble(sal) %>% View()

#### Load soil data ####
soils <- read_rds("Data/NATSGO/CONUS_gNATSGO_soilvars.rds")
#### Load soil raster ####
musoil_ras <- rast("Data/NATSGO/NATSGO_mapunit30m.tif") 

#### Extract map units from soil data into salamanders ####
sal_nad83 <- project(sal, musoil_ras)

musoil_sal <- terra::extract(musoil_ras, sal_nad83)

musoil_sal <- 
    musoil_sal %>%
    mutate(MUKEY = as.character(MUKEY)) %>%
    rename(mukey = MUKEY)

#### Concatenate soil data ####
lapply(soils, class) == "list"

soils_df <- soils[[1]]
for(i in 2:length(soils)) {
    
    if(all(class(soils[[i]]) != "list")) {
        soils_df <- 
            soils_df %>% left_join(soils[[i]], by = "mukey")
    } else {
        soil_var <- 
            soils[[i]][[1]] %>% 
            left_join(soils[[i]][[2]] %>% rename(!!names(soils[[i]][[1]])[2] := 'Value'), 
                      by = quo_name(names(soils[[i]][[1]])[2]))
        names(soil_var)[3] <- paste(names(soil_var)[2:3], collapse = "_")
        
        soils_df <- 
            soils_df %>%
            left_join(soil_var, by = "mukey")
    }
    
}

#### Match NATSGO data to our MUKEYS ####
musoil_sal <- 
    musoil_sal %>%
    left_join(soils_df, by = "mukey")

musoil_sal <- musoil_sal[!duplicated(musoil_sal$ID), ]


sal <- 
    sal %>%
    mutate(musoil_sal[, -c(1:2)])
#as_tibble(rattle) %>% View()

sal_df <- as_tibble(sal)

writeVector(sal, "Results/Amb_macopac_env.shp")
write_csv(sal_df, "Results/Amb_macopac_env.csv")
