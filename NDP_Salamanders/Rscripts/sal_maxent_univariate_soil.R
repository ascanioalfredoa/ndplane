#------------------------------------------------------------------------------#
########################### Load packages ######################################
#------------------------------------------------------------------------------#
library(terra)
library(tidyverse)
library(tidyterra)
library(maxnet)
library(ENMTools)
source("Rscripts/thin_records.R")
source("Rscripts/enmtools_helper_functions.R")
source("Rscripts/mx_ndp_function.R")
#------------------------------------------------------------------------------#
#################### Load and clean Occurrences and basemap ####################
#------------------------------------------------------------------------------#

##### Load thinned occurrences ####
sal_thin <- vect("Data/Ambystoma/Amb_macmar_thin.shp")

sal_thin %>%
    ggplot(aes(color = Cmmn_Nm)) +
    geom_spatvector()

#### Load north-american basemap ####
na <- vect("Data/Basemaps/northamerica/World_Countries_Generalized.shp")
#na <- terra::project(na, sal_thin)
na <- na %>%
    filter(ISO %in% c("US", "CA", "MX"))
na <- terra::aggregate(na)

na %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = sal_thin, aes(color = Cmmn_Nm))

#### Divide occurrences by species ####
mar <- sal_thin %>%
    filter(Cmmn_Nm == "Marbled Salamander")
spo <- sal_thin %>%
    filter(Cmmn_Nm == "Spotted Salamander")

#------------------------------------------------------------------------------#
########################### Create study areas #################################
#------------------------------------------------------------------------------#

#### Create study area around salamander occurrences ####
#### There should be two study areas at the end (marbled vs spotted)

# Create buffer
mar_sa <- buffer(mar, 1e5) #100 km buffer (10k meters) for study area
spo_sa <- buffer(spo, 1e5)

# Union/aggregate of buffers
mar_sa <- terra::aggregate(mar_sa)
spo_sa <- terra::aggregate(spo_sa)

# Intersect study areas with north american boundaries
mar_sa <- terra::intersect(mar_sa, na)
spo_sa <- terra::intersect(spo_sa, na)

# Create "total study area" with the combination of both study areas to ease some later processes down the line
all_sa <- terra::union(mar_sa, spo_sa) %>% terra::aggregate()

all_sa %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = sal_thin, aes(color = Cmmn_Nm))

#------------------------------------------------------------------------------#
######### Maybe load full salamander dataset to build TGS background ###########
#------------------------------------------------------------------------------#

# Load all ambystoma occurrences
all_am <- vect("Data/Ambystoma/Full_WGS84.shp")

# Filter by study areas
all_am <- terra::intersect(all_am, all_sa)


# # For each species, subset TGS that does not contain them
# spo_tgs <- all_am %>%
#     filter(!(Cmmn_Nm %in% "Spotted Salamander"))
# mar_tgs <- all_am %>%
#     filter(!(Cmmn_Nm %in% "Marbled Salamander"))

# If available, load pre-loaded thinned TGS
spo_tgs_t <- vect("Data/Ambystoma/spo_tgs_t.shp")
mar_tgs_t <- vect("Data/Ambystoma/mar_tgs_t.shp")

# Check thinned TGS points across the study area    
spo_sa %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = spo_tgs_t, aes(color = Cmmn_Nm)) +
    xlim(ext(all_sa)[1:2]) + ylim(ext(all_sa)[3:4])

mar_sa %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = mar_tgs_t, aes(color = Cmmn_Nm)) +
    xlim(ext(all_sa)[1:2]) + ylim(ext(all_sa)[3:4])


#Find soil files
soil_names <- list.files("../../../Ambystoma_NDP_studyarea/", pattern = ".tif$")
soil_files <- list.files("../../../Ambystoma_NDP_studyarea/", pattern = ".tif$", full.names = T)

#Load soil rasters
soil_ras <- rast(soil_files)

#Clean soil names
names(soil_ras) <- gsub("SoilRas_|_30Meter", "", names(soil_ras))

#Project all salamander study areas and points into NAD83 (matching soils)
all_sa_nad83 <- project(all_sa, soil_ras)
mar_sa_nad83 <- project(mar_sa, soil_ras)
spo_sa_nad83 <- project(spo_sa, soil_ras)
mar_nad83 <- project(mar, soil_ras)
spo_nad83 <- project(spo, soil_ras)
mar_tgs_t_nad83 <- project(mar_tgs_t, soil_ras)
spo_tgs_t_nad83 <- project(spo_tgs_t, soil_ras)


#soil_ras <- crop(soil_ras, all_sa_nad83, mask = TRUE)

#terra::writeRaster(soil_ras, "Data/NATSGO/soil_ras_crop.tif")

#soil_ras <- rast("Data/NATSGO/soil_ras_crop.tif")

#### Aggregate soil rasters --- Too heavy to run as they are ####
# dir.create("Data/NATSGO/Aggregated34")
# Sys.time()
# for(i in 1:length(names(soil_ras))) {
#     soil_ras_ag34_mean <- terra::aggregate(soil_ras[[i]], fact = 34, fun = mean, na.rm = T)
#     terra::writeRaster(soil_ras_ag34_mean, paste("Data/NATSGO/Aggregated34/", names(soil_ras_ag34_mean), ".tif", sep = ""))
#     print(i)
#     gc()
# }
# Sys.time()
# 
# dir.create("Data/NATSGO/Aggregated26")
# Sys.time()
# for(i in 1:length(names(soil_ras))) {
#     soil_ras_ag34_mean <- terra::aggregate(soil_ras[[i]], fact = 26, fun = mean, na.rm = T)
#     terra::writeRaster(soil_ras_ag34_mean, paste("Data/NATSGO/Aggregated26/", names(soil_ras_ag34_mean), ".tif", sep = ""))
#     print(i)
#     gc()
# }
# Sys.time()
# 
# dir.create("Data/NATSGO/Aggregated13")
# Sys.time()
# for(i in 1:length(names(soil_ras))) {
#     soil_ras_ag34_mean <- terra::aggregate(soil_ras[[i]], fact = 13, fun = mean, na.rm = T)
#     terra::writeRaster(soil_ras_ag34_mean, paste("Data/NATSGO/Aggregated13/", names(soil_ras_ag34_mean), ".tif", sep = ""))
#     print(i)
#     gc()
# }
# Sys.time()

soil_ras_ag34_mean <- rast(list.files("Data/NATSGO/Aggregated13", full.names = T, pattern = ".tif$"))

mar_et_n83 <- enmtools.species(range = crop(soil_ras_ag34_mean, mar_sa_nad83, mask = TRUE), 
                               presence.points = mar_nad83, 
                               background.points = mar_tgs_t_nad83, 
                               species.name = "Marbled")
gc()
spo_et_n83 <- enmtools.species(range = crop(soil_ras_ag34_mean, spo_sa_nad83, mask = TRUE), 
                               presence.points = spo_nad83, 
                               background.points = spo_tgs_t_nad83, 
                               species.name = "Spotted")
gc()
ind_res <- NULL
plot_res <- NULL

for(i in 1:length(names(soil_ras_ag34_mean))) {
    #Given the size of the soil rasters, I'll load one by one and delete from memory at the end
    #options(java.parameters = "-Xmx512m") #Default java memory
    #options()$java.parameters
    options(java.parameters = "-Xmx100g") 
    
    #### Running maxent models ####
    mar_mx <- enmtools.maxent(mar_et_n83, env = soil_ras_ag34_mean[[i]], bg.source = "points")
    gc()
    spo_mx <- enmtools.maxent(spo_et_n83, env = soil_ras_ag34_mean[[i]], bg.source = "points")
    gc()
    
    #### Calculate D and I from maxent results ####
    mx_sim <- raster.overlap(mar_mx, spo_mx) %>% do.call(cbind, .)
    
    #### Calculate niche divergence plane indices from maxent response curves ####
    ndp_sim <- maxent_ndp(mx_sp1 = mar_mx, mx_sp2 = spo_mx, cut_off = 0.95)
    
    ind_res <- rbind(ind_res, cbind(ndp_sim$indices, mx_sim))
    plot_res <- rbind(plot_res, ndp_sim$response_curves)
    
    gc()
    print(i)
}

write.csv(ind_res, "Results/NDP_salamander_values_soil_ag13.csv")
write_rds(plot_res, "Results/NDP_Salamanders_ResponseCurves_soil_ag13.rds")

sys.time()
#------------------------------------------------------------------------------#
#################### Environmental Similarity - MOP ############################
#------------------------------------------------------------------------------#

#### Check for unique environments using MOP ####
#Should this be univariate or with all the environmental vars?
# Both and see?
