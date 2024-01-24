#------------------------------------------------------------------------------#
########################### Load packages ######################################
#------------------------------------------------------------------------------#
library(terra)
library(tidyverse)
library(tidyterra)

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
################# Load and clean environmental chelsa data #####################
#------------------------------------------------------------------------------#
# Load bioclim rasters
chelsa_files <- list.files("Data/CHELSA_v2_1/", full.names = T)
chelsa_bioclim <- rast(chelsa_files)

# Load NATSGO variables?
# Do I need to resample NATSGO variables?

# Crop datasetets to the combination of study areas to reduce speed costs
chelsa_bioclim <- crop(chelsa_bioclim, all_sa)
chelsa_bioclim <- mask(chelsa_bioclim, all_sa)
#plot(chelsa_bioclim)

#Find soil files
soil_names <- list.files("../../../Ambystoma_NDP_studyarea/", pattern = ".tif$")
soil_files <- list.files("../../../Ambystoma_NDP_studyarea/", pattern = ".tif$", full.names = T)

#Load soil rasters
soil_ras <- rast(soil_files)

nad83_chelsa_sa <- project(chelsa_bioclim[[1]], crs(soil_ras))

terra::writeRaster(nad83_chelsa_sa, "Data/nad83_chelsa_sa.tif")

soil_ras_resampled <- terra::resample(x = soil_ras, y = nad83_chelsa_sa, threads = T)

terra::writeRaster(soil_ras_resampled, "Data/soilras_chelsa_resampled_bilinear_nad83.tif")
