#------------------------------------------------------------------------------#
########################### Load packages ######################################
#------------------------------------------------------------------------------#
library(terra)
library(tidyverse)
library(tidyterra)
library(maxnet)
library(ENMTools)
source("Rscripts/thin_records.R")

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


# For each species, subset TGS that does not contain them
spo_tgs <- all_am %>%
    filter(!(Cmmn_Nm %in% "Spotted Salamander"))
mar_tgs <- all_am %>%
    filter(!(Cmmn_Nm %in% "Marbled Salamander"))


#Thin TGS occurrences
spo_tgs_t <- spo_tgs %>%
    terra::split(f = "Cmmn_Nm") %>%
    lapply(FUN = function(x) thin_records(x, thin.par = 10, reps = 10)) %>%
    do.call(rbind, .)
spo_tgs_t <- thin_records(spo_tgs_t, thin.par = 10, reps = 10)

mar_tgs_t <- mar_tgs %>%
    terra::split(f = "Cmmn_Nm") %>%
    lapply(FUN = function(x) thin_records(x, thin.par = 10, reps = 10)) %>%
    do.call(rbind, .)
mar_tgs_t <- thin_records(mar_tgs_t, thin.par = 10, reps = 10)
    
all_sa %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = spo_tgs_t, aes(color = Cmmn_Nm)) +
    xlim(ext(all_sa)[1:2]) + ylim(ext(all_sa)[3:4])

# Save TGS as background points
writeVector(mar_tgs_t, filename = "Data/Ambystoma/mar_tgs_t.shp")
writeVector(spo_tgs_t, filename = "Data/Ambystoma/spo_tgs_t.shp")

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

#### Extract environmental values for occurrences
# Extract environmental values for each set of occurrences
sal_bioclim <- terra::extract(chelsa_bioclim, sal_thin)

#### Extract environmental values for background/TGS
# Extract environmental values for each set of background/TGS

#------------------------------------------------------------------------------#
############################ Running Models ####################################
#------------------------------------------------------------------------------#
#### Set up data for maxent ####
pa <- c(rep(1, length(mar)), rep(0, length(mar_tgs_t)))
envdat <- rbind(mar, mar_tgs_t)
envdat <- terra::extract(chelsa_bioclim, envdat)[, -1]

#Clean NAs and save pa and envdat again
data <- cbind(pa, envdat)
data <- data[complete.cases(data),]
pa <- data$pa
envdat <- data[, -1]

#### Run univariate ENMeval using a combination of features up to cubic ####
#Running maxnet
mdat <- data.frame(var = envdat[1], var2 = envdat[1]^2)
m <- maxnet(p = pa, data = mdat, f = maxnet.formula(p = pa, data = mdat, classes = "l"))
m1 <- maxnet(p = pa, data = envdat, f = formula(~ bio1 + I(bio1^2) - 1))

enmtools.maxent()
predict(m1, chelsa_bioclim)
# This would allow mean, variance, and skewness to be part of each response
# What should I do about the regularization multiplier? r.m.
# How do I pick the best univariate model?

#### Build loop around univariate maxent model selection to pick best model for our purposes, if possible ####
# Loop should go through each environmental variable (bio/NATSGO)

#------------------------------------------------------------------------------#
#################### Environmental Similarity - MOP ############################
#------------------------------------------------------------------------------#

#### Check for unique environments using MOP ####
#Should this be univariate or with all the environmental vars?
# Both and see?