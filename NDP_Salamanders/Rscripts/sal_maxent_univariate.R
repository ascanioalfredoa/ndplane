#------------------------------------------------------------------------------#
########################### Load packages ######################################
#------------------------------------------------------------------------------#
library(terra)
library(tidyverse)
library(tidyterra)
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
all_amt <- 
    all_am %>% 
    #terra::split(f = "Cmmn_Nm") %>%
    lapply(FUN = function(x) thin_records(x, thin.par = 10, reps = 10)) #%>%
    #do.call(rbind, .)

all_sa %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = mar_tgs, aes(color = Cmmn_Nm)) +
    xlim(ext(all_sa)[1:2]) + ylim(ext(all_sa)[3:4])

# Save TGS as background points

#------------------------------------------------------------------------------#
################# Load and clean environmental chelsa data #####################
#------------------------------------------------------------------------------#
# Load bioclim rasters
chelsa_files <- list.files("Data/CHELSA_v2_1/", full.names = T)
chelsa_bioclim <- rast(chelsa_files)

# Load NATSGO variables?

# Do I need to resample NATSGO variables?

# Crop datasetets to the combination of study areas to reduce speed costs

#### Extract environmental values for occurrences
# Extract environmental values for each set of occurrences
sal_bioclim <- terra::extract(chelsa_bioclim, sal_thin)

#### Extract environmental values for background/TGS
# Extract environmental values for each set of background/TGS

#------------------------------------------------------------------------------#
############################ Running Models ####################################
#------------------------------------------------------------------------------#

#### Run univariate ENMeval using a combination of features up to cubic ####
# This would allow mean, variance, and skewness to be part of each response
# What should I do about the regularization multiplier? r.m.
# How do I pick the best univariate model?

#### Build loop around univariate maxent model selection to pick best model for our purposes, if possible ####
# Loop should go through each environmental variable (bio/NATSGO)


#### Check for unique environments using MOP ####
#Should this be univariate or with all the environmental vars?
# Both and see?