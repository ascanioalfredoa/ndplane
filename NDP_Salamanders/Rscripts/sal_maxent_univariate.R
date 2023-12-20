#### Load packages ####
library(terra)
library(tidyverse)
library(tidyterra)
library(ENMTools)

##### Load thinned occurrences ####
sal_thin <- vect("Data/Ambystoma/Amb_macmar_thin.shp")

#### Create study area around salamander occurrences ####
#### There should be two study areas at the end (marbled vs spotted)

# Create buffer

# Union/aggregate of buffers

# Dissolve bounderies

# Create "total study area" with the combination of both study areas to ease some later processes down the line

#### Maybe load full salamander dataset to build TGS background ####
# Load all ambystoma occurrences

# Filter by study areas

# For each species, subset TGS that does not contain them
# Check whether TGS occurrences were previously thinned, if not, thin them. 

# Save TGS as background points

#### Load and clean environmental chelsa data ####
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

#### Run univariate ENMeval using a combination of features up to cubic ####
# This would allow mean, variance, and skewness to be part of each response
# What should I do about the regularization multiplier? r.m.
# How do I pick the best univariate model?

#### Build loop around univariate maxent model selection to pick best model for our purposes, if possible ####
# Loop should go through each environmental variable (bio/NATSGO)


#### Check for unique environments using MOP ####
#Should this be univariate or with all the environmental vars?
# Both and see?