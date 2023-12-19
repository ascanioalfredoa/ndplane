#### Load packages ####
library(terra)
library(tidyverse)
library(tidyterra)
library(ENMTools)

##### Load thinned occurrences ####
sal_thin <- vect("Data/Ambystoma/Amb_macmar_thin.shp")

#### Load environmental chelsa data ####
chelsa_files <- list.files("Data/CHELSA_v2_1/", full.names = T)
chelsa_bioclim <- rast(chelsa_files)

sal_bioclim <- terra::extract(chelsa_bioclim, sal_thin)

