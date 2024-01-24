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


# For each species, subset TGS that does not contain them
spo_tgs <- all_am %>%
    filter(!(Cmmn_Nm %in% "Spotted Salamander"))
mar_tgs <- all_am %>%
    filter(!(Cmmn_Nm %in% "Marbled Salamander"))

# If available, load pre-loaded thinned TGS
spo_tgs_t <- vect("Data/Ambystoma/spo_tgs_t.shp")
mar_tgs_t <- vect("Data/Ambystoma/mar_tgs_t.shp")

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

# Save TGS as background points
writeVector(mar_tgs_t, filename = "Data/Ambystoma/mar_tgs_t.shp")
writeVector(spo_tgs_t, filename = "Data/Ambystoma/spo_tgs_t.shp")

# Check thinned TGS points across the study area    
all_sa %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = spo_tgs_t, aes(color = Cmmn_Nm)) +
    xlim(ext(all_sa)[1:2]) + ylim(ext(all_sa)[3:4])

all_sa %>%
    ggplot() +
    geom_spatvector() +
    geom_spatvector(data = spo_tgs_t, aes(color = Cmmn_Nm)) +
    xlim(ext(all_sa)[1:2]) + ylim(ext(all_sa)[3:4])

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

#### Run univariate ENMTools using default enmtools.maxent configuration with appropriate background points####
#### Build loop around univariate maxent model selection to pick best model for our purposes, if possible ####
# Loop should go through each environmental variable (bio)

#### Create enmtools.species objects ####
mar_et <- enmtools.species(range = mask(chelsa_bioclim, mar_sa), 
                           presence.points = mar, 
                           background.points = mar_tgs_t, 
                           species.name = "Marbled")

spo_et <- enmtools.species(range = mask(chelsa_bioclim, spo_sa), 
                           presence.points = spo, 
                           background.points = spo_tgs_t, 
                           species.name = "Spotted")

ind_res <- NULL
plot_res <- NULL

for(i in 1:length(names(chelsa_bioclim))) {
    
#options(java.parameters = "-Xmx512m") #Default java memory
#options()$java.parameters
#options(java.parameters = "-Xmx36g") 

#### Running maxent models ####
mar_mx <- enmtools.maxent(mar_et, env = chelsa_bioclim[[i]], bg.source = "points")
spo_mx <- enmtools.maxent(spo_et, env = chelsa_bioclim[[i]], bg.source = "points")

#### Calculate D and I from maxent results ####
mx_sim <- raster.overlap(mar_mx, spo_mx) %>% do.call(cbind, .)

#### Calculate niche divergence plane indices from maxent response curves ####
ndp_sim <- maxent_ndp(mx_sp1 = mar_mx, mx_sp2 = spo_mx, cut_off = 0.95)

ind_res <- rbind(ind_res, cbind(ndp_sim$indices, mx_sim))
plot_res <- rbind(plot_res, ndp_sim$response_curves)
print(i)
gc()
}

write.csv(ind_res, "Results/NDP_salamander_values2.csv", row.names = F)
write_rds(plot_res, "Results/NDP_Salamanders_ResponseCurves2.rds")

