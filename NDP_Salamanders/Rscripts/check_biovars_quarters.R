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
all_sa <- vect("Data/Ambystoma/All_StudyArea.shp")

#### Crop WC historical data and crop to study area for all years and months ####
wc_files <- list.files("/home/ascaniaa/Downloads/wc_historical/", pattern = ".tif", recursive = T, full.names = T)

r <- rast(wc_files)
r <- terra::crop(r, y = all_sa, mask = TRUE)
gc()
Sys.time()

#Rename layers by variable, month, and year
names(r) <- gsub("wc2.1_2.5m_", "", names(r))
names(r) <- paste(names(r), unlist(lapply(str_split(wc_files, pattern = "_"), function(x) str_split(x[[9]], "-")[[1]][[1]])), sep = "_")

#### Filter only years of interest and divide variables ####
prec <- r[[grepl("prec", names(r)) & grepl(paste0(1979:2013, collapse = "|"), names(r))]]
tmax <- r[[grepl("tmax", names(r)) & grepl(paste0(1979:2013, collapse = "|"), names(r))]]
tmin <- r[[grepl("tmin", names(r)) & grepl(paste0(1979:2013, collapse = "|"), names(r))]]

#### Loop over months to get the climatic normals 1979-2013 (same as CHELSA) ####
months <- lapply(str_split(names(prec), "_"), function(x) x[[2]]) %>% unlist() %>% as.integer()
months == lapply(str_split(names(tmin), "_"), function(x) x[[2]]) %>% unlist() %>% as.integer()
months == lapply(str_split(names(tmax), "_"), function(x) x[[2]]) %>% unlist() %>% as.integer()

ntmax = ntmin = nprec <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    nprec[[as.character(i)]] <- mean(prec[[months == i]])
    ntmax[[as.character(i)]] <- mean(tmax[[months == i]])
    ntmin[[as.character(i)]] <- mean(tmin[[months == i]])
}

ntmean <- (ntmax + ntmin)/2

#### Produce list of quarters through the year ####
quarters <- list()
for(i in 1:12) {
    quarters[[i]] <- ifelse(i:(i+2) <= 12, i:(i+2), i:(i+2) - 12)
}
quarters

#### Determine wettest/dryer quarter based on precipitation ####
wet <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    wet[[as.character(i)]] <- sum(nprec[[quarters[[i]]]])
}
wet

wettest <- which.max(wet)
dryest <- which.min(wet)

plot(wettest)
plot(dryest)

#### Determine average temperature of each quarter ####
nqtemp <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    nqtemp[[as.character(i)]] <- mean(ntmean[[quarters[[i]]]])
}


#### Mean temperature of wettest quarter ####
wqtemp <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    wqtemp[[as.character(i)]] <- nqtemp[[i]][wettest == i, drop = F]
}
wqtemp <- sum(wqtemp, na.rm = T)
plot(wqtemp)

chelsa_bio8 <- rast("Data/CHELSA_v2_1/bio8.tif")
chelsa_bio8 <- crop(chelsa_bio8, all_sa, mask = T)
plot(chelsa_bio8)

#### Mean temperature of driest quarter ####
dqtemp <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    dqtemp[[as.character(i)]] <- nqtemp[[i]][dryest == i, drop = F]
}
dqtemp <- sum(dqtemp, na.rm = T)
plot(dqtemp)

chelsa_bio9 <- rast("Data/CHELSA_v2_1/bio9.tif")
chelsa_bio9 <- crop(chelsa_bio9, all_sa, mask = T)
plot(chelsa_bio9)

################################################################################
################### Warmest/Coldest quarters ###################################
################################################################################
war <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    war[[as.character(i)]] <- mean(ntmean[[quarters[[i]]]])
}
war

warmest <- which.max(war)
coldest <- which.min(war)

plot(warmest)
plot(coldest)

#### Determine precipitation of each quarter ####
nqprec <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    nqprec[[as.character(i)]] <- sum(nprec[[quarters[[i]]]])
}

#### Precipitation of the warmest quarter ####
wqprec <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    wqprec[[as.character(i)]] <- nqprec[[i]][warmest == i, drop = F]
}
wqprec <- sum(wqprec, na.rm = T)
plot(wqprec)

chelsa_bio18 <- rast("Data/CHELSA_v2_1/bio18.tif")
chelsa_bio18 <- crop(chelsa_bio18, all_sa, mask = T)
plot(chelsa_bio18)

#### Precipitation of the warmest quarter ####
cqprec <- rast(ext(prec), resolution = res(prec), nlyrs = 12)
for(i in 1:12) {
    cqprec[[as.character(i)]] <- nqprec[[i]][coldest == i, drop = F]
}
cqprec <- sum(cqprec, na.rm = T)
plot(cqprec)

chelsa_bio19 <- rast("Data/CHELSA_v2_1/bio19.tif")
chelsa_bio19 <- crop(chelsa_bio19, all_sa, mask = T)
plot(chelsa_bio19)

library(RColorBrewer)

png(filename = "Figures/Check_biovars.png", width = 9000, height = 6000, res = 150)
par(mfrow = c(3, 4), cex = 2.5)
plot(as.factor(wettest), main = "Initial Month of Wettest quarter", col = brewer.pal(12, "Paired"))
plot(as.factor(dryest), main = "Initial Month of Driest quarter", col = brewer.pal(12, "Paired"))
plot(as.factor(warmest), main = "Initial Month of Warmest Quarter", col = brewer.pal(12, "Paired"))
plot(as.factor(coldest), main = "Initial Month of Coldest Quarter", col = brewer.pal(12, "Paired"))
plot(wqtemp, main = "BIO8: Mean Temperature of Wettest Quarter (Worldclim)")
plot(dqtemp, main = "BIO9: Mean Temperature of Driest Quarter (Worldclim)")
plot(wqprec, main = "BIO18: Precipitation of Warmest Quarter (Worldclim)")
plot(cqprec, main = "BIO19: Precipitation of Coldest Quarter (Worldclim)")
plot(chelsa_bio8, main = "BIO8: Mean Temperature of Wettest Quarter (CHELSA)")
plot(chelsa_bio9, main = "BIO9: Mean Temperature of Driest Quarter (CHELSA)")
plot(chelsa_bio18, main = "BIO18: Precipitation of Warmest Quarter (CHELSA)")
plot(chelsa_bio19, main = "BIO19: Precipitation of Coldest Quarter (CHELSA)")
dev.off()
