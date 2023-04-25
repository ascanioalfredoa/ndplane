#### Set working directory ####
setwd("")

#### Load packages ####
library(virtualspecies)
library(tidyverse)
library(pracma)
library(doParallel)
source("rscripts/virtualsp_functions.R")
source("rscripts/BetaFunctions.R")

#### Get raster data ####
#worldclim <- getData("worldclim", var = "bio", res = 5)
if(!dir.exists("data/wc")) {
	dir.create("data/wc")
	P_url<- "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_5m_bio.zip"
	download.file(P_url, destfile="data/wc/wc2.1_5m_bio.zip")
	unzip("data/wc/wc2.1_5m_bio.zip", exdir = "data/wc/wc2.1_5m_bio")
}

worldclim <- list.files("data/wc/wc2.1_5m_bio", full.names = T)
worldclim <- stack(worldclim)
names(worldclim) <- gsub("wc2.1_5m|_", "", names(worldclim))
worldclim <- crop(worldclim, extent(-170, -55, 25, 70))

# Set parameter domains for changing response curve
a = b <- round(seq(10, 35, 1), 1)
alpha = gamma <- seq(1, 10, 3)
#alpha = gamma <- seq(0.1, 10, 0.5)

# Set parameter space (combination of parameters for any beta function)
par_space <- expand.grid(a = a, b = b, alpha = alpha, gamma = gamma)
par_space <- par_space[par_space$a < par_space$b, ]
#par_space$n <- 1:nrow(par_space)
gc()
write_csv(par_space, "output/virtualsp_Beta_Parameter_Space.csv")
