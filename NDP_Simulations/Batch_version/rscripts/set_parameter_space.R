#### Set working directory ####
setwd("/home/ascaniaa/Documents/NDP_Simulation_batch")

#### Generate parameter space ####
library(tidyverse)
library(pracma)
library(doParallel)
source("rscripts/BetaFunctions.R")

# Set parameter domains
#General cover of the niche space
alpha = gamma <- round(c(1, seq(2, 12, 2)), 3)
a = b <- round(seq(0, 1, 0.1), 3)
A <- expand.grid(a, b, alpha, gamma)

#Cover right end of the niche space - Almost fully exclusive niches
#a = b <- round(seq(0, 0.025, 0.005), 3)
a = b <- round(seq(0, 0.025, 0.005), 3)
B <- expand.grid(a, b, alpha, gamma)

#Cover left end of the niche space - Almost fully inclusive niches
a = round(seq(0, 0.1, 0.025), 3)
b = round(seq(0.9, 1, 0.025), 3)
C <- expand.grid(a, b, alpha, gamma)

#a = b <- round(seq(0, 1, 0.1), 3)
#alpha = gamma <- round(c(1, seq(2, 12, 2)), 3)
#alpha = gamma <- seq(0.1, 10, 0.5)

# Set parameter space (combination of parameters for any beta function)
par_space <- rbind(A, B, C)
names(par_space) <- c("a", "b", "alpha", "gamma")
par_space <- par_space[par_space$a < par_space$b, ]
gc()
write_csv(par_space, "output/Beta_Parameter_Space.csv")
