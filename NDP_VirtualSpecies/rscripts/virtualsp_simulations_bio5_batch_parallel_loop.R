#### Set working directory ####
#setwd("")

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
worldclim <- rast(worldclim)
names(worldclim) <- gsub("wc2.1_5m|_", "", names(worldclim))
worldclim <- crop(worldclim, ext(-170, -55, 25, 70))

#### Read parameter space ####
par_space <- read_csv("output/virtualsp_Beta_Parameter_Space.csv")

#### Apply beta functions to parameter space ####
sp1 = sp2 <- 1:nrow(par_space)

results_names <- c("spA_par", "spB_par", "Dissimilarity", "Exclusivity", "D", "I",
                   paste(names(par_space), "_A", sep = ""),
                   paste(names(par_space), "_B", sep = ""),
                   "ND_magnitude", "ND_angle"
)

results <- as.data.frame(array(data = 0, 
                               dim = c(1, length(results_names)), 
                               dimnames = list(1, results_names)))
result_k <- results
k = 1

#### Set up parallel run ####
no_cores <- detectCores() - 4
#cl <- makeCluster(no_cores, type="FORK") #Fork doesn't work in Windows
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

dir.create("output/virtualspecies_bio5_sim_output")

########################### LOOP ##########################################
#no_nodes <- 8
#node_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#node_i

current_time <- Sys.time()
#foreach(i = sp1[round(seq(node_i, nrow(par_space), no_nodes))], .packages = c("tidyverse", "pracma", "virtualspecies"))  %dopar% {
#foreach(i = sp1, .packages = c("tidyverse", "pracma", "virtualspecies"))  %dopar% {
foreach(i = 1810:1818, .packages = c("tidyverse", "pracma", "virtualspecies"))  %dopar% {
    #### Set working directory ####
    #setwd("")    
    source("rscripts/BetaFunctions.R")
    source("rscripts/virtualsp_functions.R")
    
    #### Terra does not allow paralellization of SpatRasters or SpatVectors
    #I'll load the rasters inside the foreach loop
    worldclim <- list.files("data/wc/wc2.1_5m_bio", full.names = T)
    worldclim <- rast(worldclim)
    names(worldclim) <- gsub("wc2.1_5m|_", "", names(worldclim))
    worldclim <- crop(worldclim, ext(-170, -55, 25, 70))
    
    #####
    file <- paste("output/virtualspecies_bio5_sim_output/BetaSimulation_Results", "_spA_", i, ".csv", sep = "")
    
    #### Generate virtual species 1 ####
    # Define function parameters
    sp1.parameters <- formatFunctions(bio5 = c(fun = 'betaFun', 
                                               p1 = par_space$a[i], 
                                               p2 = par_space$b[i], 
                                               alpha = par_space$alpha[i], 
                                               gamma = par_space$gamma[i])#,
                                      #bio7 = c(fun = 'betaFun', p1 = 70, p2 = 650, alpha = 3, gamma = 3),
                                      #bio12 = c(fun = 'betaFun', p1 = 500, p2 = 5000, alpha = 3, gamma = 3)
    )
    
    
    #vsp1 <- generateSpFromFun(raster.stack = worldclim[[c("bio5", "bio12")]],
    vsp1 <- generateSpFromFun(raster.stack = worldclim[[c("bio5")]],
                              parameters = sp1.parameters,
                              #formula = "bio5 + bio12",
                              formula = "bio5",
                              plot = F) # Change to false once loop is ready
    
    vsp1 <- normalize_suitability(vsp1)
    
    write_csv(results[0, ], file)
    if(i > 1) {
        sp2 <- sp1
        sp2 <- sp2[-c(1:(i-1))]
    }
    for(j in sp2) {
        
        if(par_space$a[j] > par_space$b[i]) {next}
        if(par_space$a[i] > par_space$b[j]) {next}
        
        #### Generate virtual species 2 ####
        # Define function parameters
        sp2.parameters <- formatFunctions(bio5 = c(fun = 'betaFun', 
                                                   p1 = par_space$a[j], 
                                                   p2 = par_space$b[j], 
                                                   alpha = par_space$alpha[j], 
                                                   gamma = par_space$gamma[j])#,
                                          #bio7 = c(fun = 'betaFun', p1 = 70, p2 = 650, alpha = 3, gamma = 3),
                                          #bio12 = c(fun = 'betaFun', p1 = 500, p2 = 5000, alpha = 3, gamma = 3) #change to 2000-5000 for sp2
        )
        
        
        #vsp2 <- generateSpFromFun(raster.stack = worldclim[[c("bio5", "bio12")]],
        vsp2 <- generateSpFromFun(raster.stack = worldclim[[c("bio5")]],
                                  parameters = sp2.parameters,
                                  #formula = "bio5 + bio12",
                                  formula = "bio5",
                                  plot = F) # Change to false once loop is ready
        
        vsp2 <- normalize_suitability(vsp2)
        
        sp_similarity <- calculate_similarity(vsp1, vsp2)
        
        responses <- betaPDF_pair(a1 = par_space$a[i],
                                  b1 = par_space$b[i],
                                  alpha1 = par_space$alpha[i],
                                  gamma1 = par_space$gamma[i],
                                  a2 = par_space$a[j],
                                  b2 = par_space$b[j],
                                  alpha2 = par_space$alpha[j],
                                  gamma2 = par_space$gamma[j],
                                  pden = 3000)
        
        A <- responses[, c(1, 2)][responses[, 2] > 0, ]
        names(A)[2] <- "y"
        A <- unique(A)
        B <- responses[, c(1, 3)][responses[, 3] > 0, ]
        names(B)[2] <- "y"
        B <- unique(B)
        
        C <- 
            A[A$x %in% B$x, ] %>%
            left_join(B[B$x %in% A$x, ], by = "x", suffix = c("_a", "_b")) %>%
            pivot_longer(cols = 2:3, names_to = "source_curve", values_to = "y") %>%
            group_by(x) %>%
            summarise(y = min(y))
        
        NicheDiss <- 1 - (trapz(C$x, C$y)/trapz(A$x, A$y) + trapz(C$x, C$y)/trapz(B$x, B$y))/2
        NicheDiss <- round(NicheDiss, 3)
        
        RNO <- 1 - (min(c(max(A$x), max(B$x))) - max(c(min(A$x), min(B$x))))/(max(c(max(A$x), max(B$x))) - min(c(min(A$x), min(B$x))))
        if(RNO > 1) RNO <- 1
        RNO <- round(RNO, 3)
        
        
        ND_magnitude <- sqrt((NicheDiss^2) + (RNO^2))
        ND_angle <- atan2(RNO, NicheDiss)*(180/pi)
        
        result_k[1, 1:4] <- c(i, j, NicheDiss, RNO)
        result_k[1, 5:6] <- sp_similarity
        result_k[1, 7:10] <- par_space[i, ]
        result_k[1, 11:14] <- par_space[j, ]
        result_k[1, 15] <- ND_magnitude
        result_k[1, 16] <- ND_angle
        
        
        
        write_csv(result_k, file, append = TRUE)
        
        #if(k %% 100 == 0) {
        #Elapsed_time <- Sys.time() - current_time
        #ETA <- (Elapsed_time/k)*length(sp2)
        
        #RTA <- ETA-Elapsed_time
        #print(paste("spA =", i, "| spB =", j, "Init_Time =", current_time, 
        #            "Elapsed_Time =", Elapsed_time, "ETA =", ETA, "RTA =", RTA))
        
        #}
        #k = k + 1
        
        gc()
    }
    file.remove(dir(tempdir(), full.names=TRUE))
}
final_time <- Sys.time()
stopCluster(cl)


