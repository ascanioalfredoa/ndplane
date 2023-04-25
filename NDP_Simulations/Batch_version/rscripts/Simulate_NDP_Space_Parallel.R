################### Niche Divergence Plane Simulations #########################

#### Set working directory ####
setwd("/home/ascaniaa/Documents/NDP_Simulation_batch")

#### Generate parameter space ####
library(tidyverse)
library(pracma)
library(doParallel)
source("rscripts/BetaFunctions.R")

#### Read parameter space ####
par_space <- read_csv("output/Beta_Parameter_Space.csv")

#### Apply beta functions to parameter space ####
sp1 = sp2 <- 1:nrow(par_space)

results_names <- c("spA_par", "spB_par", "Dissimilarity", "Exclusivity",
                   paste(names(par_space), "_A", sep = ""),
                   paste(names(par_space), "_B", sep = "")
)

results <- as.data.frame(array(data = 0, 
                               dim = c(1, length(results_names)), 
                               dimnames = list(1, results_names)))
result_k <- results
k = 1


no_cores <- detectCores()
#cl <- makeCluster(no_cores, type="FORK") #Fork doesn't work in Windows
cl <- makeCluster(no_cores)  
registerDoParallel(cl)  

dir.create("output/sim_output_v1")

########################### LOOP ##########################################
no_nodes <- 8
node_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
node_i

#foreach(i = sp1[round(seq(1, 10780, 4))], .packages = c("tidyverse", "pracma"))  %dopar% {
#foreach(i = sp1[round(seq(2, 10780, 4))], .packages = c("tidyverse", "pracma"))  %dopar% {
#foreach(i = sp1[round(seq(3, 10780, 4))], .packages = c("tidyverse", "pracma"))  %dopar% {
foreach(i = sp1[round(seq(node_i, nrow(par_space), no_nodes))], .packages = c("tidyverse", "pracma"))  %dopar% {
    current_time <- Sys.time()
    A <- betaPDF(a = par_space$a[i], b = par_space$b[i], 
                 alpha = par_space$alpha[i], gamma = par_space$gamma[i],
                 k = 1)
    
    file <- paste("output/sim_output_v1/BetaSimulation_Results", "_spA_", i, ".csv", sep = "")
    
    write_csv(results[0, ], file)
    if(i > 1) {
        sp2 <- sp1
        sp2 <- sp2[-c(1:(i-1))]
    }
    for(j in sp2) {
        B <- betaPDF(a = par_space$a[j], b = par_space$b[j], 
                     alpha = par_space$alpha[j], gamma = par_space$gamma[j],
                     k = 1)
        
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
        
        result_k[1, 1:4] <- c(i, j, NicheDiss, RNO)
        result_k[1, 5:8] <- par_space[i, ]
        result_k[1, 9:12] <- par_space[j, ]
        
        write_csv(result_k, file, append = TRUE)
        
        #if(k %% 100 == 0) {
            #Elapsed_time <- Sys.time() - current_time
            #ETA <- (Elapsed_time/k)*length(sp2)
            
            #RTA <- ETA-Elapsed_time
            #print(paste("spA =", i, "| spB =", j, "Init_Time =", current_time, 
            #            "Elapsed_Time =", Elapsed_time, "ETA =", ETA, "RTA =", RTA))
            
        #}
        #k = k + 1
        
        
    }
}

stopCluster(cl)

