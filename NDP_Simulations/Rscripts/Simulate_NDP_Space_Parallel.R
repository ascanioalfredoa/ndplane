################### Niche Divergence Plane Simulations #########################

#### Generate parameter space ####
library(tidyverse)
library(pracma)
library(doParallel)
source("Rscripts/BetaFunctions.R")

# Set parameter domains 
a = b <- round(seq(0, 1, 0.1), 3)
alpha = gamma <- round(c(0.001, 0.01, seq(0.1, 1, 0.2), seq(2, 10, 2), 100, 1000), 3)
#alpha = gamma <- seq(0.1, 10, 0.5)

# Set parameter space (combination of parameters for any beta function)
par_space <- expand.grid(a = a, b = b, alpha = alpha, gamma = gamma)
par_space <- par_space[par_space$a < par_space$b, ]
gc()

dir.create("Results")
write_csv(par_space, "Results/Beta_Parameter_Space.csv")

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
#result <- foreach(i=10:10000)
#stopCluster(cl)  

dir.create("Results/sim_output_v2")

foreach(i = sp1, .packages = c("tidyverse", "pracma"))  %dopar% {
    current_time <- Sys.time()
    A <- betaPDF(a = par_space$a[i], b = par_space$b[i], 
                 alpha = par_space$alpha[i], gamma = par_space$gamma[i],
                 k = 1)
    
    file <- paste("Results/sim_output_v2/BetaSimulation_Results", "_spA_", i, ".csv", sep = "")
    
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
        
        RNO <- 1 - (min(c(max(A$x), max(B$x))) - max(c(min(A$x), min(B$x))))/(max(c(max(A$x), max(B$x))) - min(c(min(A$x), min(B$x))))
        if(RNO > 1) RNO <- 1
        
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

