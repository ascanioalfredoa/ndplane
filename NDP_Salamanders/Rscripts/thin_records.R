thin_records <- function(data, lonlat = c("lon", "lat"), thin.par, reps = 10) {
    
    if(any(class(data) == "data.frame") | any(class(data) == "tibble")) {
        spatdat <- data %>% select(lonlat[1], lonlat[2]) %>% as.matrix()
        distmat <- terra::distance(spatdat, lonlat = T) %>% 
            as.matrix()
        
    } else if(class(data) == "matrix") {
        spatdat <- data %>% select(lonlat[1], lonlat[2])
        distmat <- terra::distance(spatdat, lonlat = T, unit = "km") %>% 
            as.matrix()
        
    } else if(class(data) == "SpatVector"){
        spatdat <- data
        distmat <- terra::distance(spatdat, unit = "km") %>% 
            as.matrix()
    }
    
    res <- 0
    
    for(i in 1:reps) {
        distlogi <- distmat < thin.par
        diag(distlogi) <- F
        sumvec <- rowSums(distlogi)
        keepvec <- sumvec >= 0
        
        
        while (any(distlogi) && sum(keepvec) > 1) {
            RemoveRec <- which(sumvec == max(sumvec))
            
            if (length(RemoveRec) > 1 && any(distlogi[RemoveRec, RemoveRec])) {
                RemoveRec <- sample(RemoveRec, 1)
            }
            
            distlogi[RemoveRec, ] <- FALSE
            distlogi[, RemoveRec] <- FALSE
            keepvec[RemoveRec] <- FALSE
            
            sumvec <- rowSums(distlogi)
            
        }
        
        if(sum(keepvec) > sum(res)) {
            res <- keepvec
        } else if(i > 0 && sum(keepvec) > sum(res)) {
            res <- list(res, keepvec)[[
                which.max(c(mean(distmat[res, res]),
                            mean(distmat[res, res])))
            ]]
        }
    }
    data[res, ]
}
