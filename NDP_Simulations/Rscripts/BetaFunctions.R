#normfunction <- expression((1/(s*sqrt(2*pi)))*exp(-(1/2)*((x-mu)/s)^2))
betafunction_max <- expression(k*((x-a)^alpha)*((b-x)^gamma))

#Constructing Beta PDF - Normalized by peak = 1 for each species response
betaPDF <- function(a = 0, b = 1, alpha = 1, gamma = 1, k = 1) {
    x = round(seq(a, b, 0.001), 3)
    y = eval(betafunction_max)
    x = x[y >= 0]
    y = y[y >= 0]
    y = y/max(y)
    return(data.frame(x, y))
}

#Constructing Beta PDF - Normalized so the area under each curve is 1
betaPDF_A1 <- function(a = 0, b = 1, alpha = 1, gamma = 1, k = 1) {
    x = round(seq(a, b, 0.001), 3)
    y = eval(betafunction_max)
    x = x[y >= 0]
    y = y[y >= 0]
    #y = y/max(y)
    return(data.frame(x, y))
}

#General betaPDF for a pair of species
betaPDF_pair <- function(a1 = 0, b1 = 1, alpha1 = 1, gamma1 = 1, 
                    a2 = 0, b2 = 1, alpha2 = 1, gamma2 = 1,
                    k = 1, pden = 3000) {
    
    betafunction_max1 <- expression(k*((x-a1)^alpha1)*((b1-x)^gamma1))
    betafunction_max2 <- expression(k*((x-a2)^alpha2)*((b2-x)^gamma2))
    
    x = round(seq(min(a1, a2), max(b1, b2), length = pden), 3)
    y1 = eval(betafunction_max1)
    y2 = eval(betafunction_max2)
    
    x = x[y1 >= 0 | y2 >= 0]
    #x1 = x[y1 >= 0]
    y1[x < a1 | x > b1] <- 0
    y1 = y1/max(y1)
    
    
    #x2 = x[y2 >= 0]
    #y2 = y2[y2 >= 0]
    y2[x < a2 | x > b2] <- 0
    y2 = y2/max(y2)
    return(data.frame(x, y1, y2))
}