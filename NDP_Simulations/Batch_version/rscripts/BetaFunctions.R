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

