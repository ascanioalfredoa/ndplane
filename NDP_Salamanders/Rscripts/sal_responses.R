#### Load packages ####
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(pracma)
library(gridExtra)
library(terra)
library(tidyterra)
library(rnaturalearth) # boundary data
library(ggpubr)
library(cowplot)

#### Load Data ####
dat <- read_csv("Results/Amb_macopac_env.csv")

#### NDP for rattlesnakes ####

# Select variable names for which we want to compare the niche divergence
var <- names(dat)
var <- var[-c(1:2)]
var <- var[!grepl("cl$|_Levels", var)]

#### Loop Structure ####
category <- "Cmmn_Nm"

res <- list()
res[["NDP"]] <- tibble(Taxa1 = NULL, Taxa2 = NULL, Var = NULL, NicheDiss = NULL, NicheEx = NULL)
res[["Responses"]] <- tibble(category = NULL, x = NULL, y = NULL)

for(i in 1:length(var)) {
    # Select out a single variable at a time and calculate overall density
    # this standardizes a number of points along the X axis for all species of interests
    envdat <- dat[!is.na(dat[[var[i]]]), ]
    
    var_dens <- density(envdat[[var[i]]], n = 3000) # setting 3000 points for every variable
    x_all <- density(envdat[[var[i]]], n = 3000)$x # Saves the points over x for which we want the density
    
    # Apply density calculation by desired groups
    
    cat_names <- envdat %>% 
        #filter(Source == dat_source) %>% 
        select(category) %>% 
        unique() %>%
        pull()
    
    ndp_data <- 
        envdat %>%
        #filter(Source == dat_source) %>%
        select(category, var[i]) %>% #Selects variable of interest
        group_by_at(category) %>% #Groups by categorization of interest
        group_map(~ .x %>% 
                      pull(1) %>% 
                      density(n = 3000)
        ) %>% # group_map executes the density function for each category
        set_names(nm = cat_names) %>%
        lapply(function(x) data.frame(x = x$x,  y = x$y)) %>% # Extracts only x and y values for each density curve %>% 
        lapply(function(x) approx(x$x, x$y, xout = x_all) %>%
                   do.call(cbind, .) %>%
                   as.data.frame() %>%
                   filter(complete.cases(.)) %>%
                   mutate(Cumden = cumsum(y)) %>%
                   mutate(Cumden = Cumden/max(Cumden)) %>%
                   filter(Cumden >= 0.025 & Cumden <= 0.975) %>%
                   mutate(y = y/max(y)) %>%
                   select(1:2)
        ) %>%  # Calculates density for each species on standard points within variable range (x_all)
        do.call(rbind, .) %>%
        rownames_to_column(var = "category") %>%
        mutate(category = str_remove(category, "\\.[0-9]*"))# %>%
    
    #ndp_data %>%
    #    ggplot(aes(x = x, y = y, color = category, fill = category)) +
    #    geom_line(size = 2)
    
    taxa_comb <- t(combn(unique(ndp_data$category), 2))
    
    res[["Responses"]] <- bind_rows(res[["Responses"]],
                                    ndp_data %>% mutate(var = var[[i]]))
    
    for(j in 1:nrow(taxa_comb)) {
        A <- ndp_data %>% filter(grepl(taxa_comb[j, 1], category))
        B <- ndp_data %>% filter(grepl(taxa_comb[j, 2], category))
        
        #### Integration by determining the lower function ####
        # Where the values of X are shared between distribution
        # And where the values of Y are the lowest between the two possibilities
        
        C <- 
            A[A$x %in% B$x, ] %>%
            left_join(B[B$x %in% A$x, ], by = "x", suffix = c("_a", "_b")) %>%
            #rename(!!(.[["Source_a"]] %>% unique()) := y_a, !!(.[["Source_b"]] %>% unique()) := y_b) %>%
            dplyr::select(2, 3, 5) %>% 
            pivot_longer(cols = 2:3, names_to = "source_curve", values_to = "y") %>%
            group_by(x) %>%
            summarise(y = min(y))
        
        # Calculate area
        NicheDiss <- 1 - (trapz(C$x, C$y)/trapz(A$x, A$y) + trapz(C$x, C$y)/trapz(B$x, B$y))/2
        
        RNO <- 1 - (min(c(max(A$x), max(B$x))) - max(c(min(A$x), min(B$x))))/(max(c(max(A$x), max(B$x))) - min(c(min(A$x), min(B$x))))
        if(RNO > 1) RNO <- 1
        
        res[["NDP"]] <- 
            bind_rows(res[["NDP"]],
                      tibble(Taxa1 = taxa_comb[j, 1], Taxa2 = taxa_comb[j, 2], 
                             Var = var[i],
                             NicheDiss = NicheDiss, NicheEx = RNO)
            )
    }
}

p3 <-
    res[["Responses"]] %>% 
    mutate(vartype = ifelse(grepl("bio", var), "bioclimatic", "soil")) %>%
    #filter(var %in% c("bio06", "bio12", "bio17", "soimoistdept_r", "awc_r", "claytotal_r")) %>%
    mutate(var = case_when(
        var == "bio01" ~ "Annual Mean Temperature (°C) [bio01]",
        var == "bio02" ~ "Mean Diurnal Range (°C) [bio02]",
        var == "bio03" ~ "Isothermality (°C) [bio03]",
        var == "bio04" ~ "Temperature Seasonality (°C/100) [bio04]",
        var == "bio05" ~ "Max Temperature of Warmest Month (°C) [bio05]",
        var == "bio06" ~ "Min Temperature of Coldest Month (°C) [bio06]",
        var == "bio07" ~ "Temperature Annual Range(°C) [bio07]",
        var == "bio08" ~ "Mean Temperature of Wettest Quarter (°C) [bio08]",
        var == "bio09" ~ "Mean Temperature of Warmest Quarter (°C) [bio09]",
        var == "bio10" ~ "Mean Temperature of Coldest Quarter (°C) [bio10]",
        var == "bio11" ~ "Mean Temperature of Coldest Quarter (°C) [bio11]",
        var == "bio12" ~ "Annual Precipitation (kg/m^2) [bio12] [bio12]",
        var == "bio13" ~ "Precipitation of Wettest Month (kg/m^2) [bio13]",
        var == "bio14" ~ "Precipitation of Driest Month (kg/m^2) [bio14]",
        var == "bio15" ~ "Precipitation Seasonality (kg/m^2) [bio15]",
        var == "bio16" ~ "Precipitation of Wettest Quarter (kg/m^2) [bio16]",
        var == "bio17" ~ "Precipitation of Driest Quarter (kg/m^2) [bio17]",
        var == "bio18" ~ "Precipitation of Warmest Quarter (kg/m^2) [bio18]",
        var == "bio19" ~ "Precipitation of Coldest Quarter (kg/m^2) [bio19]",
        var == "aws.x" ~ "Available Water Storage (mL) [aws.x]",
        var == "aws.y" ~ "Available Water Supply (mL) [aws.y]",
        var == "awc_r" ~ "Available Water Capacity (mL) [awc_r]",
        var == "ffd_r" ~ "Frost Free Days (Julian Days) [ffd_r]",                                                                            
        var =="resdept_r" ~ "Depth to Any Soil Restricted Layer (cm) [resdept_r]",
        var == "sandtotal_r" ~ "Sand Total (%) [sandtotal_r]",
        var == "claytotal_r" ~ "Clay Total (%) [claytotal_r]",
        var == "silttotal_r" ~ "Silt Total (%) [silttotal_r]",
        var == "om_r" ~ "Organic Matter (%) [om_r]",
        var == "soimoistdept_r" ~ "Depth to water table (cm) [soimoistdept_r]",
        var == "dbovendry_r" ~ "Db Oven Dry (g/mL) [dbovendry_r]",
        var == "dbthirdbar_r" ~ "Db 0.33 Bar H2O (g/mL) [dbovendry_r]"
    )) %>%
    ggplot(aes(x = x, y = y, color = category)) +
    geom_line(linewidth = 2, alpha = 0.7) +
    facet_wrap(vartype ~ var, scales = "free_x", ncol = 4) +
    ylab("Frequency") +
    xlab("Environmental gradient") +
    scale_color_manual("Species", values = colpal) +
    theme_bw() +
    theme(legend.position = "bottom") +
    theme(legend.text = element_text(size = 20)) +
    guides(color = guide_legend(override.aes = list(size=8)))


pdf("Figures/Amacopac_AllResponses.pdf", width = 18, height = 16)
p3
dev.off()

png("Figures/Amacopac_AllResponses.png", width = 10000, height = 8500, res = 600)
p3
dev.off()
