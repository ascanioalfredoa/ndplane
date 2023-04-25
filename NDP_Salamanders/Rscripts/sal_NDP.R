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

#Color pallette

#### Load spatialized occurrences ####
maps <- vect("Results/Amb_macopac_env.shp")
ggplot(maps) +
    geom_spatvector(aes(color = Cmmn_Nm))

colpal <- c("#000000", "#62A39F")

# Load basemap
basemap_sf <- ne_states(
    country = "united states of america",
    #scale = "medium", 
    returnclass = "sf") %>%
    #dplyr::select(adm0_sr) %>%
    filter(!(code_hasc %in% c("US.HI", "US.AK")))

# Plot everyting
p1 <- ggplot() + 
    geom_sf(data = basemap_sf, # grey fill behind raster
            fill = "grey90",
            col = NA) + 
    geom_spatvector(data = maps, aes(color = Cmmn_Nm), 
                    alpha = 0.7, size = 2) +
    geom_sf(data = basemap_sf, # lines on top of raster
            fill = NA,
            col = "grey20",
            lwd = 0.25) +
    scale_fill_gradientn(
        colors = colors(100),
        na.value = "transparent") +
    scale_color_manual("Species", values = colpal) +
    theme_bw() + 
    theme(legend.position = "bottom") +
    xlim(-96, -68) +
    labs(fill = NULL) +
    guides(color = guide_legend(override.aes = list(size=8)))
p1
table(maps$Cmmn_Nm)

#### Load Data ####
dat <- read_csv("Results/Amb_macopac_env.csv")
#dat <- dat %>% rename()
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

p2 <- res[["NDP"]] %>%
    mutate(comparison = paste(Taxa1, Taxa2, sep = "/"),
           vartype = c(rep("Temperature", 11), rep("Precipitation", 8),
                       rep("Soil", 12))
           #case_when(grepl("parent", Taxa1) & grepl("parent", Taxa2) ~ "parental",
           #         Taxa1 == "high_concolor" & Taxa2  ~ paste("hybrid", Taxa2, sep = "/"),
           #         Taxa2 == "hybrid" ~ paste("hybrid2", Taxa1, sep = "/")
           #)
    ) %>%
    #mutate(comparison = 
    #           case_when(grepl("parent", Taxa1) & grepl("parent", Taxa2) ~ "parental",
    #                     Taxa1 == "hybrid" ~ paste("hybrid", Taxa2, sep = "/"),
    #                     Taxa2 == "hybrid" ~ paste("hybrid2", Taxa1, sep = "/")
    #                     )
    #                     )%>%
    ggplot(aes(x = NicheEx, y = NicheDiss, label = Var, fill = vartype)) +
    xlab("Niche Exclusivity") + ylab("Niche Dissimilarity") +
    xlim(0, 1) + ylim(0, 1) +
    geom_hline(yintercept = 0.5) +
    geom_vline(xintercept = 0.5) +
    #geom_point(size = 2) + facet_wrap(~ Var) +
    geom_point(size = 4, 
               #fill = "#62A39F", 
               shape = 21, alpha = 0.8) +
    #geom_text(size = 5, nudge_x = 0.05, check_overlap = T) +
    geom_text_repel(size = 4) +
    scale_fill_manual("Variable type", 
                      values = c("#62A39F", "#000000", "#2F8745")) +
    theme_bw() +
    theme(#legend.position = "bottom",
          legend.position = c(0.87, 0.75),
          legend.background = element_rect(fill = "white", color = "black")) +
    guides(fill = guide_legend(override.aes = list(size=8)))
p2

p3 <-
    res[["Responses"]] %>% 
    mutate(vartype = ifelse(grepl("bio", var), "bioclimatic", "soil")) %>%
    filter(var %in% c("bio06", "bio12", "bio17", "soimoistdept_r", "awc_r", "claytotal_r")) %>%
    mutate(var = case_when(
        #var == "bio01" ~ "Annual Mean Temperature (°C)",
        #var == "bio02" ~ "Mean Diurnal Range (°C)",
        #var == "bio03" ~ "Isothermality (°C)",
        #var == "bio04" ~ "Temperature Seasonality (°C/100)",
        #var == "bio05" ~ "Max Temperature of Warmest Month (°C)",
        var == "bio06" ~ "Min Temperature of Coldest Month (°C) [bio06]",
        #var == "bio07" ~ "Temperature Annual Range(°C)",
        #var == "bio08" ~ "Mean Temperature of Wettest Quarter (°C)",
        #var == "bio09" ~ "Mean Temperature of Warmest Quarter (°C)",
        #var == "bio10" ~ "Mean Temperature of Coldest Quarter (°C)",
        #var == "bio11" ~ "Mean Temperature of Coldest Quarter (°C)",
        var == "bio12" ~ "Annual Precipitation (kg/m^2) [bio12]",
        #var == "bio13" ~ "Precipitation of Wettest Month (kg/m^2)",
        #var == "bio14" ~ "Precipitation of Driest Month (kg/m^2)",
        #var == "bio15" ~ "Precipitation Seasonality (kg/m^2)",
        #var == "bio16" ~ "Precipitation of Wettest Quarter (kg/m^2)",
        var == "bio17" ~ "Precipitation of Driest Quarter (kg/m^2) [bio17]",
        #var == "bio18" ~ "Precipitation of Warmest Quarter (kg/m^2)",
        #var == "bio19" ~ "Precipitation of Coldest Quarter (kg/m^2)",
        #var == "aws.x" ~ "Available Water Storage (mL)",
        #var == "aws.y" ~ "Available Water Supply (mL)",
        var == "awc_r" ~ "Available Water Capacity (mL) [awc_r]",
        #var == "ffd_r" ~ "Frost Free Days (Julian Days)",                                                                            
        #var =="resdept_r" ~ "Depth to Any Soil Restricted Layer (cm)",
        #var == "sandtotal_r" ~ "Sand Total (%)",
        var == "claytotal_r" ~ "Clay Total (%) [claytotal_r]",
        #var == "silttotal_r" ~ "Silt Total (%)",
        #var == "om_r" ~ "Organic Matter (%)",
        var == "soimoistdept_r" ~ "Depth to water table (cm) [soimoistdept_r]",
        #var == "dbovendry_r" ~ "Db Oven Dry (g/mL)",
        #var == "dbthirdbar_r" ~ "Db 0.33 Bar H2O (g/mL)"
    )) %>%
    ggplot(aes(x = x, y = y, color = category)) +
    geom_line(linewidth = 2, alpha = 0.7) +
    facet_wrap(vartype ~ var, scales = "free_x") +
    ylab("Frequency") +
    xlab("Environmental gradient") +
    scale_color_manual("Species", values = colpal) +
    theme_bw() +
    theme(legend.position = "top")
p3

gt <- arrangeGrob(p1 + theme(legend.text = element_text(size = 20)) +
                      guides(color = guide_legend(override.aes = list(size=8))), 
                  p2 + theme(axis.title = element_text(size = 20),
                             legend.text = element_text(size = 20)), 
                  p3 + theme(legend.position = "none",
                             axis.title = element_text(size = 20)),
                  ncol = 2, nrow = 2, 
                  layout_matrix = rbind(c(1, 2), 
                                        c(1, 3)))
# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
    draw_plot_label(label = c("A)", "B)", "C)"), size = 20,
                    x = c(-0.005, 0.49, 0.49), y = c(1, 1, 0.5)) # Add labels
p

pdf("Figures/Amacopac_NDP.pdf", width = 19, height = 11)
p
dev.off()

png("Figures/Amacopac_NDP.png", width = 11000, height = 6000, res = 600)
p
dev.off()
