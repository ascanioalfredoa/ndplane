#### Load packages ####
library(tidyverse)
library(terra)
library(tidyterra)
library(RColorBrewer)
library(rnaturalearth) # boundary data
library(ggrepel)
#library(gridExtra)
#library(ggpubr)
#library(cowplot)

#### Read NDP data ####
rc <- read_rds("Results/NDP_Salamanders_ResponseCurves2.rds")
rc <- rbind(rc, read_rds("Results/NDP_Salamanders_ResponseCurves_soil_ag13.rds"))

ndp <- read_csv("Results/NDP_salamander_values2.csv")
ndp <- rbind(ndp, read_csv("Results/NDP_salamander_values_soil_ag13.csv")[, -1])

#### Add type of climatic variable to variables ####
rc$vartype <- ifelse(rc$var %in% paste("bio", 1:11, sep = ""), yes = "Temperature",
       no = ifelse(rc$var %in% paste("bio", 12:19, sep = ""), yes = "Precipitation",
                   no = "Soil"))
ndp$vartype <- ifelse(ndp$var %in% paste("bio", 1:11, sep = ""), yes = "Temperature",
                     no = ifelse(ndp$var %in% paste("bio", 12:19, sep = ""), yes = "Precipitation",
                                 no = "Soil"))

#### Read salamander points ####
sal <- vect("Data/Ambystoma/Amb_macmar_thin.shp")

#### Load basemap ####
basemap_sf <- ne_states(
    country = "united states of america",
    #scale = "medium", 
    returnclass = "sf") %>%
    #dplyr::select(adm0_sr) %>%
    filter(!(code_hasc %in% c("US.HI", "US.AK")))

#### Plot US map with both species ####
#Define color pallette
colpal <- c("#000000", "#62A39F")

p1 <- ggplot() + 
    geom_sf(data = basemap_sf, # grey fill behind raster
            fill = "grey90",
            col = NA) + 
    geom_spatvector(data = sal, aes(color = Cmmn_Nm), 
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

#### Plot NDP ###
p2 <- ndp %>%
    ggplot(aes(x = exclusivity, y = dissimilarity, label = var, fill = vartype)) +
    xlim(0, 1) + ylim(0, 1) +
    xlab("Niche Exclusivity") + ylab("Niche Dissimilarity") +
    geom_hline(yintercept = 0.5) +
    geom_vline(xintercept = 0.5) +
    geom_point(size = 4, shape = 21, alpha = 0.8) +
    geom_text_repel(size = 4) +
    scale_fill_manual("Variable type", 
                      values = c("#62A39F", "#000000", "#2F8745")) +
    theme_bw() +
    theme(#legend.position = "bottom",
        legend.position = c(0.87, 0.75),
        legend.background = element_rect(fill = "white", color = "black")) +
    guides(fill = guide_legend(override.aes = list(size=8)))

p2

#### Plot response curves ####
rc %>% 
    #filter(var %in% c("")) %>%
    #mutate(var = case_when(
        #var == "bio1" ~ "Annual Mean Temperature (°C)",
        #var == "bio2" ~ "Mean Diurnal Range (°C)",
        #var == "bio3" ~ "Isothermality (°C)",
        #var == "bio4" ~ "Temperature Seasonality (°C/100)",
        #var == "bio5" ~ "Max Temperature of Warmest Month (°C)",
        #var == "bio6" ~ "Min Temperature of Coldest Month (°C) [bio06]",
        #var == "bio7" ~ "Temperature Annual Range(°C)",
        #var == "bio8" ~ "Mean Temperature of Wettest Quarter (°C)",
        #var == "bio9" ~ "Mean Temperature of Warmest Quarter (°C)",
        #var == "bio10" ~ "Mean Temperature of Coldest Quarter (°C)",
        #var == "bio11" ~ "Mean Temperature of Coldest Quarter (°C)",
        #var == "bio12" ~ "Annual Precipitation (kg/m^2) [bio12]",
        #var == "bio13" ~ "Precipitation of Wettest Month (kg/m^2)",
        #var == "bio14" ~ "Precipitation of Driest Month (kg/m^2)",
        #var == "bio15" ~ "Precipitation Seasonality (kg/m^2)",
        #var == "bio16" ~ "Precipitation of Wettest Quarter (kg/m^2)",
        #var == "bio17" ~ "Precipitation of Driest Quarter (kg/m^2) [bio17]",
        #var == "bio18" ~ "Precipitation of Warmest Quarter (kg/m^2)",
        #var == "bio19" ~ "Precipitation of Coldest Quarter (kg/m^2)",
        #var == "aws.x" ~ "Available Water Storage (mL)",
        #var == "aws.y" ~ "Available Water Supply (mL)",
        #var == "awc_r" ~ "Available Water Capacity (mL) [awc_r]",
        #var == "ffd_r" ~ "Frost Free Days (Julian Days)",                                                                            
        #var =="resdept_r" ~ "Depth to Any Soil Restricted Layer (cm)",
        #var == "sandtotal_r" ~ "Sand Total (%)",
        #var == "claytotal_r" ~ "Clay Total (%) [claytotal_r]",
        #var == "silttotal_r" ~ "Silt Total (%)",
        #var == "om_r" ~ "Organic Matter (%)",
        #var == "soimoistdept_r" ~ "Depth to water table (cm) [soimoistdept_r]",
        #var == "dbovendry_r" ~ "Db Oven Dry (g/mL)",
        #var == "dbthirdbar_r" ~ "Db 0.33 Bar H2O (g/mL)"
    #)) %>%
    ggplot(aes(x = layer, y = value, color = source)) +
    geom_line(linewidth = 2, alpha = 0.7) +
    facet_wrap(vartype ~ var, scales = "free_x") +
    ylab("Frequency") +
    xlab("Environmental gradient") +
    scale_color_manual("Species", values = colpal) +
    theme_bw() +
    theme(legend.position = "top")
p3
