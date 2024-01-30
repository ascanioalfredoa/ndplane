#### Load packages ####
library(tidyverse)
library(terra)
library(tidyterra)
library(RColorBrewer)
library(rnaturalearth) # boundary data
library(ggrepel)
library(gridExtra)
library(ggpubr)
library(cowplot)

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
    ggplot(aes(x = exclusivity, y = dissimilarity, label = var, fill = vartype, 
               shape = vartype)) +
    xlim(0, 1) + ylim(0, 1) +
    xlab("Niche Exclusivity") + ylab("Niche Dissimilarity") +
    geom_hline(yintercept = 0.5) +
    geom_vline(xintercept = 0.5) +
    geom_point(size = 6, alpha = 0.8) +
    geom_text_repel(size = 8) +
    scale_fill_manual("Variable type", 
                      values = c("#62A39F", "#000000", "#2F8745")) +
    theme_bw() +
    scale_shape_manual("Variable type", 
                      values = c(21, 22, 23)) +
    theme_bw() +
    theme(#legend.position = "bottom",
        legend.position = c(0.87, 0.75),
        legend.background = element_rect(fill = "white", color = "black")) +
    guides(fill = guide_legend(override.aes = list(size=8)))

p2

#### Plot response curves ####
p3 <- 
    rc %>% 
    filter(var %in% c("bio1", "bio3", "bio12", "bio14", "Sand", "WC3rdbar")) %>%
    mutate(var = case_when(
        var == "bio1" ~ "Annual Mean Temperature (°C) [bio1]",
        #var == "bio2" ~ "Mean Diurnal Range (°C) [bio2]",
        var == "bio3" ~ "Isothermality (°C) [bio3]",
        #var == "bio4" ~ "Temperature Seasonality (°C/100) [bio4]",
        #var == "bio5" ~ "Max Temperature of Warmest Month (°C) [bio5]",
        #var == "bio6" ~ "Min Temperature of Coldest Month (°C) [bio6]",
        #var == "bio7" ~ "Temperature Annual Range(°C) [bio7]",
        #var == "bio8" ~ "Mean Temperature of Wettest Quarter (°C) [bio8]",
        #var == "bio9" ~ "Mean Temperature of Warmest Quarter (°C) [bio9]",
        #var == "bio10" ~ "Mean Temperature of Coldest Quarter (°C) [bio10]",
        #var == "bio11" ~ "Mean Temperature of Coldest Quarter (°C) [bio11]",
        var == "bio12" ~ "Annual Precipitation (kg/m^2) [bio12] [bio12]",
        #var == "bio13" ~ "Precipitation of Wettest Month (kg/m^2) [bio13]",
        var == "bio14" ~ "Precipitation of Driest Month (kg/m^2) [bio14]",
        #var == "bio15" ~ "Precipitation Seasonality (kg/m^2) [bio15]",
        #var == "bio16" ~ "Precipitation of Wettest Quarter (kg/m^2) [bio16]",
        #var == "bio17" ~ "Precipitation of Driest Quarter (kg/m^2) [bio17]",
        #var == "bio18" ~ "Precipitation of Warmest Quarter (kg/m^2) [bio18]",
        #var == "bio19" ~ "Precipitation of Coldest Quarter (kg/m^2) [bio19]",
        #var == "AWC" ~ "Available Water Capacity (mL) [AWC]",
        #var == "AWS" ~ "Available Water Storage (mL) [AWS]",
        #var == "CaC03" ~ "Calcium Carbonate Equivalent [CaCO3]",                                                                            
        #var == "Clay" ~ "Clay Total (%)",
        #var == "Db3rdbar" ~ "Bulk Density at 1/3 bar (g/mL) [Db3rdbar]",
        #var == "Dep2AnyRes" ~ "Depth to Any Soil Restricted Layer (cm) [Dep2AnyRes]",
        #var == "Dep2WatTbl" ~ "Depth to water table (cm) [Dep2WatTbl]",
        #var == "FrostFDays" ~ "Frost Free Days (Julian Days) [FrostFDays]",                                                                            
        #var == "LiqLim" ~ "Liquid Limit [LiqLim]",
        #var == "OrgMatter" ~ "Organic Matter (%)",
        #var == "pHwater" ~ "pH",
        var == "Sand" ~ "Sand Total (%)",
        #var == "Silt" ~ "Silt Total (%)",
        #var == "WC15Bar" ~ "Water Content at 15 Bar (%) [WC15Bar]",
        var == "WC3rdbar" ~ "Water Content at 1/3 Bar (%) [WC3rdbar]"
    )) %>%
    ggplot(aes(x = layer, y = value, color = source)) +
    geom_line(linewidth = 2, alpha = 0.7) +
    facet_wrap(vartype ~ var, scales = "free_x") +
    ylab("Frequency") +
    xlab("Environmental gradient") +
    scale_color_manual("Species", values = colpal) +
    theme_bw() +
    theme(legend.position = "top")
p3

gt <- arrangeGrob(p1 + theme(legend.text = element_text(size = 14)) +
                      guides(color = guide_legend(override.aes = list(size=3))), 
                  p2 + theme(axis.title = element_text(size = 20),
                             legend.text = element_text(size = 20),
                             text = element_text(size = 20)), 
                  p3 + theme(text = element_text(size = 16),
                             legend.position = "none",
                             axis.title = element_text(size = 20)),
                  ncol = 2, nrow = 2, 
                  layout_matrix = rbind(c(1, 1, 3, 3, 3, 3, 3),
                                        c(1, 1, 3, 3, 3, 3, 3),
                                        c(1, 1, 3, 3, 3, 3, 3),
                                        c(2, 2, 2, 2, 2, 2, 2),
                                        c(2, 2, 2, 2, 2, 2, 2),
                                        c(2, 2, 2, 2, 2, 2, 2),
                                        c(2, 2, 2, 2, 2, 2, 2),
                                        c(2, 2, 2, 2, 2, 2, 2)))

# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
    draw_plot_label(label = c("A)", "B)", "C)"), size = 20,
                    x = c(-0.005, 0.30, -0.005), y = c(1, 1, 0.62)) # Add labels
#x = c(-0.005, 0.49, 0.49), y = c(1, 1, 0.5)) # Add labels
#x = c(-0.005, -0.005, -0.005), y = c(1, 2, 3)) # Add labels
p

pdf("Figures/Amacopac_NDP.pdf", width = 19, height = 15)
p
dev.off()

png("Figures/Amacopac_NDP.png", width = 11000, height = 9000, res = 600)
p
dev.off()


########### Supplementary figure with all response curves ###############
p4 <- 
    rc %>% 
    #filter(var %in% c("")) %>%
    mutate(var = case_when(
        var == "bio1" ~ "Annual Mean Temperature (°C) [bio1]",
        var == "bio2" ~ "Mean Diurnal Range (°C) [bio2]",
        var == "bio3" ~ "Isothermality (°C) [bio3]",
        var == "bio4" ~ "Temperature Seasonality (°C/100) [bio4]",
        var == "bio5" ~ "Max Temperature of Warmest Month (°C) [bio5]",
        var == "bio6" ~ "Min Temperature of Coldest Month (°C) [bio6]",
        var == "bio7" ~ "Temperature Annual Range(°C) [bio7]",
        var == "bio8" ~ "Mean Temperature of Wettest Quarter (°C) [bio8]",
        var == "bio9" ~ "Mean Temperature of Warmest Quarter (°C) [bio9]",
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
        var == "AWC" ~ "Available Water Capacity (mL) [AWC]",
        var == "AWS" ~ "Available Water Storage (mL) [AWS]",
        var == "CaC03" ~ "Calcium Carbonate Equivalent [CaCO3]",                                                                            
        var == "Clay" ~ "Clay Total (%)",
        var == "Db3rdbar" ~ "Bulk Density at 1/3 bar (g/mL) [Db3rdbar]",
        var == "Dep2AnyRes" ~ "Depth to Any Soil Restricted Layer (cm) [Dep2AnyRes]",
        var == "Dep2WatTbl" ~ "Depth to water table (cm) [Dep2WatTbl]",
        var == "FrostFDays" ~ "Frost Free Days (Julian Days) [FrostFDays]",                                                                            
        var == "LiqLim" ~ "Liquid Limit [LiqLim]",
        var == "OrgMatter" ~ "Organic Matter (%)",
        var == "pHwater" ~ "pH",
        var == "Sand" ~ "Sand Total (%)",
        var == "Silt" ~ "Silt Total (%)",
        var == "WC15Bar" ~ "Water Content at 15 Bar (%) [WC15Bar]",
        var == "WC3rdbar" ~ "Water Content at 1/3 Bar (%) [WC3rdbar]"
    )) %>%
    ggplot(aes(x = layer, y = value, color = source)) +
    geom_line(linewidth = 2, alpha = 0.7) +
    facet_wrap(vartype ~ var, scales = "free_x") +
    ylab("Frequency") +
    xlab("Environmental gradient") +
    scale_color_manual("Species", values = colpal) +
    theme_bw() +
    theme(legend.position = "top")
p4


pdf("Figures/Amacopac_AllResponses.pdf", width = 19, height = 15)
p4
dev.off()

png("Figures/Amacopac_AllResponses.png", width = 11000, height = 9000, res = 600)
p4
dev.off()


#### Correlation between empirical Dissimilarity/Exclusivity and D/I ####
library(GGally)
pairs(ndp[, c(1:2, 4:5)])

ggpairs(ndp, 
        columns = c(1, 2, 4, 5), 
        aes(label = var, fill = vartype, shape = vartype),
        upper = list(continuous = wrap("cor", size = 18)), 
        lower = list(continuous = wrap("points", size = 3)), 
        legends = T) +
    scale_fill_manual("Variable type", values = c("#62A39F", "#000000", "#2F8745")) +
    scale_shape_manual("Variable type", values = c(21, 22, 23)) +
    theme_bw() +
    theme(text = element_text(size = 18),
          legend.position = "bottom")

cor(ndp[, c(1:2, 4:5)])
cor.test(ndp$dissimilarity, ndp$D)
cor.test(ndp$dissimilarity, ndp$I)
cor.test(ndp$exclusivity, ndp$D)
cor.test(ndp$exclusivity, ndp$I)
