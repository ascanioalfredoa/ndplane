#### Load packages ####
library(tidyverse)
library(viridis)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)
library(cowplot)
#### Read data ####
vsp_sims <- read_csv("output/vsp_sims_results.csv")

#### Check range of D and I ####
range(vsp_sims$D)
range(vsp_sims$I)

round(seq(0.15, 1, 0.25), 2) %in% round(vsp_sims$D, 3) # Check that intervals appear in vsp D

round(seq(0.40, 1, 0.15), 2) %in% round(vsp_sims$I, 3) # Check that intervals appear in vsp I

vsp_sims <- 
    vsp_sims %>%
    mutate(ND_Magnitude = sqrt((Exclusivity^2) + (Dissimilarity^2)),
           ND_angle = round(atan2(Exclusivity, Dissimilarity)*(180/pi), 2))

#### Filter and plot D intervals ####
vsp_sims_D <- vsp_sims %>%
    filter(D %in% round(seq(0.15, 1, 0.25), 2))

png(filename = "Figures/vsp_D_NDP.png", width = 5000, height = 4500, res = 600)
gg_a <- 
    vsp_sims_D %>%
    ggplot(aes(x = Exclusivity, y = Dissimilarity, color = as.factor(D))) +
    geom_point(alpha = 0.5) +
    #scale_color_viridis("Similarity D", discrete = TRUE) +
    scale_color_manual("Similarity D", values = brewer.pal(4, "Set1")) +
    theme_bw() +
    theme(legend.position = "bottom")
gg_a
dev.off()


#### Filter and plot I intervals ####

vsp_sims_I <- vsp_sims %>%
    filter(I %in% round(seq(0.40, 0.95, 0.15), 2))

png(filename = "Figures/vsp_I_NDP.png", width = 5000, height = 4500, res = 600)
gg_b <- 
    vsp_sims_I %>%
    ggplot(aes(x = Exclusivity, y = Dissimilarity, color = as.factor(I))) +
    geom_point(alpha = 0.5) +
    #scale_color_viridis("Similarity I", discrete = TRUE) +
    scale_color_manual("Similarity I", values = brewer.pal(4, "Set1")) +
    theme_bw() +
    theme(legend.position = "bottom")
gg_b
dev.off()

png(filename = "Figures/vsp_D_I.png", width = 5000, height = 4500, res = 600)
gg_c <- 
    vsp_sims %>%
    select(I, D) %>%
    distinct() %>%
    ggplot(aes(x = D, y = I)) +
    geom_point() +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw()
gg_c
dev.off()

png(filename = "Figures/vsp_NDMag_I.png", width = 5000, height = 4500, res = 600)
gg_d <- 
    vsp_sims %>%
    select(I, ND_Magnitude) %>%
    distinct() %>%
    ggplot(aes(x = I, y = ND_Magnitude)) +
    geom_point() +
    xlim(0, 1) + ylim(0, sqrt(2)) +
    theme_bw()
gg_d
dev.off()

png(filename = "Figures/vsp_NDMag_D.png", width = 5000, height = 4500, res = 600)
gg_e <- 
    vsp_sims %>%
    select(D, ND_Magnitude) %>%
    distinct() %>%
    ggplot(aes(x = D, y = ND_Magnitude)) +
    geom_point() +
    xlim(0, 1) + ylim(0, sqrt(2)) +
    theme_bw()
gg_e
dev.off()


png(filename = "Figures/vsp_Dissimilarity_D.png", width = 5000, height = 4500, res = 600)
gg_f <- 
    vsp_sims %>%
    select(D, Dissimilarity) %>%
    distinct() %>%
    ggplot(aes(x = D, y = Dissimilarity)) +
    geom_point() +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw()
gg_f
dev.off()

png(filename = "Figures/vsp_Dissimilarity_I.png", width = 5000, height = 4500, res = 600)
gg_g <- 
    vsp_sims %>%
    select(I, Dissimilarity) %>%
    distinct() %>%
    ggplot(aes(x = I, y = Dissimilarity)) +
    geom_point() +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw()
gg_g
dev.off()

png(filename = "Figures/vsp_Exclusivity_D.png", width = 5000, height = 4500, res = 600)
gg_h <- 
    vsp_sims %>%
    select(D, Exclusivity) %>%
    distinct() %>%
    ggplot(aes(x = D, y = Exclusivity)) +
    geom_point() +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw()
gg_h
dev.off()

png(filename = "Figures/vsp_Exclusivity_I.png", width = 5000, height = 4500, res = 600)
gg_i <- 
    vsp_sims %>%
    select(I, Exclusivity) %>%
    distinct() %>%
    ggplot(aes(x = I, y = Exclusivity)) +
    geom_point() +
    xlim(0, 1) + ylim(0, 1) +
    theme_bw()
gg_i
dev.off()

gg_ab <- marrangeGrob(list(gg_a, gg_b), ncol = 2, nrow = 1)
gg_cde <- marrangeGrob(list(gg_c, gg_d, gg_e), ncol = 3, nrow = 1)
marrangeGrob(list(gg_ab, gg_cde), nrow = 2, ncol = 1)


gg_a <- 
    gg_a +
    theme(text = element_text(size = 20)) +
    scale_color_manual("Schoener's D",
                       values = c("#000000", "#62A39F", "#725745", "#0b40ab")) +
    guides(color = guide_legend(override.aes = list(size=8, alpha=1)))

gg_b <- 
    gg_b +
    theme(text = element_text(size = 20)) +
    scale_color_manual("Hellinger's I",
                       values = c("#000000", "#62A39F", "#725745", "#0b40ab")) +
    guides(color = guide_legend(override.aes = list(size=8, alpha=1)))

gg_c <- 
    gg_c + theme(text = element_text(size = 20)) +
    xlab("Schoener's D") + ylab("Hellinger's I")


gg_d <- 
    gg_d + theme(text = element_text(size = 20)) +
    ylab("Niche Divergence Magnitude") + xlab("Hellinger's I")


gg_e <- 
    gg_e + theme(text = element_text(size = 20)) +
    ylab("Niche Divergence Magnitude") + xlab("Schoener's D")


gg_f <- 
    gg_f + theme(text = element_text(size = 20)) +
    ylab("Niche Dissimilarity") + xlab("Schoener's D")

gg_g <- 
    gg_g + theme(text = element_text(size = 20)) +
    ylab("Niche Dissimilarity") + xlab("Hellinger's I")

gg_h <- 
    gg_h + theme(text = element_text(size = 20)) +
    ylab("Niche Exclusivity") + xlab("Schoener's D")

gg_i <- 
    gg_i + theme(text = element_text(size = 20)) +
    ylab("Niche Exclusivity") + xlab("Hellinger's I")



gt <- arrangeGrob(gg_a, gg_b, gg_c, gg_d, gg_e, 
                  layout_matrix = rbind(
                      c(1, 1, 1, 2, 2, 2),
                      c(3, 3, 4, 4, 5, 5)
                  ))
# Add labels to the arranged plots
p <- as_ggplot(gt) +                                # transform to a ggplot
    draw_plot_label(label = c("A)", "B)", "C)", "D)", "E)"), size = 20,
                    x = c(0, 0.5, 0, 0.33, 0.66), y = c(1, 1, 0.5, 0.5, 0.5)) # Add labels

png(filename = "Figures/Fig4.png", width = 11000, height = 9000, res = 600)
tiff(filename = "Figures/Fig4.tif", width = 11000, height = 9000, res = 600)
p + ggtitle("Figure 4")
#grid.arrange(gg_a, gg_b, gg_c, gg_d, gg_e, 
#             layout_matrix = rbind(
#                 c(1, 1, 1, 2, 2, 2),
#                 c(3, 3, 4, 4, 5, 5)
#             ))
dev.off()

png(filename = "Figures/FigSX_vsp_cors.png", width = 10000, height = 9000, res = 600)
grid.arrange(gg_f, gg_g, gg_h, gg_i, 
             layout_matrix = rbind(
                 c(1, 1, 1, 2, 2, 2),
                 c(3, 3, 3, 4, 4, 4)
             ))
dev.off()

vsp_sims %>%
    select(I, Exclusivity) %>%
    #distinct() %>%
    cor.test(~ I + Exclusivity, data = .)

vsp_sims %>%
    select(D, Exclusivity) %>%
    #distinct() %>%
    cor.test(~ D + Exclusivity, data = .)

vsp_sims %>%
    select(I, Dissimilarity) %>%
    #distinct() %>%
    cor.test(~ I + Dissimilarity, data = .)

vsp_sims %>%
    select(D, Dissimilarity) %>%
    #distinct() %>%
    cor.test(~ D + Dissimilarity, data = .)
