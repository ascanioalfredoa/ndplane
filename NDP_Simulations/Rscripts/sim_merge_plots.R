#### Reading files and plotting ####
library(tidyverse)
library(RColorBrewer)
sim_files <- list.files("sim_output_v1/", full.names = TRUE)


#### Merging all tables into 1 ####
sim <- read.csv(sim_files[1])
write.table(sim, "full_sim_output_v1.csv", row.names = F, sep = ",")

for(i in 2:length(sim_files)) {
  sim <- read.csv(sim_files[i])
  write.table(sim, "full_sim_output_v1.csv", append = T, row.names = F, sep = ",", col.names = F)
  if(i %% 100 == 0) print(i)
}

#### load full data ####
simdata <- read_csv("Results/full_sim_output_v1.csv")
Sys.Date()

png("Figures/sim_plot_all.png", width = 1000, height = 1000)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), 
     xlab = "Exclusivity", ylab = "Dissimilarity",
     pch = 19, cex = 0.2)
for(i in 1:length(sim_files)) {
  sim <- read.csv(sim_files[i])
  sim <- sim[complete.cases(sim), ]
  if(nrow(sim) == 0) next
  
  #mean(c(sim$a_A[1], sim$b_A[1]))
  
  points(x = sim$Exclusivity, y = sim$Dissimilarity, pch = 19, cex = 0.2)
  if(i %% 100 == 0) print(i)
}
dev.off()

sum(complete.cases(simdata))/nrow(simdata) # 82.56% of data have values for exclusivity and dissimilarity
sum(is.na(simdata$Dissimilarity))/nrow(simdata) #17.44% of NA vales were NAs in the dissimilarity calculation
sum(is.na(simdata$Exclusivity)) #All rows of the simulation have exclusivity

simdata2 <- simdata[complete.cases(simdata), ] # Removing all rows for which dissimilarity was NA

sum(simdata2$Dissimilarity < 0)
sum(simdata2$Exclusivity < 0)

test <- simdata2[sample(1000, ), ]

library(tidyverse)

#### Checking negative values of dissimilarity ####

negativediss <- 
  simdata2 %>%
  filter(Dissimilarity < 0)

nrow(negativediss)/nrow(simdata2) # 0.1% of values yielded negative

#### Checking NA values of dissimilarity ####
nadiss <- 
  simdata %>%
  filter(is.na(Dissimilarity))

nrow(nadiss)/nrow(simdata) #17.44% of all results yielded NA for dissimilarity. Why?

#View(nadiss[sample(1:nrow(nadiss), 1000), ])

table(nadiss$alpha_A, nadiss$gamma_A)
table(nadiss$alpha_B, nadiss$gamma_B)
table(nadiss$alpha_A, nadiss$alpha_B)
table(nadiss$gamma_A, nadiss$gamma_B)

na_table <- table(nadiss$alpha_A, nadiss$gamma_A, nadiss$alpha_B, nadiss$gamma_B)

str(na_table)
tail(na_table[, , , 13])

#### a) Conservatism to fully Nested case ####
nested_case1 <- 
  simdata2 %>%
  mutate(mid_A = (a_A + b_A)/2, mid_B = (a_B + b_B)/2) %>%
  filter(mid_A == mid_B & alpha_A == gamma_A & alpha_B == gamma_B &
           alpha_A >= 1 & alpha_B >= 1 & gamma_A >= 1 & gamma_B >= 1)

100*nrow(nested_case1)/nrow(simdata2) # Nested Case 1 represents 0.006% of all simulations

png("sim_plot_nested1.png", width = 2000, height = 1000, res = 300)
nested_case1 %>%
  ggplot(aes(x = Exclusivity, y = Dissimilarity), col = "red") +
  geom_point() +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw() +
  geom_vline(aes(xintercept = 0.5), color = "red") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_hline(aes(yintercept = 0.5), color = "red")
dev.off()
table(nested_case1$Exclusivity)

View(nested_case1[nested_case1$Dissimilarity == 1, ])

View(nested_case1[nested_case1$Dissimilarity > 0.5, ])


#### b) Conservatism to fully weighted divergence ####
weighted_case1 <- 
  simdata2 %>%
  #mutate(mid_A = (a_A + b_A)/2, mid_B = (a_B + b_B)/2) %>%
  filter(a_A == a_B & b_A == b_B &
           alpha_A >= 1 & alpha_B >= 1 & gamma_A >= 1 & gamma_B >= 1)

100*nrow(weighted_case1)/nrow(simdata2) # Nested Case 1 represents 0.09% of all simulations

png("sim_plot_weighted1.png", width = 2000, height = 1000, res = 300)
weighted_case1 %>%
  ggplot(aes(x = Exclusivity, y = Dissimilarity), col = "red") +
  geom_point() +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw() +
  geom_vline(aes(xintercept = 0.5), color = "red") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_hline(aes(yintercept = 0.5), color = "red")
dev.off()

#### c) Conservatism to soft and full hard divergence ####

csh_case1 <- 
  simdata2 %>%
  #mutate(mid_A = (a_A + b_A)/2, mid_B = (a_B + b_B)/2) %>%
  filter((b_A - a_A) == (b_B - a_B) & 
           alpha_A == gamma_A & alpha_B == gamma_B &
           alpha_A >= 1 & alpha_B >= 1 & gamma_A >= 1 & gamma_B >= 1)

100*nrow(csh_case1)/nrow(simdata2) # Nested Case 1 represents 0.0072% of all simulations

png("sim_plot_csh1.png", width = 2000, height = 1000, res = 300)
csh_case1 %>%
  ggplot(aes(x = Exclusivity, y = Dissimilarity), col = "red") +
  geom_point() +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw() +
  geom_vline(aes(xintercept = 0.5), color = "red") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_hline(aes(yintercept = 0.5), color = "red")
dev.off()
#### d) Conservatism and moving nested divergence to hard divergence ####

nested_case2 <- 
  simdata2 %>%
  #mutate(mid_A = (a_A + b_A)/2, mid_B = (a_B + b_B)/2) %>%
  filter(a_A < a_B & b_A > b_B & alpha_A == gamma_A & alpha_B == gamma_B &
           alpha_A >= 1 & alpha_B >= 1 & gamma_A >= 1 & gamma_B >= 1)

100*nrow(nested_case2)/nrow(simdata2) # Nested Case 1 represents 0.01% of all simulations

png("sim_plot_nested2.png", width = 2000, height = 1000, res = 300)
nested_case2 %>%
  ggplot(aes(x = Exclusivity, y = Dissimilarity), col = "red") +
  geom_point() +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw() +
  geom_vline(aes(xintercept = 0.5), color = "red") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_hline(aes(yintercept = 0.5), color = "red")
dev.off()

#### e) Conservatism and moving weighted divergence to hard divergence ####
weighted_case2 <- 
  simdata2 %>%
  #mutate(mid_A = (a_A + b_A)/2, mid_B = (a_B + b_B)/2) %>%
  filter((b_A - a_A) == (b_B - a_B) & 
           a_A <= a_B & b_A <= b_B & a_B < b_A &
           alpha_A >= 1 & alpha_B >= 1 & gamma_A >= 1 & gamma_B >= 1 &
           alpha_A < gamma_A & alpha_B > gamma_B)

100*nrow(weighted_case2)/nrow(simdata2) # Nested Case 1 represents 1.26% of all simulations

png("sim_plot_weighted2.png", width = 2000, height = 1000, res = 300)
weighted_case2 %>%
  ggplot(aes(x = Exclusivity, y = Dissimilarity), col = "red") +
  geom_point() +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw() +
  geom_vline(aes(xintercept = 0.5), color = "red") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_hline(aes(yintercept = 0.5), color = "red")
dev.off()





##### Combined plot ####
plotdata <- bind_rows(nested_case1 %>% mutate(source = "Nested"), 
                      nested_case2 %>% mutate(source = "Nested"), 
                      weighted_case1 %>% mutate(source = "Weighted"),
                      weighted_case2 %>% mutate(source = "Weighted"),
                      csh_case1 %>% mutate(source = "Soft"))



#png("Figures/Simulations_Nested_Soft_Weight.png", height = 5000, width = 6500, res = 600)
tiff("Figures/Simulations_Nested_Soft_Weight.tif", height = 5100, width = 6500, res = 600)
plotdata %>%
  #filter(source %in% c("Nested", "Soft", "Weighted")) %>%
  ggplot(aes(x = Exclusivity, y = Dissimilarity, fill = source)) +
  geom_point(size = 3, shape = 21, color = "black") +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw() +
  xlab("Niche Exclusivity") + ylab("Niche Dissimilarity") +
  scale_fill_manual("Niche Divergence Gradients (To Hard Divergence)", 
                    #values = c("#62A39F", "#000000", "#2F8745")) +
                    values = brewer.pal(8, "Set1")[c(1, 2, 7)]) +
  #values = c("#000000", "#62A39F", "#725745")) +
  theme(text = element_text(size = 20),
        legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(size=8))) +
  ggtitle("Figure 3")
dev.off()
