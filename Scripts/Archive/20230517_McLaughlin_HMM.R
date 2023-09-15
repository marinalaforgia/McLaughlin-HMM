#### HMMs McLaughlin Manuscript ####

rm(list = ls())

calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

# 12-16-2022 Notes: 
## Using quadrats, more data, sufficiently independent
## filtering out HMMs that did not converge in 150 iterations
## only using HMM estimates when we have data from 80 or more plots (20% of dataset), ideally we would have 100 plots, but 20% of the dataset seems reasonable (ref)
## Check out Ostats and Codyn
## look at how trait space varies between rare and common species
## why not look at SB data from mclaughlin from exotic forbs

# compare traits across species and traits and talk about how trait distributions are diffierent or not between rare and common species within the datasets across sites
# use traits for species we have less data for if they overlap in trait space for those we had to drop from the HMM
# murray & connor paper for hierarchical partitioning analysis 

#### Load Libraries ####
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(rdacca.hp)

# Read in Data ####

mcl.df.q <- readRDS("Scripts/HMMs/McLaughlin/20221212_mcl_HMM_quadrat.RDS") # quadrat level HMM
#mcl.df.q <- readRDS("Scripts/HMMs/McLaughlin/20230704_mcl_HMM_quadrat.RDS") # quadrat level HMM
seed.species <- read.csv("Data/20211001_Full-Species-List.csv")
seed.traits <- read.csv("Data/20230607_Seed-Traits_clean_site.csv")

rm <- c("Logfia spp./Micropus californicus", "Aegilops triuncialis", "Trifolium bifidum/gracilentum", "Torilis arvensis/nodosa", "Hesperolinon spp.", "Lactuca serriola", "Cryptantha hispidula", "Cuscuta californica", "Silene gallica", "Erodium botrys", "Erodium brachycarpum", "Trifolium dubium", "Trifolium microdon", "Trifolium microcephalum") # rm commonly mis-ID'd species, species that hadnt always been censused (cuscuta) and goatgrass which is treated at the reserve
# 
# mcl <- filter(mcl, Species_Name %in% adult.traits[adult.traits$Annual.Perennial == "Annual",]$Species_Name, !(Species_Name %in% rm)) 
# mcl <- merge(serp[,1:2], mcl, by = "Site")
# 
# mcl[mcl$Species_Name == "Geranium molle",]$Species_Name <- "Geranium dissectum" # Susan confirmed these should all be dissectum
# 
# mcl <- mcl %>%
#   dplyr::group_by(Year, Site, Quadrat, Serpentine, Species_Name) %>%
#   dplyr::summarize(Cover = sum(Cover)) # re-summarize by species after merging geraniums

mcl.df.q$Species_Name <- recode_factor(mcl.df.q$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis")

mcl.df.q[mcl.df.q$Species_Name == "Daucus pusillus",]$FunGroup <- "Native Forb"

mcl.df.q$FunGroup <- recode_factor(mcl.df.q$FunGroup, 'Exotic Forb' = "non-native forb", 'Exotic Grass' = "non-native grass", 'Native Forb' = "native forb")
mcl.df.q <- filter(mcl.df.q, !Species_Name %in% rm, iter < 150, n.plots >= 50, FunGroup != "Native Grass") 

#### Seed trait prep ####
seed.traits2 <- filter(seed.traits, grepl("MCL", seed.traits$site, fixed = TRUE))

seed.traits <- filter(seed.traits, Species %in% mcl.df.q$Species_Name, site != "CP")

seed.traits <- unique(rbind(seed.traits, seed.traits2))
rm(seed.traits2)
# AMSMEN has a suspiciously low prop.C compared to other Amsinckia species, replace with mean
seed.traits[seed.traits$Species == "Amsinckia menziesii",]$prop.C <- mean(seed.traits$prop.C, na.rm = T)

# still some missing 
seed.traits[is.na(seed.traits$prop.C),]$prop.C <- mean(seed.traits$prop.C, na.rm = T)
seed.traits[is.na(seed.traits$prop.N),]$prop.N <- mean(seed.traits$prop.N, na.rm = T)

seed.traits$coat.thick.per.width <- seed.traits$both.thick/seed.traits$width

seed.traits$nat.inv <- recode_factor(seed.traits$nat.inv, invasive = "non-native")

seed.traits$FunGroup <- paste(seed.traits$nat.inv, seed.traits$group, sep = " ")

#### Q1: C/S Tradeoff ####
##### Figure ####
fig1 <- ggplot(mcl.df.q, aes(x = s, y = c, col = FunGroup)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point(aes(shape = FunGroup), size = 3) +
  #geom_text(aes(label = Code)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.83, 0.87),
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.box.background = element_rect(colour = "black")
  ) +
  scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "probability of temporal dispersal",
       y = "probability of spatial dispersal")

#ggsave("McLaughlin/Figures/C-S.jpeg", fig1, height = 5, width = 6, units = "in", dpi = 600)

##### Models ####
par(mfrow = c(2,2))

cor.test(mcl.df.q$s, mcl.df.q$c)

cor.test(mcl.df.q[mcl.df.q$FunGroup == "non-native forb",]$s, 
         mcl.df.q[mcl.df.q$FunGroup == "non-native forb",]$c)

cor.test(mcl.df.q[mcl.df.q$FunGroup == "native forb",]$s, 
         mcl.df.q[mcl.df.q$FunGroup == "native forb",]$c)

Q1sc.m <- lm(s ~ c, data = mcl.df.q)
plot(Q1sc.m)
summary(Q1sc.m)

Q1s.m <- lm(s ~ FunGroup, data = mcl.df.q)
plot(Q1s.m)
pairs(emmeans(Q1s.m, ~ FunGroup), adjust = "BH")

Q1c.m <- lm(c ~ FunGroup, data = mcl.df.q)
plot(Q1c.m)
pairs(emmeans(Q1c.m, ~ FunGroup), adjust = "BH")

par(mfrow = c(1,1))

#### Q2: Forb PCA ####

##### Prep ####
# trait <- c("Species", "new.code", "family", "group", "nat.inv", "FunGroup", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "E.S", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.width", "wing.loading", "ldd.natural", "ldd.all", "disp.cat.all", "disp.cat.nat")
# 
# trait.pca.forb <- seed.traits[seed.traits$group == "forb", trait]
# 
# #hist(log(trait.pca.forb$wing.loading))
# trait.pca.forb$wing.loading <- log(trait.pca.forb$wing.loading)
# 
# #hist(log(trait.pca.forb$cn))
# #trait.pca.forb$cn <- log(trait.pca.forb$cn)
# 
# #hist(log(trait.pca.forb$prop.C))
# trait.pca.forb$prop.C <- log(trait.pca.forb$prop.C)
# 
# #hist(trait.pca.forb$prop.N)
# #hist(log(trait.pca.forb$prop.N))
# trait.pca.forb$prop.N <- log(trait.pca.forb$prop.N)
# 
# #hist(log(trait.pca.forb$coat.perm.perc))
# trait.pca.forb$coat.perm.perc <- log(trait.pca.forb$coat.perm.perc)
# 
# #hist(log(trait.pca.forb$coat.thick.per.width))
# trait.pca.forb$coat.thick.per.width <- log(trait.pca.forb$coat.thick.per.width)
# 
# #hist(log(trait.pca.forb$morph.mass.mg))
# trait.pca.forb$morph.mass.mg <- log(trait.pca.forb$morph.mass.mg)
# 
# #hist(log(trait.pca.forb$chem.mass.mg))
# trait.pca.forb$chem.mass.mg <- log(trait.pca.forb$chem.mass.mg)
# 
# #hist(trait.pca.forb$size.mm)
# trait.pca.forb$size.mm <- log(trait.pca.forb$size.mm)
# 
# #hist(trait.pca.forb$set.time.mpsec)
# 
# #hist(log(trait.pca.forb$E.S))
# trait.pca.forb$E.S <- log(trait.pca.forb$E.S)
# 
# ##### PCA ####
# # PROP N AND E:S removed because they just muddied the figure and explain no variation
# trait.pca <- c("chem.mass.mg", "shape", "prop.C", "set.time.mpsec", "coat.perm.perc", "coat.thick.per.width", "size.mm", "ldd.natural", "prop.N", "E.S")
# 
# pca.forb <- prcomp(trait.pca.forb[, trait.pca], scale = T, center = T)
# summary(pca.forb)
# screeplot(pca.forb)
# (ev <- pca.forb$sdev^2)
# trait.pca.forb <- cbind(trait.pca.forb, pca.forb$x[,1:4])
# 
# autoplot(pca.forb, x = 1, y = 2, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", shape = "FunGroup", size = 2.5) +
#   theme_bw() +
#   theme(
#     #panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#     legend.title = element_blank(),
#     #legend.position = c(.2,.15),
#     legend.position = "bottom",
#     legend.text = element_text(size = 14),
#     plot.margin = unit(c(1,0.5,1,0.5), "cm"),
#     #legend.box.background = element_rect(colour = "black", linewidth = 0.5),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14)
#   )  +
#   scale_color_manual(values = c("#7570b3", "#d95f02"))

#ggsave("McLaughlin/Figures/PCA.jpeg", fig.pca, height = 5, width = 6, units = "in", dpi = 600)

# autoplot(pca.forb, x = 3, y = 4, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "nat.inv", shape = "nat.inv", size = 2.5) +
# theme_classic() +
# theme(
#   panel.border = element_rect(colour = "black", fill = NA, size = 1),
#   legend.position = "none"
#  )  +
# scale_color_manual(values = c("#d95f02", "#7570b3"))
# 
# autoplot(pca.forb, x = 5, y = 6, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "nat.inv", shape = "nat.inv", size = 2.5) +
# theme_classic() +
# theme(
#   panel.border = element_rect(colour = "black", fill = NA, size = 1),
#   legend.position = "none"
#  )  +
# scale_color_manual(values = c("#d95f02", "#7570b3"))

##### Nicer PCA fig ####
#' pca_load <- 
#'   as_tibble(pca.forb$rotation, rownames = 'variable') %>% 
#'   #we can rename the variables so they look nicer on the figure
#'   mutate(variable = dplyr::recode(variable,
#'                   'set.time.mpsec' = 'settling time',
#'                   'chem.mass.mg' = 'mass',
#'                   'coat.thick.per.width' = 'coat thickness',
#'                   #'prop.N' = 'N',
#'                   'prop.C' = 'C',
#'                   'ldd.natural' = 'dispersal',
#'                   'size.mm' = 'size',
#'                   #'E.S' = 'E:S',
#'                   'coat.perm.perc' = 'coat permeability',
#'   
#' ))
#' 
#' fig.pca <- ggplot(trait.pca.forb, aes(x = PC1, y = PC2)) +
#'   geom_point(aes(colour = FunGroup, shape = FunGroup), size = 2.5) +
#'   theme_light() + 
#'   geom_segment(data = pca_load, 
#'                aes(x = 0, y = 0, 
#'                    xend = PC1*5,
#'                    yend = PC2*5),
#'                arrow = arrow(length = unit(1/2, 'picas'))) +
#'   annotate('text', 
#'            x = (pca_load$PC1*5.8), 
#'            y = ifelse(pca_load$PC2 > 0, pca_load$PC2*5.4, pca_load$PC2*7),
#'            label = pca_load$variable,
#'            size = 4)  +
#'   theme(
#'     #panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#'     legend.title = element_blank(),
#'     #legend.position = c(.2,.15),
#'     legend.position = "bottom",
#'     legend.text = element_text(size = 14),
#'     plot.margin = unit(c(1,0.5,1,0.5), "cm"), 
#'     #legend.box.background = element_rect(colour = "black", linewidth = 0.5),
#'     axis.text = element_text(size = 12),
#'     axis.title = element_text(size = 14)  
#'   )  +
#'   scale_color_manual(values = c("#7570b3", "#d95f02")) +
#'   labs(x = "PC1 (26.38%)", y = "PC2 (19.44)")

#ggsave("Manuscript/McLaughlin/Figures/PCA.jpeg", fig.pca, height = 5, width = 6, units = "in", dpi = 600)

# #### PCA: nonnatives ###
# trait.pca <- c("chem.mass.mg", "shape", "prop.C", "set.time.mpsec", "coat.perm.perc", "coat.thick.per.width", "size.mm", "ldd.natural", "prop.N", "E.S")
# 
# pca.i <- prcomp(trait.pca.forb[trait.pca.forb$nat.inv == "non-native", trait.pca], scale = T, center = T)
# trait.pca.forb.i <- cbind(trait.pca.forb[trait.pca.forb$nat.inv == "non-native",], pca.i$x[,1:4])
# 
# tmp1 <- autoplot(pca.i, x = 1, y = 2, data = trait.pca.forb.i, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", shape = "FunGroup", size = 2.5) +
#   theme_bw() +
#   theme(
#     #panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#     legend.title = element_blank(),
#     #legend.position = c(.2,.15),
#     legend.position = "bottom",
#     legend.text = element_text(size = 14),
#     plot.margin = unit(c(1,0.5,1,0.5), "cm"),
#     #legend.box.background = element_rect(colour = "black", linewidth = 0.5),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14)
#   )  +
#   scale_color_manual(values = c("#7570b3", "#d95f02"))
# 
# pca.n <- prcomp(trait.pca.forb[trait.pca.forb$nat.inv == "native", trait.pca], scale = T, center = T)
# trait.pca.forb.n <- cbind(trait.pca.forb[trait.pca.forb$nat.inv == "native",], pca.n$x[,1:4])
# 
# tmp2 <- autoplot(pca.n, x = 1, y = 2, data = trait.pca.forb.n, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", shape = "FunGroup", size = 2.5) +
#   theme_bw() +
#   theme(
#     #panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
#     legend.title = element_blank(),
#     #legend.position = c(.2,.15),
#     legend.position = "bottom",
#     legend.text = element_text(size = 14),
#     plot.margin = unit(c(1,0.5,1,0.5), "cm"),
#     #legend.box.background = element_rect(colour = "black", linewidth = 0.5),
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14)
#   )  +
#   scale_color_manual(values = c("#7570b3", "#d95f02"))
# 
# ggarrange(tmp1, tmp2)
# 
# mcl.df.q.forb.PCi <- merge(mcl.df.q.forb[,-13], trait.pca.forb.i, by.x = "Species_Name", by.y = "Species", all = F)
# 
# ggplot(mcl.df.q.forb.PCi, aes(x = PC1, y = s)) +
#   geom_smooth(method = "lm") +
#   geom_point()
# 
# ggplot(mcl.df.q.forb.PCi, aes(x = PC1, y = c)) +
#   geom_smooth(method = "lm") +
#   geom_point()
# 
# ggplot(mcl.df.q.forb.PCi, aes(x = PC2, y = s)) +
#   geom_smooth(method = "lm") +
#   geom_point()
# 
# ggplot(mcl.df.q.forb.PCi, aes(x = PC2, y = c)) +
#   geom_smooth(method = "lm") +
#   geom_point()

#### Q3: HMMs & PCs ####
mcl.df.q.forb <- filter(mcl.df.q, FunGroup == "native forb" | FunGroup == "non-native forb")

mcl.df.q.forb.PC <- merge(mcl.df.q.forb[,-13], trait.pca.forb, by.x = "Species_Name", by.y = "Species", all = F)


###### Figure ####

## PC1 ##
a <- ggplot(mcl.df.q.forb.PC, aes(x = PC1, y = c, col = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2, linetype = 2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 14)  
  ) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(-.005,0.26)) +
  labs(y = "spatial dispersal")

b <- ggplot(mcl.df.q.forb.PC, aes(x = PC1, y = s, col = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2, linetype = 2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)    
  ) +
  scale_color_manual(values = c("#7570b3", "#d95f02"))  +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_y_continuous(breaks = c(0.4,0.6,0.8), limits = c(0.2,0.85)) +
  labs(y = "temporal dispersal")

## PC 2 ##
c <- ggplot(mcl.df.q.forb.PC, aes(x = PC2, y = c, col = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$FunGroup == "non-native forb",],method = "lm", aes(fill = FunGroup), linetype = 2, alpha = 0.2, fill = "#d95f02") + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank()   
  ) +
  scale_color_manual(values = c("#7570b3", "#d95f02"))  +
  #scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_y_continuous(breaks = c(0.0, 0.1, 0.2), limits = c(-.005,0.26)) +
  scale_linetype_manual(values = c(1,2)) +
  labs(y = "spatial dispersal") 

d <- ggplot(mcl.df.q.forb.PC, aes(x = PC2, y = s, col = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$FunGroup == "native forb",],method = "lm", aes(fill = FunGroup), linetype = 1, alpha = 0.2, fill = "#7570b3") + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank()
  ) +
  scale_color_manual(values = c("#7570b3", "#d95f02"))  +
  #scale_fill_manual(values = c("#7570b3", "#d95f02")) +  
  scale_y_continuous(breaks = c(0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_linetype_manual(values = c(1,2)) +
  labs(y = "temporal dispersal") 

PC.fig <- ggarrange(a,c,b,d, ncol = 2, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)"))

#annotate_figure(figure, left = textGrob("Common y-axis", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),bottom = textGrob("Common x-axis", gp = gpar(cex = 1.3)))

#ggsave("McLaughlin/Figures/HMM-PC.jpeg", PC.fig, height = 10, width = 11.75, units = "in", dpi = 600)


PCA.panel <- ggarrange(fig.pca, 
                       ggarrange(a,c,b,d, ncol = 2, nrow = 2, 
                                 labels = c("(b)","(c)","(d)","(e)"),
                                 widths = c(1, 0.85), 
                                 heights = c(0.9,1), 
                                 label.x = c(0.16,0.04,0.16,0.04), 
                                 label.y = 0.96), 
                       widths = c(0.8, 1), 
                       ncol = 2, 
                       labels = "(a)", 
                       label.x = 0.12, 
                       label.y = 0.92)

ggsave("Manuscript/McLaughlin/Figures/PCA-panel.jpeg", PCA.panel, height = 6, width = 12, units = "in", dpi = 600)

### PC3 ##
# ggplot(mcl.df.q.forb.PC, aes(x = PC3, y = c, col = nat.inv, shape = nat.inv)) +
#   geom_smooth(method = "lm") +
#   #geom_text(aes(label = Code)) +
#   geom_point(size = 2) +
#   theme_bw() +
#   theme(
#     legend.position = "none"
#   ) +
#   scale_color_manual(values = c("#d95f02", "#7570b3"))  +
#   labs(y = "probability of spatial dispersal")
# 
# ggplot(mcl.df.q.forb.PC, aes(x = PC3, y = s, col = nat.inv, shape = nat.inv)) +
#   geom_smooth(method = "lm") +
#   #geom_text(aes(label = Code)) +
#   geom_point(size = 2) +
#   theme_bw() +
#   theme(
#     legend.position = "none"
#   ) +
#   scale_color_manual(values = c("#d95f02", "#7570b3")) +
#   labs(y = "probability of temporal dispersal")
# 
# ## PC4 ##
# ggplot(mcl.df.q.forb.PC, aes(x = PC4, y = c, col = nat.inv, shape = nat.inv)) +
#   geom_smooth(method = "lm") +
#   #geom_text(aes(label = Code)) +
#   geom_point(size = 2) +
#   theme_bw() +
#   theme(
#     legend.position = "none"
#   ) +
#   scale_color_manual(values = c("#d95f02", "#7570b3"))  +
#   labs(y = "probability of spatial dispersal")
# 
# ggplot(mcl.df.q.forb.PC, aes(x = PC4, y = s, col = nat.inv, shape = nat.inv)) +
#   geom_smooth(method = "lm") +
#   #geom_text(aes(label = Code)) +
#   geom_point(size = 2) +
#   theme_bw() +
#   theme(
#     legend.position = "none"
#   ) +
#   scale_color_manual(values = c("#d95f02", "#7570b3")) +
#   labs(y = "probability of temporal dispersal")

##### Models ####
## by PC ##
# Note 5/23: some outliers and violation of model assumptions, explore more

# cortests PC1
cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$s,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$PC1)

cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$c,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$PC1)

cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$s,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$PC1)

cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$c,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$PC1)

# COR TEST PC2
cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$s,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$PC2)

cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$c,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$PC2)

cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$s,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$PC2)

cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$c,
         mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",]$PC2)
 
##
par(mfrow = c(2,2))
m.s.n.PC1 <- lm(s ~ PC1, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",])
plot(m.s.n.PC1)
summary(m.s.n.PC1)

m.s.i.PC1 <- lm(s ~ PC1, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",])
plot(m.s.i.PC1)
summary(m.s.i.PC1)

m.c.n.PC1 <- lm(c ~ PC1, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",])
plot(m.c.n.PC1)
summary(m.c.n.PC1)

m.c.i.PC1 <- lm(c ~ PC1, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",])
plot(m.c.i.PC1)
summary(m.c.i.PC1)

m.s.n.PC2 <- lm(s ~ PC2, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",])
plot(m.s.n.PC2)
summary(m.s.n.PC2)

m.s.i.PC2 <- lm(s ~ PC2, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",])
plot(m.s.i.PC2)
summary(m.s.i.PC2)

m.c.n.PC2 <- lm(c ~ PC2, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",])
plot(m.c.n.PC2)
summary(m.c.n.PC2)

m.c.i.PC2 <- lm(c ~ PC2, data = mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "non-native",])
plot(m.c.i.PC2)
summary(m.c.i.PC2)

par(mfrow = c(1,1))

# Hierarchical Partitioning ####
mcl.df.q.trait <- merge(mcl.df.q[,-13], seed.traits, by.x = "Species_Name", by.y = "Species", all = F)

## Output
#the relative importance of any individual predictor can be simply estimated as its unique contribution to the total model R2 plus its average shared contributions with the other predictors.

# Unique: part of variation ONLY explained by that variable 

# Average.share: total shared/number of predictors sharing; Shared fractions represent the variation in the response variable explained by the correlation of the predictors involved (i.e. correlation in the predictive space). Larger shared fractions imply that more multicollinearity is present in the model. #Negative common (shared) variation is possible when predictors act as suppressors of other predictors

# Individual = Relative importance of individual predictor = unique + shared portion where shared equals = total shared/number of predictors sharing; sum of this adds up to R2

# I.perc = adds up up to 100%

# I *THINK* RDA works best with normalized data
#https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html

##### Prep ####

# should only the traits that can be directly related to survival go in? like height and settling time can only really be indirectly related

survival <- c("chem.mass.mg",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.width")

colonization <- c("shape", "chem.mass.mg", "set.time.mpsec","size.mm", "ldd.natural","prop.C", "prop.N")

# pick traits without correlations but RDA might be able to handle them together

# how is this different from a hierarchical model? order matters
# conceptually it's fairly similar 

##### All plants ####

###### S: All plants ####

HPA.s <- rdacca.hp(
  dv = mcl.df.q.trait$s,
  iv = mcl.df.q.trait[,colnames(mcl.df.q.trait) %in% survival],
  method = "RDA",
  type = "adjR2"
)

HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation # 18% variation

plot(HPA.s) # size and shape

###### C: All plants ####

HPA.c <- rdacca.hp(
  dv = mcl.df.q.trait$c,
  iv = mcl.df.q.trait[ ,colnames(mcl.df.q.trait) %in% colonization],
  method = "RDA",
  type = "adjR2"
)

HPA.c$Hier.part

HPA.c$Total_explained_variation #23%

plot.rdaccahp(HPA.c)

##### All Forbs ####

###### S: All Forbs ####

# HPA.s <- rdacca.hp(
#   dv = mcl.df.q.forb.PC$s,
#   iv = mcl.df.q.forb.PC[,colnames(mcl.df.q.forb.PC) %in% survival],
#   method = "RDA",
#   type = "adjR2",
#   var.part = T
# )
# 
# HPA.s$Hier.part
# #HPA.s$Var.part
# HPA.s$Total_explained_variation
# 
# plot.rdaccahp(HPA.s) # size and shape

###### C: All Forbs ####

# HPA.c <- rdacca.hp(
#   dv = mcl.df.q.forb.PC$c,
#   iv = mcl.df.q.forb.PC[ ,colnames(mcl.df.q.forb.PC) %in% colonization],
#   method = "RDA",
#   type = "adjR2"
# )
# 
# HPA.c$Hier.part
# 
# HPA.c$Total_explained_variation
# 
# plot.rdaccahp(HPA.c)

##### Natives #####

###### S: Natives #####

test <- rdacca.hp(
  dv = mcl.df.q.trait[mcl.df.q.trait$nat.inv == "native",]$s,
  iv = as.data.frame(scale(mcl.df.q.trait[mcl.df.q.trait$nat.inv == "native", colnames(mcl.df.q.trait) %in% survival])),
  
)

test$Total_explained_variation

test$Hier.part #20%
#test$Var.part

plot.rdaccahp(test) # shape, size

###### C: Natives #####
test <- rdacca.hp(
  dv = mcl.df.q.trait[mcl.df.q.trait$nat.inv == "native",]$c,
  iv = as.data.frame(scale(mcl.df.q.trait[mcl.df.q.trait$nat.inv == "native",colnames(mcl.df.q.trait) %in% colonization])),
)

test$Total_explained_variation #14%
test$Hier.part
plot.rdaccahp(test) # ldd.natural, size


##### Invasives #####

###### S: Invasives #####
test <- rdacca.hp(
  dv = mcl.df.q.trait[mcl.df.q.trait$nat.inv == "non-native",]$s,
  iv = as.data.frame(scale(mcl.df.q.trait[mcl.df.q.trait$nat.inv == "non-native",colnames(mcl.df.q.trait) %in% survival])),
)

test$Total_explained_variation
test$Hier.part
plot.rdaccahp(test) # coat thick per width, prop.N

###### C: Invasives #####

test <- rdacca.hp(
  dv = mcl.df.q.trait[mcl.df.q.trait$nat.inv == "non-native",]$c,
  iv = as.data.frame(scale(mcl.df.q.trait[mcl.df.q.trait$nat.inv == "non-native", colnames(mcl.df.q.trait) %in% colonization])),
  method = "RDA",
  type = "adjR2"
)

test$Total_explained_variation
test$Hier.part
plot.rdaccahp(test) # size, set time

##### Grasses #####

###### S: Invasives #####
test <- rdacca.hp(
  dv = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s,
  iv = as.data.frame(scale(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",colnames(mcl.df.q.trait) %in% survival])),
)

test$Total_explained_variation
test$Hier.part
plot.rdaccahp(test) # coat thick per width, prop.N

###### C: Invasives #####

test <- rdacca.hp(
  dv = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$c,
  iv = as.data.frame(scale(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",colnames(mcl.df.q.trait) %in% colonization])),
)

test$Total_explained_variation
test$Hier.part
plot.rdaccahp(test) # size, set time

#### Q4: HMMs & Traits ####

ggplot(mcl.df.q.forb.PC, aes(x = s, y = c, col = appendage)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point(aes(shape = FunGroup), size = 3) +
  #geom_text(aes(label = Code)) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.box.background = element_rect(colour = "black")
  ) +
  #scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "probability of temporal dispersal",
       y = "probability of spatial dispersal")

# Traits identified in Hierarchical partitioning 
## Natives: s (shape, size), c (disp potential, size)
## Non-natives: s (coat thick, prop N), c (shape, size)

trait.pca.forb.ns <- seed.traits[seed.traits$group == "forb", trait]
mcl.df.q.forb.PC.ns <- merge(mcl.df.q.forb[,-13], trait.pca.forb.ns, by.x = "Species_Name", by.y = "Species", all = F)

###### Figures ####
##### .__MASS ####
h <- ggplot(mcl.df.q.forb.PC.ns, aes(x = log(chem.mass.mg), y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(y = "spatial dispersal", x = "log mass (mg)")

 #with grasses
 # ggplot(mcl.df.q.trait, aes(x = log(chem.mass.mg), y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
 #  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
 #  #geom_text(aes(label = Code)) +
 #  geom_point(size = 2.5) +
 #  theme_bw() +
 #  theme(
 #    axis.text = element_text(size = 12),
 #    axis.title = element_text(size = 14),
 #    legend.position = "none",
 #    legend.title = element_blank(),
 #    legend.text = element_text(size = 12)
 #  ) +
 #  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
 #  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
 #  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
 #  scale_shape_manual(values = c( 19, 17, 15)) +
 #  labs(y = "spatial dispersal", x = "log mass (mg)")
 
cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$chem.mass.mg) #ns

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$chem.mass.mg)

i <- ggplot(mcl.df.q.forb.PC.ns, aes(x = log(chem.mass.mg), y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.4,0.6,0.8), limits = c(0.2,0.95)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(y = "temporal dispersal", x = "log mass (mg)")

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$chem.mass.mg)

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$chem.mass.mg) #ns

# ggplot(mcl.df.q.forb.PC.ns, aes(x = log(chem.mass.mg), y = g, col = FunGroup, group = FunGroup, shape = FunGroup)) +
#   geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
#   #geom_text(aes(label = Code)) +
#   geom_point(size = 2.5) +
#   theme_bw() +
#   theme(
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14),
#     legend.position = "none",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12)
#   ) +
#   #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
#   scale_fill_manual(values = c("#7570b3", "#d95f02")) +
#   scale_color_manual(values = c("#7570b3", "#d95f02")) +
#   labs(y = "spatial dispersal")

##### .__SHAPE ####
# MICDOU - shape is tricky... does or does it not include the pappus? 
a <- ggplot(mcl.df.q.forb.PC.ns, aes(x = shape, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$nat.inv == "native",], method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(y = "spatial dispersal")

 # ggplot(mcl.df.q.trait, aes(x = shape, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
 #  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
 #  #geom_text(aes(label = Code)) +
 #  geom_point(size = 2.5) +
 #  theme_bw() +
 #  theme(
 #    axis.text = element_text(size = 12),
 #    axis.title = element_text(size = 14),
 #    legend.position = "none",
 #    legend.title = element_blank(),
 #    legend.text = element_text(size = 12)
 #  ) +
 #  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
 #  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
 #  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
 #  scale_shape_manual(values = c( 19, 17, 15)) +
 #  labs(y = "spatial dispersal")
 
cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$shape) 

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb" & mcl.df.q.forb.PC.ns$Code != "MICDOU",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb" & mcl.df.q.forb.PC.ns$Code != "MICDOU",]$shape) #ns

f <- ggplot(mcl.df.q.forb.PC.ns, aes(x = shape, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$nat.inv == "native",], method = "lm", fill = "#7570b3", alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.4,0.6,0.8), limits = c(0.2,0.95)) +
  #scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(y = "temporal dispersal")

 ggplot(mcl.df.q.trait, aes(x = shape, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "temporal dispersal")
 
cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$shape) #ns

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$shape) #sig

##### .__SIZE ####
b <- ggplot(mcl.df.q.forb.PC.ns, aes(x = log(size.mm), y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup, linetype = FunGroup), alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none", 
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +  
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  scale_linetype_manual(values = c(1,2)) +
  labs(y = "spatial dispersal", x = "log size (mm)")

 ggplot(mcl.df.q.trait, aes(x = log(size.mm), y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "temporal dispersal")
 
  ggplot(mcl.df.q.trait, aes(x = log(size.mm), y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")
  
cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$size.mm) #marg

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$size.mm) #sig

g <- ggplot(mcl.df.q.forb.PC.ns, aes(x = log(size.mm), y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",], method = "lm", aes(fill = FunGroup), linetype = 2,alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.4,0.6,0.8), limits = c(0.2,0.95)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_linetype_manual(values = c(1,2)) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(y = "temporal dispersal", x = "log size (mm)")

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$size.mm) #marg

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$size.mm) #ns

##### .__COAT ####

c <- ggplot(mcl.df.q.forb.PC.ns, aes(x = coat.thick.per.width, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$nat.inv == "non-native",], method = "lm", fill = "#d95f02", alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)  
    ) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_y_continuous(breaks = c(0.4,0.6,0.8), limits = c(0.2,0.95)) +
  labs(x = "seed coat to width ratio")

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$coat.thick.per.width) #sig

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$coat.thick.per.width) #ns

ggplot(mcl.df.q.forb.PC, aes(x = both.thick, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm") +
  geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(y = "temporal dispersal")

  ggplot(mcl.df.q.trait, aes(x = both.thick, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")
  
ggplot(mcl.df.q.trait, aes(x = coat.thick.per.width, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")
  
##### .__N ####
d <- ggplot(mcl.df.q.forb.PC.ns, aes(x = prop.N, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$nat.inv == "non-native",], method = "lm", fill = "#d95f02", alpha = 0.2, linetype = 1) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.4,0.6,0.8), limits = c(0.2,0.95), ) +
  scale_x_continuous(breaks = c(0.02,0.04,0.06), labels = scales::percent) +
  #scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(y = "temporal dispersal", x = "percent N")

ggplot(mcl.df.q.trait, aes(x = prop.N, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$prop.N) #marg

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$s, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$prop.N) #ns

###### .__LDD ####
e <- ggplot(mcl.df.q.forb.PC.ns, aes(x = ldd.natural, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$nat.inv == "native",], method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(x = "dispersal potential")

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$ldd.natural) #ns


cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$ldd.natural) #sig

ggplot(mcl.df.q.trait, aes(x = ldd.natural, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")

##### .__set time ####
ggplot(mcl.df.q.forb.PC.ns, aes(x = set.time.mpsec, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(x = "settling speed (m/s)")

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$set.time.mpsec) #ns


cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$set.time.mpsec) #sig

ggplot(mcl.df.q.trait, aes(x = set.time.mpsec, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")

##### .__C ####
ggplot(mcl.df.q.forb.PC.ns, aes(x = prop.C, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(x = "proportion C")


cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$prop.C) #ns


cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$c, 
         mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$prop.C) #sig

ggplot(mcl.df.q.trait, aes(x = prop.C, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")

##### .__ coat perm ####
ggplot(mcl.df.q.forb.PC.ns, aes(x = log(coat.perm.perc), y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3", "#d95f02")) +
  scale_color_manual(values = c("#7570b3", "#d95f02")) +
  labs(x = "seed coat permeability")

cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$s, 
         log(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "non-native forb",]$coat.perm.perc)) #ns


cor.test(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$s, 
         log(mcl.df.q.forb.PC.ns[mcl.df.q.forb.PC.ns$FunGroup == "native forb",]$coat.perm.perc)) #sig

ggplot(mcl.df.q.trait, aes(x = log(coat.perm.perc), y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")

ggplot(mcl.df.q.trait, aes(x = height.cm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")

leg <- get_legend(
  ggplot(mcl.df.q.forb.PC.ns, aes(x = ldd.natural, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
    geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 14)
    ) +
    geom_point(size = 2.5) + 
    scale_fill_manual(values = c("#7570b3", "#d95f02")) +
    scale_color_manual(values = c("#7570b3", "#d95f02")))

#### .__Panel Plot ####
# with mass
#trait.panel <- ggarrange(h,a,b,e,leg,i,f,g,c,d, ncol = 5, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)","" , "(e)", "(f)", "(g)", "(h)", "(i)"), widths = c(1, 0.85, 0.85, 0.85, 0.85), heights = c(1,1), label.x = c(0.2,0.04,0.04,0.04,0.04,0.2,0.04,0.04,0.04,0.04), label.y = 0.96)

trait.panel <- ggarrange(a,b,e,leg,f,g,c,d, ncol = 4, nrow = 2, labels = c("(a)", "(b)", "(c)","" , "(e)", "(f)", "(g)", "(h)"), widths = c(1, 0.85, 0.85, 0.85), heights = c(1,1), label.x = c(0.2,0.04,0.04,0.04,0.2,0.04,0.04,0.04), label.y = 0.96)

ggsave("Manuscript/McLaughlin/Figures/trait-panel.jpeg", trait.panel, height = 6, width = 10, units = "in", dpi = 600)


### RA Plot  ####
mcl.mean <- mcl %>%
  dplyr::group_by(Species_Name) %>%
  dplyr::summarize(cover.mean = mean(Cover, na.rm = T),
            cover.se = calcSE(Cover),
            cover.var = var(Cover, na.rm = T))

mcl2 <- merge(mcl.df.q, mcl.mean, by = "Species_Name")

ggplot(mcl2[mcl2$iter < 150 & mcl2$n.plots >= 80,], aes(x = s, y = c, size = cover.mean)) +
  geom_point() # variation here, but looks like more prevalent things typically have higher dispersal; more rare species hanging out at seed bank, need bet hedging or strong intra specific density dependent
# there might be colonizaiton attempts but harsh serp and harsh competition
# get back to spatial dispersal literature, venable and brown paper 

ggplot(mcl2[mcl2$iter < 150 & mcl2$n.plots >= 80,], aes(x = s, y = c, size = cover.var)) +
  geom_point() 

#### Spider Plots ####
# Library
library(fmsb)

## Colonization 
colonization <- c("c", "chem.mass.mg",  "shape", "set.time.mpsec", "height.cm", "size.mm")

spider.c.n <- cor(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native", colnames(mcl.df.q.forb.PC) %in% colonization])

trait.p <- data.frame()

for(i in colonization[-1]){
  tmp <- cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$c, mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",i])
  trait.p <- rbind(trait.p, c(trait = i, p.value = tmp$p.value))
}

colnames(trait.p) <- c("trait", "p.value")

spider.c.i <- cor(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive", colnames(mcl.df.q.forb.PC) %in% colonization])

trait.p <- data.frame()

for(i in colonization[-1]){
  tmp <- cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive",]$c, mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive",i])
  trait.p <- rbind(trait.p, c(trait = i, p.value = tmp$p.value))
}

colnames(trait.p) <- c("trait", "p.value")

## Survival
survival <- c("s", "chem.mass.mg",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.size")

spider.s.n <- cor(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native", colnames(mcl.df.q.forb.PC) %in% survival])

spider.s.i <- cor(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive", colnames(mcl.df.q.forb) %in% survival])

trait.p <- data.frame()

for(i in survival[-1]){
  tmp <- cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",]$s, mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native",i])
  trait.p <- rbind(trait.p, c(trait = i, p.value = tmp$p.value))
}

colnames(trait.p) <- c("trait", "p.value")

trait.p <- data.frame()

for(i in survival[-1]){
  tmp <- cor.test(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive",]$s, mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive",i])
  trait.p <- rbind(trait.p, c(trait = i, p.value = tmp$p.value))
}

colnames(trait.p) <- c("trait", "p.value")

# t <- cor.test(mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native",]$s, mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native",]$size.mm)$p.value
# 
# 
# 
spider.s.i <- as.data.frame(spider.s.i)
spider.s.i$nat.inv <- "invasive"

spider.s.n <- as.data.frame(spider.s.n)
spider.s.n$nat.inv <- "native"

spider.s.i <- as.data.frame(spider.s.i)
spider.s.i$nat.inv <- "invasive"
spider.s.i$trait <- row.names(spider.s.i)

spider.s.n <- as.data.frame(spider.s.n)
spider.s.n$nat.inv <- "native"
spider.s.n$trait <- row.names(spider.s.n)

cor.s <- as.data.frame(rbind(spider.s.i[-1,c(1,9,10)], spider.s.n[-1,c(1,9,10)]))
colnames(cor.s)[1] <- "s.cor"

### Correlations ####
ggplot(seed.traits[seed.traits$family != "Poaceae",], aes(x = size.mm, y = shape)) +
geom_point() +
geom_smooth(method = "lm") # not actually correlated 

ggplot(seed.traits[seed.traits$family != "Poaceae",], aes(x = size.mm, y = width)) +
geom_point() +
geom_smooth(method = "lm") # not actually correlated 

ggplot(seed.traits[seed.traits$family != "Poaceae",], aes(x = size.mm, y = shape)) +
geom_point() +
geom_smooth(method = "lm") # not actually correlated 

# 
# ggplot(seed.traits, aes(x = both.thick, y = chem.mass.mg)) +
# geom_point() +
# geom_smooth(method = "lm")
# 
# ggplot(seed.traits, aes(x = both.thick, y = size.mm)) +
# geom_point() +
# geom_smooth(method = "lm")
# 
# ggplot(seed.traits, aes(x = coat.thick.per.mass, y = chem.mass.mg)) +
# geom_point() +
# geom_smooth(method = "lm")



