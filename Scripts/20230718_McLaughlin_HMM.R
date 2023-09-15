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
## look at how trait space varies between rare and common species
## why not look at SB data from mclaughlin from exotic forbs

# compare traits across species and traits and talk about how trait distributions are different or not between rare and common species within the datasets across sites


#### Load Libraries ####
library(ggplot2)
library(ggfortify)
library(tidyverse)
library(ggpubr)
library(emmeans)
library(rdacca.hp)

# Read in Data ####

mcl.df.q <- readRDS("Data/20221212_mcl_HMM_quadrat.RDS") # quadrat level HMM
seed.species <- read.csv("Data/20211001_Full-Species-List.csv")
seed.traits <- read.csv("Data/20230801_Seed-Traits_clean_site.csv")

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

mcl.df.q$Species_Name <- recode_factor(mcl.df.q$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis", 'Gastridium phleoides' = "Gastridium ventricosum")

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

seed.traits$nat.inv <- recode_factor(seed.traits$nat.inv, invasive = "non-native")

seed.traits$FunGroup <- paste(seed.traits$nat.inv, seed.traits$group, sep = " ")

seed.traits$new.mass <- ifelse(seed.traits$group == "grass", seed.traits$morph.mass.mg, seed.traits$chem.mass.mg)

# Normalize seed.traits
#hist(log(seed.traits$wing.loading))
seed.traits$wing.loading <- log(seed.traits$wing.loading)

#hist(log(seed.traits$coat.perm.perc))
seed.traits$coat.perm.perc <- log(seed.traits$coat.perm.perc)

#hist(log(seed.traits$morph.mass.mg))
seed.traits$morph.mass.mg <- log(seed.traits$morph.mass.mg)

hist(exp(seed.traits$chem.mass.mg))
hist(seed.traits$chem.mass.mg)
seed.traits$chem.mass.mg <- log(seed.traits$chem.mass.mg)

#hist(log(seed.traits$new.mass))
seed.traits$new.mass <- log(seed.traits$new.mass)

#hist(log(seed.traits$size.mm))
seed.traits$size.mm <- log(seed.traits$size.mm)

#### Figure 1: Trade-off ####
##### Figure ####
fig1.df <- mcl.df.q
fig1.df$s1 <- NA
fig1.df$c1 <- NA
fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$s1 <- fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$s
fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$c1 <- fig1.df[fig1.df$Species_Name == "Centaurea solstitialis",]$c

fig1 <- ggplot(fig1.df, aes(x = s, y = c, col = FunGroup, shape = FunGroup)) +
  geom_smooth(inherit.aes = F, aes(x = s, y = c), method = "lm", col = "black", se = F) + 
  geom_point(size = 3) +
  geom_point(aes(x = s1, y = c1, col = FunGroup, shape = FunGroup), size = 3) +
  #geom_jitter(aes(shape = FunGroup), size = 3, width = 0.01) +
  #geom_text(aes(label = Code)) +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linewidth = 1, fill = NA),
    legend.position = "right",
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    #legend.box.background = element_rect(colour = "black")
  ) +
  scale_color_manual(values = c("#d95f02", "#1b9e77", "#7570b3")) +
  scale_shape_manual(values = c(17, 15, 19)) +
  labs(x = "probability of temporal dispersal",
       y = "probability of spatial dispersal")

ggsave("Manuscript/McLaughlin/Figures/C-S.tiff", fig1, height = 5, width = 7, units = "in", dpi = 600)

##### Models ####
par(mfrow = c(2,2))

cor.test(mcl.df.q$s, mcl.df.q$c)

cor.test(mcl.df.q[mcl.df.q$FunGroup == "non-native forb",]$s, 
         mcl.df.q[mcl.df.q$FunGroup == "non-native forb",]$c)

cor.test(mcl.df.q[mcl.df.q$FunGroup == "native forb",]$s, 
         mcl.df.q[mcl.df.q$FunGroup == "native forb",]$c)

Q1s.m <- lm(s ~ FunGroup, data = mcl.df.q)
plot(Q1s.m)
pairs(emmeans(Q1s.m, ~ FunGroup), adjust = "BH")

Q1c.m <- lm(c ~ FunGroup, data = mcl.df.q)
plot(Q1c.m)
pairs(emmeans(Q1c.m, ~ FunGroup), adjust = "BH")

par(mfrow = c(1,1))


#### Figure 2: Cor Plots ####
mcl.df.q.trait <- merge(mcl.df.q[,-13], seed.traits, by.x = "Species_Name", by.y = "Species", all = F)

##### Spatial ####
colonization <- c("c", "new.mass",  "shape", "set.time.mpsec", "height.cm", "size.mm", "ldd.natural", "wing.loading")

trait.c <- data.frame(trait = character(), p.value = numeric(), cor = numeric(), FunGroup = character())

for(i in colonization[-1]){
  for(j in unique(mcl.df.q.trait$FunGroup)) {
    
    tmp <- cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$c, mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i])
    
    #tmp2 <- SpearmanRho(x = mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$c, 
    #                   y = mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i], 
    #                   conf.level = 0.95)
    
    trait.c <- rbind(trait.c, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = j, ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
    
    #trait.c <- rbind(trait.c, data.frame(trait = i, p.value = tmp$p.value, cor = tmp2[1], FunGroup = j, ci.lower = tmp2[2], ci.upper = tmp2[3]))
  }
    tmp <- cor.test(mcl.df.q.trait$c, mcl.df.q.trait[,i])
    
  trait.c <- rbind(trait.c, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = "all", ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
}  

#trait.c$trait <- factor(trait.c$trait, levels = trait.c[trait.c$FunGroup == "all",]$trait[order(trait.c[trait.c$FunGroup == "all",]$cor)])

trait.c$trait <- factor(trait.c$trait, levels = trait.c[trait.c$FunGroup == "all",]$trait[order(abs(trait.c[trait.c$FunGroup == "all",]$cor))])

trait.c$FunGroup <- factor(trait.c$FunGroup, levels = c("non-native grass", "non-native forb", "native forb", "all"))

c <- ggplot(trait.c, aes(x = cor, y = trait, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_rect(mapping = aes(xmin = -1, xmax = 1, 
  #                         ymin = 2.57, ymax = 3.45), col = "grey86", fill =    #"grey86") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7), aes(fill = FunGroup)) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper), width = 0.1, position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.line = element_blank(),
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(linewidth = 1, fill = NA)
  ) +
  scale_x_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1), labels = c(-1,-0.5,0,0.5,1)) +
  scale_shape_manual(values = c(15, 17, 19, 23)) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  #scale_y_discrete(labels = c("wing\n loading","settling\n speed", "mass","height", "shape", "size", "dispersal\n potential")) + #old order
  #scale_y_discrete(labels = c("mass", "height", "settling\n speed",  "size", "wing\n loading", "shape",  "dispersal\n potential")) + 
  labs(x = "correlation", title = "Spatial Dispersal")  


##### Temporal ####
survival <- c("s", "new.mass",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.width")

trait.s <- data.frame(trait = character(), p.value = numeric(), cor = numeric(), FunGroup = character(), ci.lower = numeric(), ci.upper = numeric())

for(i in survival[-1]){
  for(j in unique(mcl.df.q.trait$FunGroup)) {
      
    tmp <- cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$s, mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i])
    
    # tmp2 <- SpearmanRho(x = mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,]$s, 
    #                    y = mcl.df.q.trait[mcl.df.q.trait$FunGroup == j,i], 
    #                    conf.level = 0.95)
    
    trait.s <- rbind(trait.s, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = j, ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
    
    #trait.s <- rbind(trait.s, data.frame(trait = i, p.value = tmp$p.value, cor = tmp2[1], FunGroup = j, ci.lower = tmp2[2], ci.upper = tmp2[3]))
  }
  tmp <- cor.test(mcl.df.q.trait$s, mcl.df.q.trait[,i])
    
  trait.s <- rbind(trait.s, data.frame(trait = i, p.value = tmp$p.value, cor = tmp$estimate, FunGroup = "all", ci.lower = tmp$conf.int[1], ci.upper = tmp$conf.int[2]))
}  

trait.s$FunGroup <- factor(trait.s$FunGroup, levels = c("non-native grass", "non-native forb", "native forb", "all"))

trait.s$trait <- factor(trait.s$trait, levels = trait.s[trait.s$FunGroup == "all",]$trait[order(abs(trait.s[trait.s$FunGroup == "all",]$cor))])
 
s <- ggplot(trait.s, aes(x = cor, y = trait, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  # geom_rect(mapping = aes(xmin = -1, xmax = 1, 
  #                         ymin = 3.55, ymax = 4.45), col = "grey86", fill = "grey86") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(size = 2.5, position = position_dodge(width = 0.7), aes(fill = FunGroup)) +
  geom_errorbar(aes(xmin = ci.lower, xmax = ci.upper), width = 0.1, position = position_dodge(width = 0.7)) +
  theme_classic() +
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.line = element_blank(),
    axis.title.x = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(linewidth = 1, fill = NA)
  ) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_shape_manual(values = c(15, 17, 19, 23)) +
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3","black")) +
  #scale_y_discrete(labels = c("size", "shape", "%C", "mass", "coat\n permeability", "E:S", "%N", "coat\n thickness")) +
  #scale_y_discrete(labels = c("mass", "%C", "coat\n thickness", "coat\n permeability", "size", "%N", "shape")) +
  labs(x = "correlation", title = "Temporal Dispersal")  

cor.plot <- ggarrange(s,c, widths = c(0.7,1), labels = c("(a)", "(b)"))
  
ggsave("Manuscript/McLaughlin/Figures/cor-plot.png", cor.plot, height = 6, width = 9, units = "in", dpi = 600)


#### Figure 3: HMM & traits  ####

##### Shape ####
a <- ggplot(mcl.df.q.trait, aes(x = shape, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", aes(col = FunGroup), se = F) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = shape, y = c), method = "lm", col = "black", se = F, linetype = 2, size = 1.3) +  
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.background = element_rect(color = "white")
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.038,0.28)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal")
 
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$shape)  #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$shape) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$shape) #ns

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$shape) #marg


f <- ggplot(mcl.df.q.trait, aes(x = shape, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", aes(col = FunGroup), se = F) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = shape, y = s), method = "lm", col = "black", se = F, size = 1.3) +  
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "temporal dispersal")
 
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$shape) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$shape) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$shape) #ns

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$shape)

##### Size ####
b <- ggplot(mcl.df.q.trait, aes(x = log(size.mm), y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", aes(col = FunGroup), linetype = 2, se = F) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = log(size.mm), y = s), method = "lm", col = "black", se = F, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "log size (mm)")
 
  
cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$size.mm) #marg

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$size.mm) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$size.mm) #ns

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$size.mm)

g <- ggplot(mcl.df.q.trait, aes(x = log(size.mm), y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "non-native grass",], method = "lm", aes(col = FunGroup), se = F, linetype = 2) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = log(size.mm), y = c), method = "lm", col = "black", se = F, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),    
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.038,0.28)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_linetype_manual(values = c(1,2)) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "log size (mm)")

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$size.mm) #marg

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$size.mm) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$size.mm) #sig

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$size.mm)

##### Coat thick ####

c <- ggplot(mcl.df.q.trait, aes(x = coat.thick.per.width, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",], method = "lm", fill = "#d95f02", se = F) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = coat.thick.per.width, y = s), method = "lm", col = "black", se = F, linetype = 2, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "seed coat to width ratio")

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$coat.thick.per.width) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$coat.thick.per.width) #ns

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$coat.thick.per.width) #marg

# ggplot(mcl.df.q.trait, aes(x = both.thick, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
#   geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
#   #geom_text(aes(label = Code)) +
#   geom_point(size = 2) +
#  theme_classic() +
#   theme(
#     axis.text = element_text(size = 12),
#     axis.title = element_text(size = 14),
#     legend.position = "none",
#     legend.title = element_blank(),
#     legend.text = element_text(size = 12)
#   ) +
#   #scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.25)) +
#   scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
#   scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
#   scale_shape_manual(values = c( 19, 17, 15)) +
#   labs(y = "temporal dispersal")
  
  
##### N ####
d <- ggplot(mcl.df.q.trait, aes(x = prop.N, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup != "native forb",], method = "lm", aes(fill = FunGroup), se = F) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = prop.N, y = s), method = "lm", col = "black", se = F, linetype = 2, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "%N")

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$prop.N) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$prop.N) #ns

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$prop.N) #sig

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$prop.N)

##### Disp ####
e <- ggplot(mcl.df.q.trait, aes(x = ldd.natural, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", fill = "#7570b3", se = F) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = ldd.natural, y = c), method = "lm", col = "black", se = F, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.038,0.28)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "dispersal potential")


cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$ldd.natural) #ns


cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$ldd.natural) #sig

##### Wing loading ####
o <- ggplot(mcl.df.q.trait, aes(x = log(wing.loading), y = c, col = FunGroup, shape = FunGroup)) +
  geom_smooth(data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",], method = "lm", fill = "#7570b3", se = F) + #alpha = 0.2) +
  geom_smooth(inherit.aes = F, aes(x = log(wing.loading), y = c), method = "lm", col = "black", se = F, size = 1.3) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
  theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  #scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "log wing loading")

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",]$wing.loading) #ns


cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",]$wing.loading) #sig

cor.test(mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$c, 
         mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",]$wing.loading) #ns

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$wing.loading) #sig

##### Panel Plot ####
# leg <- get_legend(
#   ggplot(mcl.df.q.trait, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
#     theme_classic() + 
#     geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2, se = F) +
#     theme(
#       legend.title = element_blank(),
#       legend.text = element_text(size = 14)
#     ) +
#     geom_point(size = 2) + 
#     scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
#     scale_color_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
#     scale_shape_manual(values = c( 19, 17, 15))
#   )
# 
# trait.panel <- ggarrange(a,g,e,o,f,b,c,d, ncol = 4, nrow = 2, labels = c("(a)", "(b)", "(c)" , "(d)", "(e)", "(f)", "(g)", "(h)"), widths = c(1, 0.85, 0.85, 0.85), heights = c(1,1), label.x = c(0.2,0.04,0.04,0.04,0.2,0.04,0.04,0.04), label.y = 0.96, common.legend = TRUE,legend="bottom")

trait.panel.s <- ggarrange(a,g,e,o, ncol = 4, nrow = 1, labels = c("(a)", "(b)", "(c)" , "(d)"), widths = c(1, 0.85, 0.85, 0.85), heights = 1, label.x = c(0.2,0.04,0.04,0.04), label.y = 0.96) + 
  bgcolor("white") +
  border(color = "white")

trait.panel.s <- annotate_figure(trait.panel.s, top = text_grob("Spatial Dispersal", 
               color = "black", face = "bold", size = 16)) + 
  bgcolor("white") +
  border(color = "white")

trait.panel.t <- ggarrange(f,b,c,d, ncol = 4, nrow = 1, labels = c("(e)", "(f)", "(g)", "(h)"), widths = c(1, 0.85, 0.85, 0.85), heights = c(1,1), label.x = c(0.2,0.04,0.04,0.04), label.y = 0.96, common.legend = TRUE, legend = "bottom") + 
  bgcolor("white") +
  border(color = "white")

trait.panel.t <- annotate_figure(trait.panel.t, top = text_grob("Temporal Dispersal", color = "black", face = "bold", size = 16)) + 
  bgcolor("white") +
  border(color = "white")

trait.panel <- ggarrange(trait.panel.s, trait.panel.t, ncol = 1, nrow = 2, heights = c(0.9, 1))
                           
ggsave("Manuscript/McLaughlin/Figures/trait-panel.png", trait.panel, height = 7, width = 11, units = "in", dpi = 600)

#### Linear models ####
par(mfrow = c(2,2))

##### Shape ####
m.shape.s <- lm(s ~ shape, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.shape.s)
summary(m.shape.s)
  
m.shape.c <- lm(c ~ shape, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.shape.c)
summary(m.shape.c)

##### Size ####
m.size.s <- lm(s ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.size.s)
summary(m.size.s)
  
m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.size.c)
summary(m.size.c)

m.size.c <- lm(c ~ size.mm, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",])
plot(m.size.c)
summary(m.size.c)

##### Coat ####
m.coat.s <- lm(s ~ coat.thick.per.width, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",])
plot(m.coat.s)
summary(m.coat.s)

##### N ####
m.N.c.EF <- lm(s ~ prop.N, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native forb",])
plot(m.N.c.EF)
summary(m.N.c.EF)

m.N.c.EG <- lm(s ~ prop.N, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "non-native grass",])
plot(m.N.c.EG)
summary(m.N.c.EG)

1##### Disp ####
m.disp.c <- lm(c ~ ldd.natural, data = mcl.df.q.trait[mcl.df.q.trait$FunGroup == "native forb",])
plot(m.disp.c)
summary(m.disp.c)

par(mfrow = c(1,1))

#### Figure S1 ####

##### Mass ####
h <- ggplot(mcl.df.q.trait, aes(x = new.mass, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(y = "spatial dispersal", x = "log mass (mg)")

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$new.mass)

i <- ggplot(mcl.df.q.trait, aes(x = chem.mass.mg, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  geom_text(aes(label = Code)) +
  #geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  labs(y = "temporal dispersal", x = "log mass (mg)")

ggplot(mcl.df.q.trait, aes(x = new.mass, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) + 
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  #scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  labs(y = "temporal dispersal", x = "log mass (mg)")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$new.mass)

cor.test(mcl.df.q.trait[mcl.df.q.trait$group == "grass",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$group == "grass",]$chem.mass.mg)

m1 <- lm(s ~ morph.mass.mg, mcl.df.q.trait[mcl.df.q.trait$group == "grass",])
plot(m1)

cor.test(mcl.df.q.trait[mcl.df.q.trait$group == "grass" & mcl.df.q.trait$new.code != "GASVEN",]$s, 
         mcl.df.q.trait[mcl.df.q.trait$group == "grass" & mcl.df.q.trait$new.code != "GASVEN",]$new.mass)


library(caret)


#fit a regression model and use LOOCV to evaluate performance
model <- train(s ~ new.mass, data = mcl.df.q.trait[mcl.df.q.trait$group == "grass",], method = "lm", trControl = trainControl(method = "LOOCV"))

print(model)

model <- train(s ~ new.mass, data = mcl.df.q.trait[mcl.df.q.trait$group == "grass" & mcl.df.q.trait$new.code != "GASVEN",], method = "lm", trControl = trainControl(method = "LOOCV"))

print(model)

##### Set speed ####
j <- ggplot(mcl.df.q.trait, aes(x = set.time.mpsec, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "settling speed (m/s)")

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$set.time.mpsec) #ns

##### C ####
k <- ggplot(mcl.df.q.trait, aes(x = prop.C, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_x_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "%C")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$prop.C) #ns

##### Coat perm ####
l <- ggplot(mcl.df.q.trait, aes(x = coat.perm.perc, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "log coat permeability")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$coat.perm.perc) #ns

##### Height ####
m <- ggplot(mcl.df.q.trait, aes(x = height.cm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2), limits = c(0.01,0.3)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "height (cm)")

cor.test(mcl.df.q.trait$c, 
         mcl.df.q.trait$height.cm) #ns

##### ES ####
n <- ggplot(mcl.df.q.trait, aes(x = E.S, y = s, col = FunGroup, group = FunGroup, shape = FunGroup)) +
  #geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2) +
  #geom_text(aes(label = Code)) +
  geom_point(size = 2) +
 theme_classic() +
  theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(breaks = c(0.2, 0.4,0.6,0.8), limits = c(0.2,0.85)) +
  scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
  scale_color_manual(values = c("#7570b3", "#d95f02", "#1b9e77")) +
  scale_shape_manual(values = c( 19, 17, 15)) +
  labs(x = "embryo:seed")

cor.test(mcl.df.q.trait$s, 
         mcl.df.q.trait$coat.perm.perc)

##### Panel Plot ####
leg <- get_legend(
  ggplot(mcl.df.q.trait, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup, shape = FunGroup)) +
    theme_classic() + 
    geom_smooth(method = "lm", aes(fill = FunGroup), alpha = 0.2, se = F) +
    theme(
      legend.title = element_blank(),
      legend.text = element_text(size = 12)
    ) +
    geom_point(size = 2) + 
    scale_fill_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
    scale_color_manual(values = c("#7570b3","#d95f02", "#1b9e77")) +
    scale_shape_manual(values = c( 19, 17, 15))
  )

s3.appen.panel <- ggarrange(h,j,m,i,k,l, ncol = 3, nrow = 2, labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), widths = c(1, 0.85, 0.85), heights = c(1,1,1), label.x = c(0.2,0.04,0.04,0.2,0.04,0.04), label.y = 0.96, common.legend = T, legend = "bottom") +
  bgcolor("white") +
  border(color = "white")

ggsave("Manuscript/McLaughlin/Figures/s3-appen-panel.jpeg", s3.appen.panel, height = 5.5, width = 8, units = "in", dpi = 600)

#### Correlation check ####
mcl.sim <- readRDS("Scripts/Validation/McL-sim-154.RDS")

mcl.sim <- as.data.frame(mcl.sim)

mcl.sim <- filter(mcl.sim, iter < 150)

sim.df <- data.frame(cor = NA, p.value = NA)

for(i in 1:50000){
  tmp <- slice_sample(mcl.sim, n = 59, replace = F)
  tmp.test <- cor.test(tmp$s, tmp$c)
  sim.df[i,]$cor <- tmp.test$estimate
  sim.df[i,]$p.value <- tmp.test$p.value
}


ggplot(sim.df, aes(x = cor)) +
  geom_histogram(bins = 20, col = "black", fill = "white") +
  theme_classic() +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 1),
    axis.line = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  geom_vline(xintercept = -0.553, col = "red") +
  labs(x = "simulated correlation coefficient")

mean(sim.df$cor)
nrow(sim.df[sim.df$cor <= -0.553,])/nrow(sim.df) #p value < 0.0001

#### model validation ####
#mcl.sim.sp <- readRDS("Scripts/Validation/McL-sim-species-validation.RDS")
mcl.sim.sp <- readRDS("Scripts/Validation/McL-sim-species-validation-50.RDS")
mcl.sim.sp <- as.data.frame(mcl.sim.sp)

mcl.sim.sp <- mcl.sim.sp %>% 
  dplyr::mutate(across(p0:iter.sd, ~as.numeric(.x)))

colnames(mcl.sim.sp)[2:7] <- paste(colnames(mcl.sim.sp)[2:7], "est", sep=".")

mcl.sim.sp <- merge(mcl.sim.sp, mcl.df.q, by = "Species_Name")

#mcl.sim.sp <- filter(mcl.sim.sp, mean.years > 3.1)

a <- ggplot(mcl.sim.sp, aes(y = s.est, x = s)) +
  #geom_text(aes(label = Code)) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point() +
  theme_classic() +
    theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  labs(x = "real s", y = "mean estimated s")
  

b <- ggplot(mcl.sim.sp, aes(y = c.est, x = c)) +
  #geom_text(aes(label = Code)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point() +
  theme_classic() +
    theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  labs(x = "real c", y = "mean estimated c")

c <- ggplot(mcl.sim.sp, aes(y = g, x = g.est)) +
  #geom_text(aes(label = Code)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  geom_point() +
  theme_classic() +
    theme(
    panel.border = element_rect(linewidth = 1, fill = NA),
    axis.line = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 12),
  ) +
  labs(x = "real g", y = "mean estimated g")

ggarrange(a,b,c, ncol = 3, nrow = 1, labels = c("(a)", "(b)", "(c)"))

ggplot(mcl.sim.sp, aes(y = n.plots, x = s.sd)) +
  geom_text(aes(label = Code)) 

ggplot(mcl.sim.sp, aes(y = n.plots, x = c.sd)) +
  geom_text(aes(label = Code)) 

ggplot(mcl.sim.sp, aes(y = n.plots, x = g.sd)) +
  geom_text(aes(label = Code)) 

ggplot(mcl.sim.sp, aes(y = mean.years, x = s.sd)) +
  geom_text(aes(label = Code)) 

ggplot(mcl.sim.sp, aes(y = mean.years, x = c.sd)) +
  geom_text(aes(label = Code)) 

ggplot(mcl.sim.sp, aes(y = n.plots, x = g.sd)) +
  geom_text(aes(label = Code)) 

ggplot(mcl.sim.sp, aes(y = mean.years, x = iter.est)) +
  geom_text(aes(label = Code)) 

ggplot(mcl.sim.sp, aes(y = mean.years, x = iter.sd)) +
  geom_text(aes(label = Code)) 

# to me this confirms a few things: trifolium microdon and microcephalum were confused, as were EROBOT and EROBRA; also problem species with high standard deviations: Anagallis, Lepidium, Gastridium, 

#### Model boots ####
boot <- readRDS("Scripts/Validation/McL-sim-species-boot.RDS")
boot <- as.data.frame(boot)

boot <- boot %>% 
  dplyr::mutate(across(p0:iter.sd, ~as.numeric(.x)))

colnames(boot)[2:7] <- paste(colnames(boot)[2:7], "mean", sep=".")

boot <- merge(boot, mcl.df.q, by = "Species_Name")


ggplot(boot, aes(y = s.mean, x = iter.mean)) +
  geom_point() +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(ymin = s.mean - s.sd, ymax = s.mean + s.sd))

ggplot(boot, aes(y = c.mean, x = n.plots, col = mean.years)) +
  geom_point() +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(ymin = c.mean - c.sd, ymax = c.mean + c.sd)) #LEPNIT and APHOCC look terrible

ggplot(boot, aes(y = g.mean, x = iter.mean)) +
  geom_point() +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(ymin = g.mean - g.sd, ymax = g.mean + g.sd))

ggplot(boot, aes(y = iter.mean, x = Code)) +
  geom_point() +
  geom_text(aes(label = Code)) +
  geom_errorbar(aes(ymin = iter.mean - iter.sd, ymax = iter.mean + iter.sd))

#### seed survival ####
surv <- read.csv("Data/Seed-bag-survival_Year1.csv")
surv$s.perc <- surv$n.viable/100

surv.sum <- surv %>%
  group_by(Species) %>%
  summarize(mean.s = mean(s.perc))

surv.sum <- merge(surv.sum, mcl.df.q[,c(1,5)], by.x = "Species", by.y = "Species_Name", all.y = F, all.x = F)

ggplot(surv.sum, aes(y = mean.s, x = s)) +
  geom_point() +
  geom_smooth(method = "lm")
