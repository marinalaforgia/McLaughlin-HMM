# Seed traits at McLaughlin Analysis #
rm(list = ls())

source("/Users/Marina/Documents/UC-Davis/Quantitative-Resources/Scripts/Better_Pairs.R")

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

adult.traits <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Modified_CombinedFiles/McL_80SitesSpeciesTraits_012615.csv") # adult traits
mcl <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Data-Storage-Files/Core_Community_Data2019.csv") # community data
serp <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Data-Storage-Files/McL_Abiotic-Data.csv") # serpentine data
mcl.df.t <- readRDS("HMMs/McLaughlin/20221212_mcl_HMM_transect.RDS") # transect level HMMs
mcl.df.q <- readRDS("HMMs/McLaughlin/20221212_mcl_HMM_quadrat.RDS") # quadrat level HMM
seed.species <- read.csv("Data/20211001_Full-Species-List.csv")
seed.traits <- read.csv("Data/20230328_Seed-Traits_cleaning.csv")
calipc <- read.csv("Collections/cal-ipc_inventory.csv")


rm <- c("Logfia spp./Micropus californicus", "Aegilops triuncialis", "Trifolium bifidum/gracilentum", "Torilis arvensis/nodosa", "Hesperolinon spp.", "Lactuca serriola", "Cryptantha hispidula", "Cuscuta californica", "Silene gallica", "Festuca octoflora") # rm commonly mis-ID'd species, species that hadnt always been censused (cuscuta) and goatgrass which is treated at the reserve

mcl <- filter(mcl, Species_Name %in% adult.traits[adult.traits$Annual.Perennial == "Annual",]$Species_Name, !(Species_Name %in% rm)) 
mcl <- merge(serp[,1:2], mcl, by = "Site")

mcl[mcl$Species_Name == "Geranium molle",]$Species_Name <- "Geranium dissectum" # Susan confirmed these should all be dissectum

mcl <- mcl %>%
  dplyr::group_by(Year, Site, Quadrat, Serpentine, Species_Name) %>%
  dplyr::summarize(Cover = sum(Cover)) # re-summarize by species after merging geraniums

mcl.df.q$Species_Name <- recode_factor(mcl.df.q$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis")

mcl$Species_Name <- recode_factor(mcl$Species_Name, 'Lysimachia arvensis' = "Anagallis arvensis")

#mcl.df.q <- merge(mcl.df.q, calipc[,c(1,3)], by.x = "Species_Name", by.y = "Scientific.name", all.x = T, all.y = F)
rm(serp, adult.traits, rm)

#### Seed trait prep ####
#rm <- c("wing.loading.c", "area.mm2.m", "area.mm2.c")
#seed.traits <- seed.traits[seed.traits$Species != "Festuca octoflora", !colnames(seed.traits) %in% rm]

seed.traits <- filter(seed.traits, Species %in% mcl.df.q$Species_Name)
seed.traits[is.na(seed.traits$prop.C),]$prop.C <- mean(seed.traits$prop.C, na.rm = T)
seed.traits[is.na(seed.traits$prop.N),]$prop.N <- mean(seed.traits$prop.N, na.rm = T)
seed.traits[is.na(seed.traits$cn),]$cn <- mean(seed.traits$cn, na.rm = T)

seed.traits$coat.thick.per.size <- seed.traits$both.thick/seed.traits$size.mm
seed.traits$coat.thick.per.mass <- seed.traits$both.thick/seed.traits$chem.mass.mg

# all species
hist(scale(log(seed.traits$wing.loading)))
seed.traits$wing.loading <- scale(log(seed.traits$wing.loading))

hist(scale(log(seed.traits$cn)))
seed.traits$cn <- scale(log(seed.traits$cn))

hist(scale(seed.traits$prop.C))
seed.traits$prop.C <- scale(scale(seed.traits$prop.C))

hist(scale(log(seed.traits$prop.N)))
seed.traits$prop.N <- scale(log(seed.traits$prop.N))

hist(scale(log(seed.traits$coat.perm.perc)))
seed.traits$coat.perm.perc <- scale(log(seed.traits$coat.perm.perc))


hist(scale(log(seed.traits$both.thick/seed.traits$size.mm)))
seed.traits$coat.thick.per.size <- scale(log(seed.traits$both.thick/seed.traits$size.mm))

hist(scale(log(seed.traits$both.thick/seed.traits$chem.mass.mg)))
seed.traits$coat.thick.per.mass <- scale(log(seed.traits$both.thick/seed.traits$chem.mass.mg))


hist(scale(log(seed.traits$morph.mass.mg)))
seed.traits$morph.mass.mg <- scale(log(seed.traits$morph.mass.mg))

hist(log(seed.traits$chem.mass.mg))
seed.traits$chem.mass.mg <- log(seed.traits$chem.mass.mg)

hist(scale(log(seed.traits$size.mm)))
seed.traits$size.mm <- scale(log(seed.traits$size.mm))

hist(seed.traits$set.time.mpsec)
seed.traits$set.time.mpsec <- scale(seed.traits$set.time.mpsec)

hist(sqrt(seed.traits$shape))
seed.traits$shape <- scale(sqrt(seed.traits$shape))

hist(scale(log(seed.traits$both.thick)))
seed.traits$both.thick <- scale(log(seed.traits$both.thick))

hist(scale(log(seed.traits$height.cm)))
seed.traits$height.cm <-  scale(log(seed.traits$height.cm))
  
#hist(scale(log(seed.traits$E.S)))
#seed.traits$E.S <- log(seed.traits$E.S)



# ## forbs all the same
# forb.traits <- filter(seed.traits, group == "forb")
# 
# hist(log(forb.traits$wing.loading))
# forb.traits$wing.loading <- log(forb.traits$wing.loading)
# 
# hist(log(forb.traits$cn))
# forb.traits$cn <- log(forb.traits$cn)
# 
# hist(forb.traits$prop.C)
# 
# hist(log(forb.traits$prop.N))
# forb.traits$prop.N <- log(forb.traits$prop.N)
# 
# hist(log(forb.traits$coat.perm.perc))
# forb.traits$coat.perm.perc <- log(forb.traits$coat.perm.perc)
# 
# 
# hist(log(forb.traits$both.thick/forb.traits$size.mm))
# forb.traits$coat.thick.per.size <- log(forb.traits$both.thick/forb.traits$size.mm)
# 
# hist(log(forb.traits$both.thick/forb.traits$chem.mass.mg))
# forb.traits$coat.thick.per.mass <- log(forb.traits$both.thick/forb.traits$chem.mass.mg)
#forb.traits$coat.thick.per.mass <- forb.traits$both.thick/forb.traits$chem.mass.mg
# 
# hist(log(forb.traits$morph.mass.mg))
# forb.traits$morph.mass.mg <- log(forb.traits$morph.mass.mg)
# 
# hist(log(forb.traits$chem.mass.mg))
# forb.traits$chem.mass.mg <- log(forb.traits$chem.mass.mg)
# 
# hist(log(forb.traits$size.mm))
# forb.traits$size.mm <- log(forb.traits$size.mm)
# 
# hist(forb.traits$set.time.mpsec)
# 
# hist(log(forb.traits$E.S))
# forb.traits$E.S <- log(forb.traits$E.S)

seed.traits <- merge(seed.traits, calipc[,c(1,3)], by.x = "Species", by.y = "Scientific.name", all.x = T, all.y = F)

seed.traits$FunGroup <- paste(seed.traits$nat.inv, seed.traits$group, sep = " ")
seed.traits[is.na(seed.traits$Rating),]$Rating <- "None"
seed.traits[seed.traits$nat.inv == "native",]$Rating <- "Native"


trait <- c("Species", "code", "family", "group", "nat.inv", "FunGroup", "Rating", 
           "chem.mass.mg",  "shape", "size.mm", 
           "appendage.type", "appendage", "prop.C", "prop.N", "cn", "set.time.mpsec", "height.cm", "ldd2", "coat.perm.perc", "E.S", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.size", "coat.thick.per.mass", "wing.loading")

trait.pca <- c("chem.mass.mg",  "shape", "size.mm", 
"prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "E.S", "both.thick")

#### PCA ALL species ####

tmp.pca <- seed.traits

pca <- prcomp(tmp.pca[, trait.pca], scale = T)
screeplot(pca)
summary(pca)


all.pca <- cbind(tmp.pca, pca$x[,1:4])

# color points by c/s gradient 
# color/size by relative abundance 

autoplot(pca, x = 1, y = 2, data = all.pca, frame = F, loadings = T, loadings.label = T, label = F, col = "FunGroup", size = 2) +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) #+
  #scale_color_manual(values = c(adhesive = "#1B9E77", ant = "red4", unassisted = "darkgoldenrod3", ingestion = "#7570B3", wind = "#F17236"))

autoplot(pca, x = 3, y = 4, data = all.pca, frame = F, loadings = T, loadings.label = T, label = T, col = "FunGroup", size = 2) +
  theme_classic() +
  #geom_text(aes(label = code, col = Rating)) +
  #stat_ellipse(aes(group = group)) + 
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) 

#### PCA Forbs ####
trait <- c("Species", "code", "family", "group", "nat.inv", "FunGroup", "Rating", "morph.mass.mg", "chem.mass.mg",  "shape", "size.mm", "appendage.type", "appendage", "prop.C", "prop.N", "cn", "set.time.mpsec", "height.cm", "ldd2", "coat.perm.perc", "E.S", "both.thick", "mucilage", "fleshy.end", "starchy.end", "coat.thick.per.size", "coat.thick.per.mass", "ldd")

trait.pca <- c("chem.mass.mg",  "shape", "size.mm", 
"prop.C", "prop.N", "set.time.mpsec", "height.cm", "coat.perm.perc", "E.S", "coat.thick.per.size", "ldd")

trait.pca.forb <- seed.traits[seed.traits$group == "forb", trait]
pca.forb <- prcomp(trait.pca.forb[, trait.pca], scale = T)
summary(pca.forb)
screeplot(pca.forb)
trait.pca.forb <- cbind(trait.pca.forb, pca.forb$x[,1:4])


autoplot(pca.forb, x = 1, y = 2, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "nat.inv") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
  ) 

autoplot(pca.forb, x = 3, y = 4, data = trait.pca.forb, frame = F, loadings = T, loadings.label = T, label = F, col = "nat.inv") +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.title = element_blank()
   ) 




#### HMMs & Traits (all) ####
##### Transect Plots ####
# mcl.df.t.all <- merge(mcl.df.t, trait.pca, by.x = "Species_Name", by.y = "Species", all.x = T)
# 
# ggplot(mcl.df.t.all[mcl.df.t.all$iter < 150 & mcl.df.t.all$n.plots >= 16,], aes(x = s, y = c, col = PC1)) +
#   geom_smooth(method = "lm", col = "black") + 
#   #geom_text(aes(label = Code)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_viridis_c(na.value = "lightgrey")
# 
# ggplot(mcl.df.t.all[mcl.df.t.all$iter < 150 & mcl.df.t.all$n.plots >= 16,], aes(x = s, y = c, col = PC2)) +
#   geom_smooth(method = "lm", col = "black") + 
#   #geom_text(aes(label = Code)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_viridis_c(na.value = "lightgrey")
# 
# ggplot(mcl.df.t.all[mcl.df.t.all$iter < 150 & mcl.df.t.all$n.plots >= 16,], aes(x = s, y = c, col = PC3)) +
#   geom_smooth(method = "lm", col = "black") + 
#   #geom_text(aes(label = Code)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_viridis_c(na.value = "lightgrey")

##### Quadrats Plot ####
mcl.df.q.all <- merge(mcl.df.q, all.pca[,-34], by.x = "Species_Name", by.y = "Species", all.x = T)

# no relationship between avg presence and c or s estimate
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = avg.P.sp)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey") # no relationship

# negative relationship between avg presence and PC1, so that species that are more common are typically larger seeded (i.e. grasses... duh)
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = avg.P.sp, y = PC1)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

# no relationship
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = avg.P.sp, y = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c)) +
  geom_smooth(method = "lm", col = "black") +
  geom_text(aes(label = Code)) +
  #geom_point() +
  theme_bw() 

## PC1 ##
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80 & mcl.df.q.all$FunGroup != "Exotic Grass",], aes(x = s, y = c, col = PC1)) +
  geom_smooth(method = "lm", col = "black") +
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")
  #scale_color_gradient(high = "red", low = "blue")


ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC1)) +
  geom_smooth(method = "lm", col = "black") +
  #geom_text(aes(label = Code)) +
  geom_point() +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")
  #scale_color_gradient(high = "red", low = "blue")

summary(lm(c ~ PC1, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))

summary(lm(s ~ PC1, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))

## PC2 ##
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  #geom_text(aes(label = code)) + 
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")
  #scale_color_gradient(high = "red", low = "blue")

ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  #geom_text(aes(label = code)) + 
  theme_bw() +
  facet_wrap(~FunGroup) +
  scale_color_viridis_c(na.value = "lightgrey")
  #scale_color_gradient(high = "red", low = "blue")

## PC3 ##
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC3)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC3)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

## PC3 ##
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC4)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC4)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_text(aes(label = Code)) +
  geom_point() +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

# ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = FunGroup)) +
#   geom_smooth(method = "lm", col = "black") + 
#   #geom_text(aes(label = Code)) +
#   geom_point() +
#   theme_bw()



# ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = PC3)) +
#   geom_smooth(method = "lm", col = "black") + 
#   #geom_text(aes(label = Code)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_viridis_c(na.value = "lightgrey")


ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = FunGroup)) +
  geom_point(aes(shape = FunGroup), size = 2) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  theme_bw() +
  scale_color_viridis_d()


##### Quadrat Pairs ####
rm <- c("ACMBRA", "CERGLO", "LEPBIC", "CROSET", "GALAPA" ,"EROBRA", "GASPHL")
col.rm <- c("appendage.type", "appendage", "code", "disp", "ldd2", "Rating", "FunGroup.y")
mcl.df.q.c.all <- filter(mcl.df.q.all[,!colnames(mcl.df.q.all) %in% col.rm], !Code %in% rm, iter < 150, n.plots >= 80)

# colonization 
# + size
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = size.mm)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  facet_wrap(~FunGroup) +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

# -N
# -wing loading
# - E.S
# - coat thick per size
ggplot(mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,], aes(x = s, y = c, col = coat.thick.per.size)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  #facet_wrap(~FunGroup) +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

#survival
# + settling time
# - height
# - size
# - shape
# + N
# + wing loading
# + E.S
# + coat thick
# + coat thick per size




pairs(mcl.df.q.c.all[, c(3:5,20:25)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(mcl.df.q.c.all[, c(3:5,26:30)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(mcl.df.q.c.all[, c(3:5,31:36)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(mcl.df.q.c.all[, c(3:5,37:41)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(mcl.df.q.c.all[, c(3:5,42:46)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

##### Transect Pairs ####
# 
# # survival
# ## -shape
# ## +settling time
# ## -CN
# ## -PC2
# # colonization 
# ## +shape
# ## -E:S
# ## +PC2
# 
# mcl.df.t.c.all <- mcl.df.t.all[mcl.df.t.all$iter < 150 & mcl.df.t.all$n.plots >= 16,]
# mcl.df.t.c.all <- mcl.df.t.c.all[complete.cases(mcl.df.t.c.all),]
# 
# pairs(mcl.df.t.c.all[, c(4:5,25:30)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)
# 
# pairs(mcl.df.t.c.all[, c(4:5,31:35)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)
# 
# pairs(mcl.df.t.c.all[, c(4:5,36:39)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)

##### All models ####
m1 <- lm(s ~ FunGroup, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80 & mcl.df.q.all$FunGroup != "Native Grass",])

pairs(emmeans(m1, ~ FunGroup), adjust = "BH")

m2 <- lm(c ~ FunGroup, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80 & mcl.df.q.all$FunGroup != "Native Grass",])

pairs(emmeans(m2, ~ FunGroup), adjust = "BH")

# main model
summary(lm(s ~ c, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))

## by PC ##
summary(lm(s ~ PC1, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))
summary(lm(c ~ PC1, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))

# summary(lm(s ~ PC2, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))
# summary(lm(c ~ PC2, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))

summary(lm(s ~ PC3, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))
#summary(lm(c ~ PC3, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))

# summary(lm(s ~ PC4, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))
# summary(lm(c ~ PC4, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))


## by traits ##
summary(lm(s ~ coat.thick.per.size, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))
summary(lm(c ~ coat.thick.per.size, data = mcl.df.q.all[mcl.df.q.all$iter < 150 & mcl.df.q.all$n.plots >= 80,]))

#### HMMs & Traits (forbs) ####
##### Transect Plots ####
# mcl.df.t.forb <- filter(mcl.df.t, FunGroup == "Native Forb" | FunGroup == "Exotic Forb")
# mcl.df.t.forb <- merge(mcl.df.t.forb, trait.pca.forb, by.x = "Species_Name", by.y = "Species", all.x = T)
# 
# ggplot(mcl.df.t.forb[mcl.df.t.forb$iter < 150 & mcl.df.t.forb$n.plots > 8,], aes(x = s, y = c, col = PC1)) +
#   geom_smooth(method = "lm", col = "black") + 
#   geom_text(aes(label = Code)) +
#   #geom_point() +
#   theme_bw() +
#   scale_color_viridis_c(na.value = "lightgrey")
# 
# ggplot(mcl.df.t.forb[mcl.df.t.forb$iter < 150 & mcl.df.t.forb$n.plots > 8,], aes(x = s, y = c, col = PC2)) +
#   geom_smooth(method = "lm", col = "black") + 
#   geom_text(aes(label = Code)) +
#   #geom_point() +
#   theme_bw() +
#   scale_color_viridis_c(na.value = "lightgrey")
# 
# ggplot(mcl.df.t.forb[mcl.df.t.forb$iter < 150 & mcl.df.t.forb$n.plots > 8,], aes(x = s, y = c, col = PC3)) +
#   geom_smooth(method = "lm", col = "black") + 
#   #geom_text(aes(label = Code)) +
#   geom_point() +
#   theme_bw() +
#   scale_color_viridis_c(na.value = "lightgrey")

##### Quadrats Plots ####
mcl.df.q.forb <- filter(mcl.df.q, FunGroup == "Native Forb" | FunGroup == "Exotic Forb", iter < 150, n.plots >= 80)
mcl.df.q.forb <- merge(mcl.df.q.forb, trait.pca.forb[,-6], by.x = "Species_Name", by.y = "Species", all.x = T)

###### Average presence ####
# no relationship between avg presence and c or s estimate
ggplot(mcl.df.q.forb, aes(x = s, y = c, col = avg.P.sp)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey") # no relationship

# no relationship between avg presence and PC scores
ggplot(mcl.df.q.forb, aes(x = avg.P.sp, y = PC1)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = avg.P.sp, y = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = avg.P.sp, y = PC3)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

###### PCs ####
ggplot(mcl.df.q.forb, aes(x = PC1, y = c, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = PC1, y = s, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 


ggplot(mcl.df.q.forb, aes(x = PC2, y = c, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = PC2, y = s, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

# PC 3 isn't important
ggplot(mcl.df.q.forb, aes(x = PC3, y = c, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = PC3, y = s, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

## PC1 ##
ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC1)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC1)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~FunGroup) +
  scale_color_viridis_c(na.value = "lightgrey")

## PC2 ## shape, nitrogen
ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC2)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~FunGroup) +
  scale_color_viridis_c(na.value = "lightgrey")

## PC3 ##
ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC3)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC3)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~FunGroup) +
  scale_color_viridis_c(na.value = "lightgrey")

## PC4 ##
ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC4)) +
  geom_smooth(method = "lm", col = "black") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() +
  scale_color_viridis_c(na.value = "lightgrey")

ggplot(mcl.df.q.forb, aes(x = s, y = c, col = PC4)) +
  geom_smooth(method = "lm", col = "black") + 
  geom_text(aes(label = Code)) +
  #geom_point() +
  theme_bw() +
  facet_wrap(~FunGroup) +
  scale_color_viridis_c(na.value = "lightgrey")

###### Traits ####
# SHAPE
ggplot(mcl.df.q.forb, aes(x = shape, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw()

ggplot(mcl.df.q.forb, aes(x = shape, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 


# SIZE
ggplot(mcl.df.q.forb, aes(x = size.mm, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw()

ggplot(mcl.df.q.forb, aes(x = size.mm, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw()


# EMBRYO:SEED
ggplot(mcl.df.q.forb, aes(x = E.S, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  #geom_text(aes(label = Code)) +
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = E.S, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

# coat perm
ggplot(mcl.df.q.forb, aes(x = coat.perm.perc, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = coat.perm.perc, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

# coat thick size ratio

ggplot(mcl.df.q.forb, aes(x = coat.thick.per.size, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = coat.thick.per.size, y = s)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = both.thick, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 


# coat thick mass ratio

ggplot(mcl.df.q.forb, aes(x = coat.thick.per.mass, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 


# prop C
ggplot(mcl.df.q.forb, aes(x = prop.C, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = prop.C, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

# prop N
ggplot(mcl.df.q.forb, aes(x = prop.N, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = prop.N, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

#Height
ggplot(mcl.df.q.forb, aes(x = height.cm, y = c, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = height.cm, y = s, col = FunGroup, group = FunGroup)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

##### Forb Models ####
## by PC ##
summary(lm(s ~ PC1, data = mcl.df.q.forb))
summary(lm(c ~ PC1, data = mcl.df.q.forb))

summary(lm(s ~ PC2, data = mcl.df.q.forb))

summary(lm(c ~ PC2, data = mcl.df.q.forb))

# summary(lm(s ~ PC2, data = mcl.df.q.forb[mcl.df.q.forb$FunGroup == "Exotic Forb",]))

# summary(lm(c ~ PC2, data = mcl.df.q.forb))

summary(lm(s ~ PC3, data = mcl.df.q.forb))
summary(lm(c ~ PC3, data = mcl.df.q.forb))

summary(lm(s ~ PC4, data = mcl.df.q.forb))
summary(lm(c ~ PC4, data = mcl.df.q.forb))

## by traits ##
summary(lm(s ~ size.mm, data = mcl.df.q.forb))

summary(lm(c ~ size.mm, data = mcl.df.q.forb[mcl.df.q.forb$FunGroup == "Exotic Forb",])) # SIZE, Prop.N, drives colonization among exotic forbs

summary(lm(s ~ shape, data = mcl.df.q.forb))

summary(lm(c ~ shape*nat.inv, data = mcl.df.q.forb))

m1 <- lm(c ~ ldd2, data = mcl.df.q.forb[!is.na(mcl.df.q.forb$ldd2) & mcl.df.q.forb$ldd2 != "unknown",])
anova(m1)
pairs(emmeans(m1, ~ ldd2), adjust = "BH")

m1 <- lm(c ~ appendage.type, data = mcl.df.q.forb[!is.na(mcl.df.q.forb$appendage.type),])
anova(m1)

m1 <- lm(c ~ appendage.type, data = mcl.df.q.forb[!is.na(mcl.df.q.forb$appendage.type),])
anova(m1)

m1 <- lm(c ~ appendage, data = mcl.df.q.forb[!is.na(mcl.df.q.forb$appendage),])
anova(m1)

##### Quadrat Pairs ####
# survival seems to be linked to a high number of traits
# colonization, not much MAYBE size

# survival
# - shape
# - size
# - propC
# - C:N
# - height
# - E:S
# + coat thick per size

# colonization 
#+size
#- seed coat perm


# germ
# - seed coat perm

mcl.df.q.c.forb <- mcl.df.q.forb[mcl.df.q.forb$iter < 150 & mcl.df.q.forb$n.plots >= 80,-c(24,25,31)]
rm <- c("ACMBRA", "CERGLO", "LEPBIC", "CROSET", "GALAPA" ,"EROBRA")
mcl.df.q.c.forb <- filter(mcl.df.q.c.forb, !Code %in% rm)

# exotic forbs
## c: size, E.S
## s: N leads s

# native forbs
## c: size
## s: shape, size

# together
## c: NA
## s: shape, size, coat thick per size 
pairs(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native", c(4,5,19:22, 25:29)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "native", c(4,5,31:33, 37:40)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)


pairs(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive", c(4,5,19:22, 25:29)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

pairs(mcl.df.q.forb.PC[mcl.df.q.forb.PC$nat.inv == "invasive", c(4,5,31:33, 37:40)], 
      lower.panel=panel.lmline, upper.panel=panel.cor)

##### Transect Pairs ####
# 
# # survival
# ## -E:S
# # colonization 
# ## -CN
# ## -coat perm
# ## +PC3
# 
# mcl.df.t.c.forb <- mcl.df.t.forb[mcl.df.t.forb$iter < 150 & mcl.df.t.forb$n.plots > 8,]
# mcl.df.t.c.forb <- mcl.df.t.c.forb[complete.cases(mcl.df.t.c.forb),]
# 
# pairs(mcl.df.t.c.forb[, c(4:5,25:30)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)
# 
# pairs(mcl.df.t.c.forb[, c(4:5,31:35)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)
# 
# pairs(mcl.df.t.c.forb[, c(4:5,36:39)], 
#       lower.panel=panel.lmline, upper.panel=panel.cor)

### S/C comparison ####
mcl.df.q.all <- filter(mcl.df.q.all, iter < 150, n.plots >= 80, FunGroup != "Native Grass")

mcl.df.q.all$FunGroup <- factor(mcl.df.q.all$FunGroup, levels = c("Exotic Grass", "Native Forb", "Exotic Forb"))

plot.c <- ggplot(mcl.df.q.all, aes(x = FunGroup, y = c)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.title.y = element_text(size = 13)
  ) +
  ylab("Probability of spatial dispersal")

plot.s <- ggplot(mcl.df.q.all, aes(x = FunGroup, y = s)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.title.y = element_text(size = 13)
  ) +
  ylab("Probability of temporal dispersal")

c.s.fig <- ggarrange(plot.c, plot.s, labels = c("(a)", "(b)"), ncol = 2)

ggsave("McLaughlin/Figures/c.s.fig.jpeg", c.s.fig, height = 4, width = 9, units = "in", dpi = 600)

mcl.df.q.forb$ldd2 <- factor(mcl.df.q.forb$ldd2, levels = c("adhesive", "ingestion", "wind", "ant", "unassisted"))


ggplot(mcl.df.q.forb[!is.na(mcl.df.q.forb$ldd2) & mcl.df.q.forb$ldd2 != "unknown",], aes(x = ldd2, y = c)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.title.y = element_text(size = 13)
  ) +
  ylab("Probability of spatial dispersal")


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


# Hierarchical Partitioning ####

## Output
#the relative importance of any individual predictor can be simply estimated as its unique contribution to the total model R2 plus its average shared contributions with the other predictors.

# Unique: part of variation ONLY explained by that variable 

# Average.share: total shared/number of predictors sharing; Shared fractions represent the variation in the response variable explained by the correlation of the predictors involved (i.e. correlation in the predictive space). Larger shared fractions imply that more multicollinearity is present in the model. #Negative common (shared) variation is possible when predictors act as suppressors of other predictors

# Individual = Relative importance of individual predictor = unique + shared portion where shared equals = total shared/number of predictors sharing; sum of this adds up to R2

# I.perc = adds up up to 100%

# I *THINK* RDA works best with normalized data
#https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html

##### Prep ####
rm <- c("ACMBRA", "CERGLO", "LEPBIC", "CROSET", "GALAPA" ,"EROBRA", "GASPHL")
# these species are removed because I have c/s estimates for them but no seed traits


# should only the traits that can be directly related to survival go in? like height and settling time can only really be indirectly related

trait <- c("chem.mass.mg",  "shape", "size.mm", "prop.C", "prop.N", "set.time.mpsec", "ldd2", "coat.thick.per.size", "height.cm", "coat.perm.perc")#, "coat.thick.per.size", "coat.thick.per.mass", "ldd")


survival <- c("chem.mass.mg",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.size")

colonization <- c("chem.mass.mg",  "shape", "set.time.mpsec", "height.cm", "ldd2", "size.mm")

# pick traits without correlations but RDA might be able to handle them together

# how is this different from a hierarchical model? order matters
# conceptually it's fairly similar 

mcl.df.q.forb <- filter(mcl.df.q.forb, !Code %in% rm)


##### All Forbs ####

###### S: All Forbs ####
# spe <- mcl.df.q.forb$s
# spe <- mcl.df.q.forb$c
# env <- mcl.df.q.forb[,colnames(mcl.df.q.forb) %in% trait]
# 
# spe.rda <- rda(spe~.,env)
# vif.cca(spe.rda)
# anova(spe.rda, by = "margin")

HPA.s <- rdacca.hp(
  dv = mcl.df.q.forb$s,
  iv = mcl.df.q.forb[,colnames(mcl.df.q.forb) %in% survival],
  method = "RDA",
  type = "adjR2",
  var.part = T
)



HPA.s$Hier.part
#HPA.s$Var.part
HPA.s$Total_explained_variation #19.7%

plot.rdaccahp(HPA.s) # size and shape

# all traits
HPA.s <- rdacca.hp(
  dv = mcl.df.q.forb$s,
  iv = mcl.df.q.forb[ ,colnames(mcl.df.q.forb) %in% trait],
  method = "RDA",
  type = "adjR2"
)

HPA.s$Hier.part

HPA.s$Total_explained_variation # total explained variation goes from 60% to 11% when including invasives

plot.rdaccahp(HPA.s)

###### C: All Forbs ####

HPA.c <- rdacca.hp(
  dv = mcl.df.q.forb$c,
  iv = mcl.df.q.forb[ ,colnames(mcl.df.q.forb) %in% colonization],
  method = "RDA",
  type = "adjR2"
)

HPA.c$Hier.part

HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c) # dispersal mechanism and settling time


# trying all traits together
# colonization
HPA.c <- rdacca.hp(
  dv = mcl.df.q.forb$c,
  iv = mcl.df.q.forb[ ,colnames(mcl.df.q.forb) %in% trait],
  method = "RDA",
  type = "adjR2"
)

HPA.c$Hier.part
HPA.c$Total_explained_variation

plot.rdaccahp(HPA.c)


##### Natives #####

###### S: Natives #####

test <- rdacca.hp(
  dv = mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native",]$s,
  iv = mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native", colnames(mcl.df.q.forb) %in% survival],
  
)

test$Total_explained_variation #22.7%

test$Hier.part
#test$Var.part

plot.rdaccahp(test) # shape, size, prop C

###### C: Natives #####
test <- rdacca.hp(
  dv = mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native",]$c,
  iv = mcl.df.q.forb[ mcl.df.q.forb$nat.inv == "native",colnames(mcl.df.q.forb) %in% colonization],
)

test$Total_explained_variation #14%
test$Hier.part
plot.rdaccahp(test) # dispersal mechanism and mass



##### Invasives #####

###### S: Invasives #####
test <- rdacca.hp(
  dv = mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive",]$s,
  iv = mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive",colnames(mcl.df.q.forb) %in% survival],
)

test$Total_explained_variation
test$Hier.part
plot.rdaccahp(test) # coat thickness, prop N, shape

###### C: Invasives #####

test <- rdacca.hp(
  dv = mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive",]$c,
  iv = mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive", colnames(mcl.df.q.forb) %in% colonization],
)

test$Total_explained_variation

plot.rdaccahp(test) # shape and height

## Panel Plot #####
# natives: 
## S: shape, size
## C: long distance dispersal, size

# invasives
## S: seed coat, prop N
## C: shape, height 

## TEMPORAL ##

ggplot(mcl.df.q.forb, aes(x = coat.thick.per.size, y = s, col = nat.inv)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  theme_bw() +
  labs(x = "Coat thickness to size ratio", 
       y = "Probability of Temporal Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

ggplot(mcl.df.q.forb, aes(x = prop.N, y = s, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  labs(x = "Proportion Nitrogen", 
       y = "Probability of Temporal Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

# amsinckia seems suspiciously low here, not sure what's going on but the other two amsinckia species have twice as much carbon as this one.
ggplot(mcl.df.q.forb[mcl.df.q.forb$Species_Name != "Amsinckia menziesii",], aes(x = prop.C, y = s, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  theme_bw() +
  geom_point() +
  labs(x = "Proportion Carbon", 
       y = "Probability of Temporal Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

ggplot(mcl.df.q.forb, aes(x = shape, y = s, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  theme_bw() +
  geom_point() +
  labs(x = "Seed Shape", 
       y = "Probability of Temporal Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

ggplot(mcl.df.q.forb, aes(x = size.mm, y = s, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  theme_bw() +
  geom_point() +
  labs(x = "Seed Size (mm)", 
       y = "Probability of Temporal Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

## SPATIAL ##
mcl.df.q.forb.sum <- mcl.df.q.forb %>%
  dplyr::group_by(ldd2) %>%
  dplyr::summarize(c.mean = mean(c, na.rm = T),
                   c.se = calcSE(c))

ggplot(mcl.df.q.forb[mcl.df.q.forb$ldd != 0,], aes(x = ldd, y = c)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() 

ggplot(mcl.df.q.forb, aes(x = shape, y = c, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  labs(x = "Shape", 
       y = "Probability of Spatial Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

ggplot(mcl.df.q.forb, aes(x = size.mm, y = c, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  labs(x = "Size (mm)", 
       y = "Probability of Spatial Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

ggplot(mcl.df.q.forb, aes(x = height.cm, y = c, col = nat.inv)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_bw() +
  labs(x = "Height.cm", 
       y = "Probability of Spatial Dispersal") +
  scale_color_manual(values = c("#238A8DFF", "#482677FF"))

#### Spider Plots ####
# Library



# Library
library(fmsb)
 
 
# Create data: note in High school for several students
set.seed(99)
data <- as.data.frame(matrix( sample( 0:20 , 15 , replace=F) , ncol=5))
colnames(data) <- c("math" , "english" , "biology" , "music" , "R-coding" )
rownames(data) <- paste("mister" , letters[1:3] , sep="-")
 
# To use the fmsb package, I have to add 2 lines to the dataframe: the max and min of each variable to show on the plot!
data <- rbind(rep(20,5) , rep(0,5) , data)
 
# plot with default options:
radarchart(data)

## Colonization 
colonization <- c("c", "chem.mass.mg",  "shape", "set.time.mpsec", "height.cm", "size.mm")

spider.c.n <- cor(mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native", colnames(mcl.df.q.forb) %in% colonization])

spider.c.i <- cor(mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive", colnames(mcl.df.q.forb) %in% colonization])

## Survival
survival <- c("s", "chem.mass.mg",  "shape", "size.mm", "prop.C", "prop.N", "coat.perm.perc", "coat.thick.per.size")

spider.s.n <- cor(mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native", colnames(mcl.df.q.forb) %in% survival])

spider.s.i <- cor(mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive", colnames(mcl.df.q.forb) %in% survival])

trait.p <- data.frame()

for(i in survival[-1]){
  tmp <- cor.test(mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native",]$s, mcl.df.q.forb[mcl.df.q.forb$nat.inv == "native",i])
  trait.p <- rbind(trait.p, c(trait = i, p.value = tmp$p.value))
}

colnames(trait.p) <- c("trait", "p.value")

trait.p <- data.frame()

for(i in survival[-1]){
  tmp <- cor.test(mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive",]$s, mcl.df.q.forb[mcl.df.q.forb$nat.inv == "invasive",i])
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
# 
# spider.s <- rbind(spider.s.n, spider.s.i)
# spider.s <- spider.s[,c(9,2:8)]
# 
# 
# rownames(spider.s) <- spider.s[,1]
# 
# spider.s <- rbind(max = c(1,1,1,1,1,1,1), min = c(-1,-1,-1,-1,-1,-1,-1), spider.s[,-1])
# 
# radarchart(spider.s,
#            axistype = 1,
#            caxislabels = seq(-1, 1, 0.5))
# 
# 
# 
# spider.c.i <- as.data.frame(spider.c.i)
# spider.c.i$nat.inv <- "invasive"
# 
# spider.c.n <- as.data.frame(spider.c.n)
# spider.c.n$nat.inv <- "native"
# 
# spider.c <- rbind(spider.c.n[1,], spider.c.i[1,])
# spider.c <- spider.c[,c(7,2:6)]
# 
# 
# rownames(spider.c) <- spider.c[,1]
# 
# spider.c <- rbind(max = c(1,1,1,1,1,1), min = c(-1,-1,-1,-1,-1,-1), spider.c[,-1])
# 
# radarchart(spider.c,
#            axistype = 1,
#            caxislabels = seq(-1, 1, 0.5))
####
spider.s.i <- as.data.frame(spider.s.i)
spider.s.i$nat.inv <- "invasive"
spider.s.i$trait <- row.names(spider.s.i)

spider.s.n <- as.data.frame(spider.s.n)
spider.s.n$nat.inv <- "native"
spider.s.n$trait <- row.names(spider.s.n)

cor.s <- as.data.frame(rbind(spider.s.i[-1,c(1,9,10)], spider.s.n[-1,c(1,9,10)]))
colnames(cor.s)[1] <- "s.cor"


##### Normalize ####
# 
# hist(trait.pca.forb$ldd)
# 
# hist(log(trait.pca.forb$wing.loading))
# trait.pca.forb$wing.loading <- log(trait.pca.forb$wing.loading)
# 
# hist(log(trait.pca.forb$cn))
# trait.pca.forb$cn <- log(trait.pca.forb$cn)
# 
# hist(trait.pca.forb$prop.C)
# 
# hist(log(trait.pca.forb$prop.N))
# trait.pca.forb$prop.N <- log(trait.pca.forb$prop.N)
# 
# hist(log(trait.pca.forb$coat.perm.perc))
# trait.pca.forb$coat.perm.perc <- log(trait.pca.forb$coat.perm.perc)
# 
# hist(log(trait.pca.forb$coat.thick.per.size))
# trait.pca.forb$coat.thick.per.size <- log(trait.pca.forb$coat.thick.per.size)
# 
# hist(log(trait.pca.forb$coat.thick.per.mass))
# trait.pca.forb$coat.thick.per.mass <- log(trait.pca.forb$coat.thick.per.mass)
# 
# hist(log(trait.pca.forb$morph.mass.mg))
# trait.pca.forb$morph.mass.mg <- log(trait.pca.forb$morph.mass.mg)
# 
# hist(log(trait.pca.forb$chem.mass.mg))
# trait.pca.forb$chem.mass.mg <- log(trait.pca.forb$chem.mass.mg)
# 
# hist(log(trait.pca.forb$size.mm))
# trait.pca.forb$size.mm <- log(trait.pca.forb$size.mm)
# 
# hist(trait.pca.forb$set.time.mpsec)
# 
# hist(log(trait.pca.forb$E.S))
# trait.pca.forb$E.S <- log(trait.pca.forb$E.S)
# 
# hist(log(trait.pca.forb$both.thick))
# trait.pca.forb$both.thick <- log(trait.pca.forb$both.thick)
